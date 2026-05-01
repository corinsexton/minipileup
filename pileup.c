// This piece of code is modified from samtools/bam2depth.c
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "ksort.h"
#include "ketopt.h"

#define VERSION "1.4-r20"

void *bed_read(const char *fn);
void *bed_read_vcf(const char *fn);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void bed_destroy(void *_h);
hts_reglist_t *bed_to_reglist(const void *_h, sam_hdr_t *hdr, int *nregs);
typedef struct { int tid; int beg, end; } bed_site_t;
bed_site_t *bed_get_sorted_sites(const void *_h, sam_hdr_t *hdr, int *n_out);

typedef struct {     // auxiliary data structure
	htsFile *fp;     // the file handler (BAM, CRAM, or SAM)
	hts_itr_t *itr;  // NULL if a region not specified
	sam_hdr_t *h;
	int min_mapQ, min_len; // mapQ filter; length filter
	int min_supp_len;
	int proper_only;
	int indexed_bed; // 1 when itr was built from bed regions via sam_itr_regions
	void *bed;       // bedidx if not NULL
} aux_t;

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret = aux->itr? bam_itr_next(aux->fp, aux->itr, b) : sam_read1(aux->fp, aux->h, b);
	if (ret < 0) return ret;
	if (b->core.tid < 0) b->core.flag |= BAM_FUNMAP;
	// htslib >=1.0 no longer auto-filters these in bam_plp_push (old BAM_PLP_MASK); do it here
	if (b->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP))
		b->core.flag |= BAM_FUNMAP;
	if (aux->proper_only && (b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR))
		b->core.flag |= BAM_FUNMAP;
	if (!(b->core.flag&BAM_FUNMAP)) {
		if ((int)b->core.qual < aux->min_mapQ) {
			b->core.flag |= BAM_FUNMAP;
		} else if (aux->min_len > 0 || aux->min_supp_len > 0 || aux->bed) {
			int k, qlen = 0, tlen = 0;
			const char *chr = aux->h->target_name[b->core.tid];
			const uint32_t *cigar = bam_get_cigar(b);
			for (k = 0; k < b->core.n_cigar; ++k) { // compute the query length in the alignment
				int op = bam_cigar_op(cigar[k]);
				int oplen = bam_cigar_oplen(cigar[k]);
				if ((bam_cigar_type(op)&1) && op != BAM_CSOFT_CLIP)
					qlen += oplen;
				if (bam_cigar_type(op)&2)
					tlen += oplen;
			}
			if (qlen < aux->min_len) b->core.flag |= BAM_FUNMAP;
			if (qlen < aux->min_supp_len && (b->core.flag&BAM_FSUPPLEMENTARY)) b->core.flag |= BAM_FUNMAP;
			if (aux->bed && !aux->indexed_bed && !(b->core.flag&BAM_FUNMAP) && !bed_overlap(aux->bed, chr, b->core.pos, b->core.pos + tlen))
				b->core.flag |= BAM_FUNMAP;
		}
	}
	return ret;
}

typedef struct {
	uint32_t is_skip:1, is_rev:1, b:4, q:8, is_del:1, k:17; // b=base, q=quality, k=allele id
	int indel; // <0: deletion; >0: insertion
	uint64_t hash;
	uint64_t pos; // i<<32|j: j-th read of the i-th sample
} allele_t;

#define allele_lt(a, b) ((a).hash < (b).hash || ((a).hash == (b).hash && (a).indel < (b).indel))
KSORT_INIT(allele, allele_t, allele_lt)

#define uint64_lt(a, b) ((a) < (b))
KSORT_INIT(uint64, uint64_t, uint64_lt)

static inline int is_proper_same_chrom(const bam_pileup1_t *p)
{
	uint32_t f = p->b->core.flag;
	return (f & BAM_FPAIRED) && (f & BAM_FPROPER_PAIR) &&
	       !(f & BAM_FMUNMAP) && p->b->core.tid == p->b->core.mtid;
}

static inline uint64_t pair_fingerprint(const bam_pileup1_t *p)
{
	int32_t a = p->b->core.pos, b = p->b->core.mpos;
	return a <= b ? (uint64_t)(uint32_t)a << 32 | (uint32_t)b
	              : (uint64_t)(uint32_t)b << 32 | (uint32_t)a;
}

// fp_used encodes two bits per slot:
//   bit 0 (value 1): skip decision — 0=keep (first of pair), 1=skip (second of pair)
//   bit 1 (value 2): consumed — set when the slot is assigned to a read in the main loop
//
// Returns 1 (skip this read) or 0 (count this read).
// Finds the first unconsumed slot for fp and returns its skip decision.
static inline int fp_seen_before(uint64_t *fp_buf, uint8_t *fp_used, int n_fp, uint64_t fp)
{
	int lo = 0, hi = n_fp;
	while (lo < hi) {
		int mid = (lo + hi) >> 1;
		if (fp_buf[mid] < fp) lo = mid + 1;
		else hi = mid;
	}
	if (lo >= n_fp || fp_buf[lo] != fp) return 0;
	// advance past already-consumed slots for this fingerprint
	while (lo < n_fp && fp_buf[lo] == fp && (fp_used[lo] & 2)) ++lo;
	if (lo >= n_fp || fp_buf[lo] != fp) return 0;
	int skip = fp_used[lo] & 1; // pre-computed skip decision for this slot
	fp_used[lo] |= 2;           // mark consumed
	return skip;
}

static inline allele_t pileup2allele(const bam_pileup1_t *p, int min_baseQ, uint64_t pos, int ref, int trim_len, int del_as_allele)
{ // collect allele information given a pileup1 record
	allele_t a;
	int i;
	const uint8_t *seq = bam_get_seq(p->b);
	a.k = (1<<17) - 1; // this will be set in count_alleles()
	a.q = bam_get_qual(p->b)[p->qpos];
	a.is_rev = bam_is_rev(p->b);
	if (del_as_allele) {
		a.is_skip = (p->is_refskip || a.q < min_baseQ);
		a.is_del = p->is_del;
	} else {
		a.is_skip = (p->is_del || p->is_refskip || a.q < min_baseQ);
		a.is_del = 0;
	}
	if (p->qpos < trim_len || p->b->core.l_qseq - p->qpos < trim_len) a.is_skip = 1;
	a.indel = p->indel;
	a.b = a.hash = bam_seqi(seq, p->qpos);
	a.pos = pos;
	if (a.is_del) a.hash = (uint64_t)-1;
	else if (p->indel > 0) // compute the hash for the insertion
		for (i = 0; i < p->indel; ++i)
			a.hash = (a.hash<<4) + a.hash + bam_seqi(seq, p->qpos + i + 1);
	a.hash = a.hash << 1 >> 1;
	if (p->indel != 0 || a.b != ref || ref == 15 || a.is_del) // the highest bit tells whether it is a reference allele or not
		a.hash |= 1ULL<<63;
	return a;
}

static inline void print_allele(FILE *out, const bam_pileup1_t *p, int l_ref, const char *ref, int pos, int max_del, int is_vcf, int del_as_allele)
{ // print the allele. The format depends on is_vcf.
	const uint8_t *seq = bam_get_seq(p->b);
	int i, rest = max_del;
	if (del_as_allele && p->is_del) {
		fputc('*', out);
		return;
	}
	fputc(seq_nt16_str[bam_seqi(seq, p->qpos)], out);
	if (p->indel > 0) {
		if (!is_vcf) fprintf(out, "+%d", p->indel);
		for (i = 1; i <= p->indel; ++i)
			fputc(seq_nt16_str[bam_seqi(seq, p->qpos + i)], out);
	} else if (p->indel < 0) {
		if (!is_vcf) {
			fprintf(out, "%d", p->indel);
			for (i = 1; i <= -p->indel; ++i)
				fputc(pos + i < l_ref? toupper(ref[pos+i]) : 'N', out);
		} else rest -= -p->indel, pos += -p->indel;
	}
	if (is_vcf)
		for (i = 1; i <= rest; ++i)
			fputc(pos + i < l_ref? toupper(ref[pos+i]) : 'N', out);
}

typedef struct {
	int n_a, n_alleles, max_del; // n_a: #reads used to compute quality sum; max_del: max deletion length
	int tot_dp, max_dp, n_cnt, max_cnt;
	allele_t *a; // allele of each read, of size $n_a
	int *cnt_strand, *cnt_supp; // cnt_strand: count of supporting reads on both strands; cnt_supp: sum of both strands
	int *support, *support_strand; // support across entire $a. It points to the last "row" of cnt_supp/cnt_strand
	int len, max_len;
	char *seq;
	uint64_t *fp_buf;  // pair fingerprint buffer for overlap deduplication
	uint8_t  *fp_used; // which fingerprint slots have been claimed this column
	int max_fp;
} paux_t;

static void count_alleles(paux_t *pa, int n)
{
	allele_t *a = pa->a;
	int i;
	a[0].k = 0; // the first allele is given allele id 0
	pa->max_del = a[0].indel < 0? -a[0].indel : 0;
	for (i = pa->n_alleles = 1; i < pa->n_a; ++i) {
		if (a[i].indel != a[i-1].indel || a[i].hash != a[i-1].hash) // change of allele
			++pa->n_alleles;
		a[i].k = pa->n_alleles - 1;
		pa->max_del = pa->max_del > -a[i].indel? pa->max_del : -a[i].indel; // max deletion
	}
	// collect per-BAM counts
	pa->n_cnt = pa->n_alleles * (n + 1);
	if (pa->n_cnt > pa->max_cnt) { // expand the arrays if necessary
		pa->max_cnt = pa->n_cnt;
		kroundup32(pa->max_cnt);
		pa->cnt_strand = (int*)realloc(pa->cnt_strand, pa->max_cnt * 2 * sizeof(int));
		pa->cnt_supp = (int*)realloc(pa->cnt_supp, pa->max_cnt * sizeof(int));
	}
	memset(pa->cnt_strand, 0, pa->n_cnt * 2 * sizeof(int));
	pa->support_strand = pa->cnt_strand + pa->n_alleles * n * 2;
	memset(pa->cnt_supp, 0, pa->n_cnt * sizeof(int));
	pa->support = pa->cnt_supp + pa->n_alleles * n; // points to the last row of cnt_supp
	for (i = 0; i < pa->n_a; ++i) { // compute counts and sums of qualities
		int j = (a[i].pos>>32)*pa->n_alleles + a[i].k;
		pa->cnt_strand[j<<1|a[i].is_rev]++;
		pa->cnt_supp[j]++;
		pa->support[a[i].k]++;
		pa->support_strand[a[i].k<<1|a[i].is_rev]++;
	}
}

// Arguments passed to each worker thread (or used directly for single-threaded run).
typedef struct {
	// files
	int n_files;
	char **filenames;     // argv + o.ind; read-only
	const char *fname_ref;
	// shared read-only state
	sam_hdr_t *h;
	void *bed;
	bed_site_t *sites;
	int nsites;
	// site range this thread is responsible for: [si_beg, si_end)
	int si_beg, si_end;
	// for the nsites==0 (streaming) path
	int tid, beg, end;
	// filter params
	int baseQ, mapQ, min_len, min_supp_len, proper_only;
	int is_vcf, var_only, show_2strand, trim_len, del_as_allele;
	int min_support, min_support_strand;
	int dedup_overlap; // -D: deduplicate overlapping paired-end reads
	double min_af;
	// output: written here by the thread; main flushes to stdout in order
	FILE *out;
	char *out_buf;
	size_t out_len;
} thread_arg_t;

static void process_sites(thread_arg_t *a)
{
	int i, j, n = a->n_files;
	int tid = -1, pos;
	int last_tid = -1;
	int l_ref = 0;
	char *ref = NULL;
	bam_mplp_t mplp;
	FILE *out = a->out;

	int *n_plp = (int*)calloc(n, sizeof(int));
	const bam_pileup1_t **plp = (const bam_pileup1_t**)calloc(n, sizeof(const bam_pileup1_t*));
	paux_t pa;
	memset(&pa, 0, sizeof(paux_t));

	// each thread opens its own fai to avoid sharing a seekable file handle
	faidx_t *fai = a->fname_ref? fai_load(a->fname_ref) : NULL;

	// open per-thread file handles; use shared indices if pre-loaded
	aux_t **data = (aux_t**)calloc(n, sizeof(aux_t*));
	hts_idx_t **idx = (hts_idx_t**)calloc(n, sizeof(hts_idx_t*));
	for (i = 0; i < n; ++i) {
		data[i] = (aux_t*)calloc(1, sizeof(aux_t));
		data[i]->fp = hts_open(a->filenames[i], "r");
		if (a->fname_ref) hts_set_fai_filename(data[i]->fp, a->fname_ref);
		data[i]->min_mapQ = a->mapQ;
		data[i]->min_len  = a->min_len;
		data[i]->min_supp_len = a->min_supp_len;
		data[i]->proper_only = a->proper_only;
		data[i]->bed = a->bed;
		sam_hdr_t *htmp = sam_hdr_read(data[i]->fp); // must read to advance past header
		sam_hdr_destroy(htmp);
		data[i]->h = a->h; // use shared header (read-only)
		if (a->nsites > 0) {
			idx[i] = sam_index_load(data[i]->fp, a->filenames[i]);
		} else if (a->tid >= 0) { // -r region, no bed: set up region iterator
			hts_idx_t *tmp = sam_index_load(data[i]->fp, a->filenames[i]);
			if (tmp) {
				data[i]->itr = bam_itr_queryi(tmp, a->tid, a->beg, a->end);
				hts_idx_destroy(tmp);
			}
		}
	}

	int site_beg = a->beg, site_end = a->end;
	int niter = a->nsites > 0? a->si_end : a->si_beg + 1; // si_beg+1 for the nsites==0 single pass
	for (int si = a->si_beg; si < niter; ++si) {
		if (a->nsites > 0) {
			int site_tid_val = a->sites[si].tid;
			site_beg = a->sites[si].beg;
			site_end = a->sites[si].end;
			if (last_tid != site_tid_val) {
				if (fai) { free(ref); ref = fai_fetch(fai, a->h->target_name[site_tid_val], &l_ref); }
				last_tid = site_tid_val;
				pa.len = 0;
			}
			for (i = 0; i < n; ++i) {
				if (data[i]->itr) { bam_itr_destroy(data[i]->itr); data[i]->itr = NULL; }
				if (idx[i]) data[i]->itr = bam_itr_queryi(idx[i], site_tid_val, site_beg, site_end);
			}
		}
		mplp = bam_mplp_init(n, read_bam, (void**)data);
		while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) {
			if (pos < site_beg || pos >= site_end) continue;
			if (a->bed && !bed_overlap(a->bed, a->h->target_name[tid], pos, pos + 1)) continue;
			for (i = pa.tot_dp = 0; i < n; ++i) pa.tot_dp += n_plp[i];
			if (last_tid != tid) {
				if (fai) { // switch of chromosomes
					free(ref);
					ref = fai_fetch(fai, a->h->target_name[tid], &l_ref);
				}
				last_tid = tid; pa.len = 0;
			}
			if (pa.tot_dp) {
			int k, r = 15, shift = 0, qual, n_fp = 0;
			allele_t *al;
			if (pa.tot_dp + 1 > pa.max_dp) { // expand array
				pa.max_dp = pa.tot_dp + 1;
				kroundup32(pa.max_dp);
				pa.a = (allele_t*)realloc(pa.a, pa.max_dp * sizeof(allele_t));
			}
			if (a->dedup_overlap && pa.tot_dp + 1 > pa.max_fp) {
				pa.max_fp = pa.tot_dp + 1;
				kroundup32(pa.max_fp);
				pa.fp_buf  = (uint64_t*)realloc(pa.fp_buf,  pa.max_fp * sizeof(uint64_t));
				pa.fp_used = (uint8_t*)realloc(pa.fp_used,  pa.max_fp * sizeof(uint8_t));
			}
			al = pa.a;
			// collect alleles; ref is fetched per full chromosome so offset is just pos
			r = (ref && pos < l_ref)? seq_nt16_table[(int)ref[pos]] : 15;
			if (a->dedup_overlap) {
				// pre-pass: build sorted fingerprint table of all proper paired reads in this column
				for (i = 0; i < n; ++i)
					for (j = 0; j < n_plp[i]; ++j)
						if (is_proper_same_chrom(&plp[i][j]))
							pa.fp_buf[n_fp++] = pair_fingerprint(&plp[i][j]);
				ks_introsort(uint64, n_fp, pa.fp_buf);
				// assign alternating skip decisions within each run of identical fingerprints:
				// even positions (0, 2, 4, ...) = keep (bit 0 = 0)
				// odd positions  (1, 3, 5, ...) = skip (bit 0 = 1)
				uint64_t prev_fp = (uint64_t)-1;
				int fp_parity = 0;
				for (int fi = 0; fi < n_fp; ++fi) {
					if (pa.fp_buf[fi] != prev_fp) { prev_fp = pa.fp_buf[fi]; fp_parity = 0; }
					pa.fp_used[fi] = fp_parity++ & 1; // bit 1 (consumed) starts at 0
				}
			}
			for (i = pa.n_a = 0; i < n; ++i)
				for (j = 0; j < n_plp[i]; ++j) {
					al[pa.n_a] = pileup2allele(&plp[i][j], a->baseQ, (uint64_t)i<<32 | j, r, a->trim_len, a->del_as_allele);
					if (a->dedup_overlap && !al[pa.n_a].is_skip && is_proper_same_chrom(&plp[i][j]))
						if (fp_seen_before(pa.fp_buf, pa.fp_used, n_fp, pair_fingerprint(&plp[i][j])))
							al[pa.n_a].is_skip = 1;
					if (!al[pa.n_a].is_skip) ++pa.n_a;
				}
			if (pa.n_a == 0) continue; // no reads are good enough; zero effective coverage
			// count alleles
			ks_introsort(allele, pa.n_a, pa.a);
			count_alleles(&pa, n);
			// squeeze out weak alleles
			for (i = k = 0; i < pa.n_a; ++i)
				if (pa.support[al[i].k] >= a->min_support && pa.support[al[i].k] >= pa.n_a * a->min_af
					&& pa.support_strand[al[i].k<<1] >= a->min_support_strand && pa.support_strand[al[i].k<<1|1] >= a->min_support_strand)
				{
					al[k++] = al[i];
				}
			if (k < pa.n_a) {
				if (k == 0) continue; // no alleles are good enough
				pa.n_a = k;
				count_alleles(&pa, n);
			}

			if (a->var_only && pa.n_alleles == 1 && al[0].hash>>63 == 0) continue; // var_only mode, but no ALT allele; skip
			if (a->var_only && pa.n_alleles <= 2 && a->del_as_allele) {
				int n_ref = 0, n_del = 0;
				for (i = 0; i < pa.n_a; ++i)
					if (al[i].is_del) ++n_del;
					else if (al[i].hash>>63 == 0) ++n_ref;
				if (n_ref + n_del == pa.n_a) continue;
			}
			// print VCF or allele summary
			fputs(a->h->target_name[tid], out); fprintf(out, "\t%d", pos+1);
			if (a->is_vcf) {
				fputs("\t.\t", out);
				for (i = 0; i <= pa.max_del; ++i) // print the reference allele up to the longest deletion
					fputc(ref && pos + i < l_ref? ref[pos + i] : 'N', out);
				fputc('\t', out);
			} else fprintf(out, "\t%c\t", ref && pos < l_ref? ref[pos] : 'N');
			// print alleles
			if (!a->is_vcf || al[0].hash>>63) { // print if there is no reference allele
				print_allele(out, &plp[al[0].pos>>32][(uint32_t)al[0].pos], l_ref, ref, pos, pa.max_del, a->is_vcf, a->del_as_allele);
				if (pa.n_alleles > 1) fputc(',', out);
			}
			for (i = k = 1; i < pa.n_a; ++i)
				if (al[i].indel != al[i-1].indel || al[i].hash != al[i-1].hash) {
					print_allele(out, &plp[al[i].pos>>32][(uint32_t)al[i].pos], l_ref, ref, pos, pa.max_del, a->is_vcf, a->del_as_allele);
					if (++k != pa.n_alleles) fputc(',', out);
				}
			if (a->is_vcf && pa.n_alleles == 1 && al[0].hash>>63 == 0) fputc('.', out); // print placeholder if there is only the reference allele
			// compute and print qual
			for (i = !(al[0].hash>>63), qual = 0; i < pa.n_alleles; ++i)
				qual = qual > pa.support[i]? qual : pa.support[i];
			if (a->is_vcf) fprintf(out, "\t%d\t.\t.\tGT:%s", qual, a->show_2strand? "ADF:ADR" : "AD");
			// print counts
			shift = (a->is_vcf && al[0].hash>>63); // in VCF, if there is no ref allele, we need to shift the allele number
			for (i = k = 0; i < n; ++i, k += pa.n_alleles) {
				int max1 = 0, max2 = 0, a1 = -1, a2 = -1, *sum_q = &pa.cnt_supp[k];
				// estimate genotype
				for (j = 0; j < pa.n_alleles; ++j)
					if (sum_q[j] > max1) max2 = max1, a2 = a1, max1 = sum_q[j], a1 = j;
					else if (sum_q[j] > max2) max2 = sum_q[j], a2 = j;
				if (max1 == 0 || (a->min_support > 0 && max1 < a->min_support)) a1 = a2 = -1;
				else if (max2 == 0 || (a->min_support > 0 && max2 < a->min_support)) a2 = a1;
				// turn the genotype to homozygous if min_af is set and the minor allele does not have high enough frequency
				if (a->min_af > 0.0 && a->min_af < 0.5 && a1 >= 0 && a2 >= 0 && a1 != a2 && max2 < (max1 + max2) * a->min_af)
					a1 = a2;
				// print genotypes
				if (a1 < 0) fprintf(out, "\t./.:");
				else fprintf(out, "\t%d/%d:", a1 + shift, a2 + shift);
				// print counts
				if (a->show_2strand) {
					if (shift) fputs("0,", out);
					for (j = 0; j < pa.n_alleles; ++j) {
						if (j) fputc(',', out);
						fprintf(out, "%d", pa.cnt_strand[(k+j)<<1]);
					}
					fputc(':', out);
					if (shift) fputs("0,", out);
					for (j = 0; j < pa.n_alleles; ++j) {
						if (j) fputc(',', out);
						fprintf(out, "%d", pa.cnt_strand[(k+j)<<1|1]);
					}
				} else {
					if (shift) fputs("0,", out);
					for (j = 0; j < pa.n_alleles; ++j) {
						if (j) fputc(',', out);
						fprintf(out, "%d", pa.cnt_supp[k+j]);
					}
				}
			} // ~for(i)
			fputc('\n', out);
		} // ~if(pa.tot_dp)
		} // ~while()
		bam_mplp_destroy(mplp);
	} // ~for(si)

	// cleanup
	free(n_plp); free(plp);
	free(pa.cnt_strand); free(pa.cnt_supp); free(pa.a); free(pa.seq);
	free(pa.fp_buf); free(pa.fp_used);
	for (i = 0; i < n; ++i) {
		hts_close(data[i]->fp);
		if (data[i]->itr) bam_itr_destroy(data[i]->itr);
		if (idx[i]) hts_idx_destroy(idx[i]);
		free(data[i]);
	}
	free(data); free(idx);
	if (ref) free(ref);
	if (fai) fai_destroy(fai);
}

static void *worker(void *arg) { process_sites((thread_arg_t*)arg); return NULL; }

int main(int argc, char *argv[])
{
	int i, n, tid, beg, end;
	int baseQ = 0, mapQ = 0, min_len = 0, min_support = 1, min_support_strand = 0, min_supp_len = 0;
	int is_vcf = 0, var_only = 0, show_2strand = 0, trim_len = 0, del_as_allele = 0, proper_only = 0, dedup_overlap = 0, nthreads = 1;
	int nsites = 0;
	bed_site_t *sites = NULL;
	double min_af = 0.0;
	char *reg = 0, *chr_end; // specified region
	char *fname = 0; // reference fasta
	faidx_t *fai = 0;
	sam_hdr_t *h = 0; // header of the 1st input
	void *bed = 0;
	ketopt_t o = KETOPT_INIT;

	// parse the command line
	while ((n = ketopt(&o, argc, argv, 1, "r:q:Q:l:f:p:vcCS:s:b:x:T:ea:yVPt:D", 0)) >= 0) {
		if (n == 'f') { fname = o.arg; fai = fai_load(fname); }
		else if (n == 'b') { if (bed) bed_destroy(bed); bed = bed_read(o.arg); }
		else if (n == 'x') { if (bed) bed_destroy(bed); bed = bed_read_vcf(o.arg); }
		else if (n == 'l') min_len = atoi(o.arg); // minimum query length
		else if (n == 'r') reg = strdup(o.arg);   // parsing a region requires a BAM header
		else if (n == 'Q') baseQ = atoi(o.arg);   // base quality threshold
		else if (n == 'q') mapQ = atoi(o.arg);    // mapping quality threshold
		else if (n == 's') min_support = atoi(o.arg);
		else if (n == 'a') min_support_strand = atoi(o.arg);
		else if (n == 'S') min_supp_len = atoi(o.arg);
		else if (n == 'v') var_only = 1;
		else if (n == 'c') is_vcf = var_only = 1;
		else if (n == 'C') show_2strand = 1;
		else if (n == 'T') trim_len = atoi(o.arg);
		else if (n == 'e') del_as_allele = 1;
		else if (n == 'p') min_af = atof(o.arg);
		else if (n == 'P') proper_only = 1;
		else if (n == 'D') dedup_overlap = 1;
		else if (n == 't') nthreads = atoi(o.arg);
		else if (n == 'y') mapQ = 20, baseQ = 20, min_support = 5, min_support_strand = 2, is_vcf = var_only = show_2strand = 1;
		else if (n == 'V') {
			puts(VERSION);
			return 0;
		}
	}
	if (min_support < 1) min_support = 1;
	if (is_vcf && fai == 0) {
		fprintf(stderr, "[E::%s] with option -c, the reference genome must be provided.\n", __func__);
		return 1;
	}
	if (o.ind == argc) {
		fprintf(stderr, "Usage: minipileup2 [options] in1.bam/cram [in2.bam/cram [...]]\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  General:\n");
		fprintf(stderr, "    -f FILE      reference genome FASTA (required for CRAM and VCF output) [null]\n");
		fprintf(stderr, "    -v           show variants only\n");
		fprintf(stderr, "    -c           output in the VCF format (force -v)\n");
		fprintf(stderr, "    -C           show count of each allele on both strands\n");
		fprintf(stderr, "    -e           use '*' to mark deleted bases\n");
		fprintf(stderr, "    -y           variant calling mode (-vcC -a2 -s5 -q30 -Q20)\n");
		fprintf(stderr, "    -V           print version number\n");
		fprintf(stderr, "  Alignment filter:\n");
		fprintf(stderr, "    -r STR       region in format of 'ctg:start-end' [null]\n");
		fprintf(stderr, "    -b FILE      BED or position list file to include [null]\n");
		fprintf(stderr, "    -x FILE      VCF/BCF of target sites; restrict pileup to these positions [null]\n");
		fprintf(stderr, "    -t INT       number of threads [1]\n");
		fprintf(stderr, "    -P           only consider properly paired reads for paired-end reads\n");
		fprintf(stderr, "    -D           deduplicate overlapping read pairs (count each fragment once)\n");
		fprintf(stderr, "    -q INT       minimum mapping quality [%d]\n", mapQ);
		fprintf(stderr, "    -l INT       minimum alignment length [%d]\n", min_len);
		fprintf(stderr, "    -S INT       minimum supplementary alignment length [0]\n");
		fprintf(stderr, "  Site filter:\n");
		fprintf(stderr, "    -Q INT       minimum base quality [%d]\n", baseQ);
		fprintf(stderr, "    -T INT       skip bases within INT-bp from either end of a read [0]\n");
		fprintf(stderr, "    -s INT       drop alleles with depth<INT [%d]\n", min_support);
		fprintf(stderr, "    -a INT       drop alleles with depth<INT on either strand [%d]\n", min_support_strand);
		fprintf(stderr, "    -p FLOAT     drop an allele if the allele fraction is below FLOAT [%g]\n", min_af);
		return 1;
	}

	// open the first file to read the shared header; threads will reopen all files
	n = argc - o.ind;
	beg = 0; end = 1<<30; tid = -1;
	{
		htsFile *fp0 = hts_open(argv[o.ind], "r");
		if (fname) hts_set_fai_filename(fp0, fname);
		h = sam_hdr_read(fp0);
		if (reg) {
			chr_end = (char*)hts_parse_reg(reg, &beg, &end);
			if (chr_end) {
				char c = *chr_end; *chr_end = 0;
				tid = bam_name2id(h, reg);
				*chr_end = c;
			}
		}
		hts_close(fp0);
	}

	if (bed) sites = bed_get_sorted_sites(bed, h, &nsites);

	// print VCF header (main thread only, before workers start)
	if (is_vcf) {
		puts("##fileformat=VCFv4.2");
		printf("##source=minipileup-%s\n", VERSION);
		if (fai) {
			printf("##reference=%s\n", fname);
			int ni, nn = faidx_nseq(fai);
			for (ni = 0; ni < nn; ni++) {
				const char *seq = faidx_iseq(fai, ni);
				int len = faidx_seq_len(fai, seq);
				printf("##contig=<ID=%s,length=%d>\n", seq, len);
			}
		}
		puts("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
		if (show_2strand) {
			puts("##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand\">");
			puts("##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand\">");
		} else puts("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">");
		fputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", stdout);
		for (i = 0; i < n; ++i) printf("\t%s", argv[o.ind+i]);
		putchar('\n');
	}

	// decide actual thread count: only parallelize when we have discrete sites to split
	int nt = (nthreads > 1 && nsites > 1)? nthreads : 1;
	if (nt > nsites && nsites > 0) nt = nsites;

	thread_arg_t *targs = (thread_arg_t*)calloc(nt, sizeof(thread_arg_t));
	for (i = 0; i < nt; ++i) {
		targs[i].n_files    = n;
		targs[i].filenames  = argv + o.ind;
		targs[i].fname_ref  = fname;
		targs[i].h          = h;
		targs[i].bed        = bed;
		targs[i].sites      = sites;
		targs[i].nsites     = nsites;
		targs[i].si_beg     = (nsites > 0)? i * nsites / nt       : 0;
		targs[i].si_end     = (nsites > 0)? (i+1) * nsites / nt   : 1;
		targs[i].tid        = tid;
		targs[i].beg        = beg;
		targs[i].end        = end;
		targs[i].baseQ      = baseQ;
		targs[i].mapQ       = mapQ;
		targs[i].min_len    = min_len;
		targs[i].min_supp_len = min_supp_len;
		targs[i].proper_only  = proper_only;
		targs[i].is_vcf       = is_vcf;
		targs[i].var_only     = var_only;
		targs[i].show_2strand = show_2strand;
		targs[i].trim_len     = trim_len;
		targs[i].del_as_allele = del_as_allele;
		targs[i].min_support  = min_support;
		targs[i].min_support_strand = min_support_strand;
		targs[i].dedup_overlap = dedup_overlap;
		targs[i].min_af       = min_af;
		if (nt > 1)
			targs[i].out = open_memstream(&targs[i].out_buf, &targs[i].out_len);
		else
			targs[i].out = stdout;
	}

	if (nt == 1) {
		process_sites(&targs[0]);
	} else {
		pthread_t *threads = (pthread_t*)calloc(nt, sizeof(pthread_t));
		for (i = 0; i < nt; ++i)
			pthread_create(&threads[i], NULL, worker, &targs[i]);
		for (i = 0; i < nt; ++i) {
			pthread_join(threads[i], NULL);
			fclose(targs[i].out);
			if (targs[i].out_len > 0)
				fwrite(targs[i].out_buf, 1, targs[i].out_len, stdout);
			free(targs[i].out_buf);
		}
		free(threads);
	}
	free(targs);

	sam_hdr_destroy(h);
	if (fai) fai_destroy(fai);
	free(sites);
	free(reg);
	if (bed) bed_destroy(bed);

	fprintf(stderr, "[M::%s] Version: %s\n", __func__, VERSION);
	fprintf(stderr, "[M::%s] CMD:", __func__);
	for (i = 0; i < argc; ++i)
		fprintf(stderr, " %s", argv[i]);
	fprintf(stderr, "\n");
	return 0;
}

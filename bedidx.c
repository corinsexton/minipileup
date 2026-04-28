#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include "htslib/sam.h"

#include "ksort.h"
KSORT_INIT_GENERIC(uint64_t)

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 8192)

typedef struct {
	int n, m;
	uint64_t *a;
	int *idx;
} bed_reglist_t;

#include "khash.h"
KHASH_MAP_INIT_STR(reg, bed_reglist_t)

#define LIDX_SHIFT 13

typedef kh_reg_t reghash_t;

int *bed_index_core(int n, uint64_t *a, int *n_idx)
{
	int i, j, m, *idx;
	m = *n_idx = 0; idx = 0;
	for (i = 0; i < n; ++i) {
		int beg, end;
		beg = a[i]>>32 >> LIDX_SHIFT; end = ((uint32_t)a[i]) >> LIDX_SHIFT;
		if (m < end + 1) {
			int oldm = m;
			m = end + 1;
			kroundup32(m);
			idx = (int*)realloc(idx, m * sizeof(int));
			for (j = oldm; j < m; ++j) idx[j] = -1;
		}
		if (beg == end) {
			if (idx[beg] < 0) idx[beg] = i;
		} else {
			for (j = beg; j <= end; ++j)
				if (idx[j] < 0) idx[j] = i;
		}
		*n_idx = end + 1;
	}
	return idx;
}

void bed_index(void *_h)
{
	reghash_t *h = (reghash_t*)_h;
	khint_t k;
	for (k = 0; k < kh_end(h); ++k) {
		if (kh_exist(h, k)) {
			bed_reglist_t *p = &kh_val(h, k);
			if (p->idx) free(p->idx);
			ks_introsort(uint64_t, p->n, p->a);
			p->idx = bed_index_core(p->n, p->a, &p->m);
		}
	}
}

int bed_overlap_core(const bed_reglist_t *p, int beg, int end)
{
	int i, min_off;
	if (p->n == 0) return 0;
	min_off = (beg>>LIDX_SHIFT >= p->n)? p->idx[p->n-1] : p->idx[beg>>LIDX_SHIFT];
	if (min_off < 0) { // TODO: this block can be improved, but speed should not matter too much here
		int n = beg>>LIDX_SHIFT;
		if (n > p->n) n = p->n;
		for (i = n - 1; i >= 0; --i)
			if (p->idx[i] >= 0) break;
		min_off = i >= 0? p->idx[i] : 0;
	}
	for (i = min_off; i < p->n; ++i) {
		if ((int)(p->a[i]>>32) >= end) break; // out of range; no need to proceed
		if ((int32_t)p->a[i] > beg && (int32_t)(p->a[i]>>32) < end)
			return 1; // find the overlap; return
	}
	return 0;
}

int bed_overlap(const void *_h, const char *chr, int beg, int end)
{
	const reghash_t *h = (const reghash_t*)_h;
	khint_t k;
	if (!h) return 0;
	k = kh_get(reg, h, chr);
	if (k == kh_end(h)) return 0;
	return bed_overlap_core(&kh_val(h, k), beg, end);
}

void *bed_read(const char *fn)
{
	reghash_t *h = kh_init(reg);
	gzFile fp;
	kstream_t *ks;
	int dret;
	kstring_t *str;
	// read the list
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, str, &dret) >= 0) { // read the chr name
		int beg = -1, end = -1;
		bed_reglist_t *p;
		khint_t k = kh_get(reg, h, str->s);
		if (k == kh_end(h)) { // absent from the hash table
			int ret;
			char *s = strdup(str->s);
			k = kh_put(reg, h, s, &ret);
			memset(&kh_val(h, k), 0, sizeof(bed_reglist_t));
		}
		p = &kh_val(h, k);
		if (dret != '\n') { // if the lines has other characters
			if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
				beg = atoi(str->s); // begin
				if (dret != '\n') {
					if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
						end = atoi(str->s); // end
						if (end < beg) end = -1;
					}
				}
			}
		}
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n'); // skip the rest of the line
		if (end < 0 && beg > 0) end = beg, beg = beg - 1; // if there is only one column
		if (beg >= 0 && end > beg) {
			if (p->n == p->m) {
				p->m = p->m? p->m<<1 : 4;
				p->a = (uint64_t*)realloc(p->a, p->m * 8);
			}
			p->a[p->n++] = (uint64_t)beg<<32 | end;
		}
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	bed_index(h);
	return h;
}

void *bed_read_vcf(const char *fn)
{
	reghash_t *h = kh_init(reg);
	gzFile fp;
	kstream_t *ks;
	int dret;
	kstring_t *str;
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) return 0;
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, str, &dret) >= 0) {
		char *line = str->s;
		if (str->l == 0 || line[0] == '#') continue; // skip VCF header lines
		// parse CHROM (field 0) and POS (field 1, 1-based)
		char *tab1 = strchr(line, '\t');
		if (!tab1) continue;
		*tab1 = '\0';
		char *pos_field = tab1 + 1;
		char *tab2 = strchr(pos_field, '\t');
		if (tab2) *tab2 = '\0';
		if (!isdigit((unsigned char)pos_field[0])) continue;
		int pos = atoi(pos_field);
		if (pos <= 0) continue;
		int beg = pos - 1, end = pos; // convert 1-based VCF POS to 0-based half-open
		bed_reglist_t *p;
		khint_t k = kh_get(reg, h, line);
		if (k == kh_end(h)) {
			int ret;
			char *s = strdup(line);
			k = kh_put(reg, h, s, &ret);
			memset(&kh_val(h, k), 0, sizeof(bed_reglist_t));
		}
		p = &kh_val(h, k);
		if (p->n == p->m) {
			p->m = p->m? p->m<<1 : 4;
			p->a = (uint64_t*)realloc(p->a, p->m * 8);
		}
		p->a[p->n++] = (uint64_t)beg<<32 | end;
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	bed_index(h);
	return h;
}

typedef struct { int tid; int beg, end; } bed_site_t;

static int bed_site_cmp(const void *a, const void *b)
{
	const bed_site_t *x = (const bed_site_t*)a;
	const bed_site_t *y = (const bed_site_t*)b;
	if (x->tid != y->tid) return x->tid - y->tid;
	return x->beg - y->beg;
}

// Returns a sorted array of all (tid, beg, end) intervals from the bed hash.
// Chromosomes absent from hdr are silently skipped.  Caller frees the result.
bed_site_t *bed_get_sorted_sites(const void *_h, sam_hdr_t *hdr, int *n_out)
{
	const reghash_t *h = (const reghash_t*)_h;
	khint_t k;
	int n = 0, m = 0, i;
	bed_site_t *sites = NULL;
	for (k = 0; k < kh_end(h); ++k) {
		if (!kh_exist(h, k)) continue;
		const char *chr = kh_key(h, k);
		int tid = bam_name2id(hdr, chr);
		if (tid < 0) continue;
		const bed_reglist_t *p = &kh_val(h, k);
		for (i = 0; i < p->n; ++i) {
			if (n == m) {
				m = m ? m << 1 : 16;
				sites = (bed_site_t*)realloc(sites, m * sizeof(bed_site_t));
			}
			sites[n].tid = tid;
			sites[n].beg = (int)(p->a[i] >> 32);
			sites[n].end = (int)(uint32_t)p->a[i];
			++n;
		}
	}
	if (n > 0) qsort(sites, n, sizeof(bed_site_t), bed_site_cmp);
	*n_out = n;
	return sites;
}

void bed_destroy(void *_h)
{
	reghash_t *h = (reghash_t*)_h;
	khint_t k;
	for (k = 0; k < kh_end(h); ++k) {
		if (kh_exist(h, k)) {
			free(kh_val(h, k).a);
			free(kh_val(h, k).idx);
			free((char*)kh_key(h, k));
		}
	}
	kh_destroy(reg, h);
}

// Build an hts_reglist_t array from the bed hash for use with sam_itr_regions().
// The returned array's reg fields borrow pointers from the bed hash (not copied),
// so the bed hash must outlive the iterator. intervals[] arrays are heap-allocated
// and freed by hts_reglist_free() when the iterator is destroyed.
hts_reglist_t *bed_to_reglist(const void *_h, sam_hdr_t *hdr, int *nregs)
{
	const reghash_t *h = (const reghash_t*)_h;
	khint_t k;
	int n = 0, i;
	hts_reglist_t *reglist;

	for (k = 0; k < kh_end(h); ++k)
		if (kh_exist(h, k) && kh_val(h, k).n > 0) ++n;

	*nregs = 0;
	if (n == 0) return 0;

	reglist = (hts_reglist_t*)calloc(n, sizeof(hts_reglist_t));
	if (!reglist) return 0;
	n = 0;
	for (k = 0; k < kh_end(h); ++k) {
		if (!kh_exist(h, k)) continue;
		const bed_reglist_t *p = &kh_val(h, k);
		if (p->n == 0) continue;
		reglist[n].reg = kh_key(h, k);
		reglist[n].count = (uint32_t)p->n;
		reglist[n].intervals = (hts_pair_pos_t*)malloc(p->n * sizeof(hts_pair_pos_t));
		if (!reglist[n].intervals) {
			for (i = 0; i < n; ++i) free(reglist[i].intervals);
			free(reglist);
			return 0;
		}
		reglist[n].min_beg = HTS_POS_MAX;
		reglist[n].max_end = 0;
		for (i = 0; i < p->n; ++i) {
			hts_pos_t beg = (hts_pos_t)(p->a[i] >> 32);
			hts_pos_t end = (hts_pos_t)(uint32_t)p->a[i];
			reglist[n].intervals[i].beg = beg;
			reglist[n].intervals[i].end = end;
			if (beg < reglist[n].min_beg) reglist[n].min_beg = beg;
			if (end > reglist[n].max_end) reglist[n].max_end = end;
		}
		++n;
	}
	*nregs = n;
	return reglist;
}

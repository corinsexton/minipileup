#!/bin/bash
# Edge-case tests: minipileup2 single-thread vs multi-thread must be identical.
#
# Cases covered:
#   1. Thread count sweep: t=2,4,8,16,50,300 (300 >> 202 sites = idle threads)
#   2. Flag variations: -e, -P, -T5, -s2, -a1, -p0.05 each with t=1 vs t=4
#   3. Single-sample input (n=1 CRAM)
#   4. Rare chromosomes: chrX and chrY only
#   5. Single chromosome: chr1 only (dense, 20 sites)

set -euo pipefail

CRAM1="${1:-/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/SMHT007/SMHT007-3AK/illuminaNovaseq_bulkWgs/seq_data/SMAFI9SM15L9.cram}"
CRAM2="${2:-/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/SMHT007/SMHT007-3AK/illuminaNovaseq_bulkWgs/seq_data/SMAFIIRDMCE5.cram}"
VCF="/n/data1/hms/dbmi/park-smaht_dac/DATA/Globus_DataShare/SNV_Indel/DAC/somatic_snvIndel/p25/filtered_calls/SMHT007/SMAFIF98HJA2.vcf.gz"
REF=~/hg38_no_alt.fa
MP2=/n/data1/hms/dbmi/park/corinne/smaht/minipileup/minipileup2

TESTDIR=$(mktemp -d)
echo "Working directory: $TESTDIR"
echo ""

PASS=0; FAIL=0; FAIL_CASES=""

# ---------------------------------------------------------------------------
# Helper: run one comparison and report PASS/FAIL.
# Usage: run_case LABEL FLAGS_T1 FLAGS_TN NTHREADS CRAM_ARGS VCF_PATH
# ---------------------------------------------------------------------------
run_case() {
    local label="$1"
    local flags="$2"         # shared flags (applied to both t=1 and t=N)
    local nt="$3"            # thread count for multi-thread run
    local cram_args="$4"     # space-separated CRAM paths
    local vcf_path="$5"

    local outdir="$TESTDIR/${label//[^a-zA-Z0-9_]/_}"
    mkdir -p "$outdir"

    # t=1
    # shellcheck disable=SC2086
    "$MP2" -f "$REF" $flags -x "$vcf_path" $cram_args \
        2>"$outdir/log_t1.txt" | grep -v '^#' > "$outdir/out_t1.txt" || true

    # t=N
    # shellcheck disable=SC2086
    "$MP2" -t "$nt" -f "$REF" $flags -x "$vcf_path" $cram_args \
        2>"$outdir/log_tN.txt" | grep -v '^#' > "$outdir/out_tN.txt" || true

    # Normalize both outputs
    normalize() {
        local infile=$1 outfile=$2
        awk '
        {
            if ($9 !~ /ADF/ && $9 !~ /AD/) next
            split($9, fmt, ":")
            adf_idx = adr_idx = ad_idx = 0
            for (i=1; i<=length(fmt); i++) {
                if (fmt[i]=="ADF") adf_idx=i
                if (fmt[i]=="ADR") adr_idx=i
                if (fmt[i]=="AD")  ad_idx=i
            }
            n_samples = NF - 9
            for (s=1; s<=n_samples; s++) {
                split($(9+s), val, ":")
                if (adf_idx && adr_idx) {
                    split(val[adf_idx], adf, ",")
                    split(val[adr_idx], adr, ",")
                    print $1, $2, $4, $5, s, adf[1], adf[2], adr[1], adr[2]
                } else if (ad_idx) {
                    split(val[ad_idx], ad, ",")
                    print $1, $2, $4, $5, s, ad[1], ad[2], 0, 0
                }
            }
        }' "$infile" | sort > "$outfile"
    }

    normalize "$outdir/out_t1.txt" "$outdir/norm_t1.txt"
    normalize "$outdir/out_tN.txt" "$outdir/norm_tN.txt"

    local n1 nN
    n1=$(wc -l < "$outdir/norm_t1.txt")
    nN=$(wc -l < "$outdir/norm_tN.txt")

    # Python diff
    result=$(python3 - "$outdir" "$nt" <<'PYEOF'
import sys, os
def load(p):
    d = {}
    with open(p) as f:
        for line in f:
            p2 = line.strip().split()
            if len(p2) >= 9:
                d[(p2[0],p2[1],p2[2],p2[3],p2[4])] = tuple(p2[5:9])
    return d
testdir, nt = sys.argv[1], sys.argv[2]
t1 = load(os.path.join(testdir, 'norm_t1.txt'))
tn = load(os.path.join(testdir, 'norm_tN.txt'))
only_t1 = set(t1) - set(tn)
only_tn = set(tn) - set(t1)
shared   = set(t1) & set(tn)
match    = sum(1 for k in shared if t1[k] == tn[k])
if len(only_t1)==0 and len(only_tn)==0 and match==len(shared):
    print(f"PASS entries={len(shared)}")
else:
    msgs = []
    if only_t1:
        ex = sorted(only_t1)[:3]
        msgs.append(f"only_t1={len(only_t1)} e.g. {ex[0]}")
    if only_tn:
        ex = sorted(only_tn)[:3]
        msgs.append(f"only_tN={len(only_tn)} e.g. {ex[0]}")
    if match < len(shared):
        mismatches = [(k, t1[k], tn[k]) for k in sorted(shared) if t1[k]!=tn[k]][:3]
        msgs.append(f"mismatch={len(shared)-match} e.g. {mismatches[0]}")
    print("FAIL " + "; ".join(msgs))
PYEOF
)

    local status="${result%% *}"
    local detail="${result#* }"
    printf "  %-55s t=1(%4d lines) t=%s(%4d lines)  %s\n" \
        "$label" "$n1" "$nt" "$nN" "$result"

    if [[ "$status" == "PASS" ]]; then
        ((PASS++)) || true
    else
        ((FAIL++)) || true
        FAIL_CASES="$FAIL_CASES\n  $label: $detail"
    fi
}

# ---------------------------------------------------------------------------
# Prepare VCFs for each chromosome filter
# ---------------------------------------------------------------------------
echo "Preparing per-chromosome VCFs..."
ALL_VCF="$TESTDIR/all.vcf.gz"
zcat "$VCF" | bgzip > "$ALL_VCF"
tabix -p vcf "$ALL_VCF"
NSITES=$(bcftools view -H "$ALL_VCF" | wc -l)
echo "  Total sites: $NSITES"

# VCF with only chrX+chrY
CHRXY_VCF="$TESTDIR/chrXY.vcf.gz"
bcftools view -r chrX,chrY "$ALL_VCF" | bgzip > "$CHRXY_VCF"
tabix -p vcf "$CHRXY_VCF"
N_XY=$(bcftools view -H "$CHRXY_VCF" | wc -l)

# VCF with only chr1
CHR1_VCF="$TESTDIR/chr1.vcf.gz"
bcftools view -r chr1 "$ALL_VCF" | bgzip > "$CHR1_VCF"
tabix -p vcf "$CHR1_VCF"
N_CHR1=$(bcftools view -H "$CHR1_VCF" | wc -l)
echo "  chrX+chrY sites: $N_XY  |  chr1 sites: $N_CHR1"
echo ""

BASE_FLAGS="-c -C -Q 20 -q 30 -s 0"
TWO_CRAMS="$CRAM1 $CRAM2"

# ---------------------------------------------------------------------------
echo "=== Case 1: Thread count sweep (t=2,4,8,16,50,300) ==="
echo "  (300 threads >> $NSITES sites exercises idle-thread path)"
for nt in 2 4 8 16 50 300; do
    run_case "t=1_vs_t=${nt}_all_sites_base_flags" "$BASE_FLAGS" "$nt" "$TWO_CRAMS" "$ALL_VCF"
done
echo ""

# ---------------------------------------------------------------------------
echo "=== Case 2: Flag variations (t=1 vs t=4) ==="

# -e: treat deletions as '*' allele
run_case "flags_del_as_allele(-e)" "-c -C -Q 20 -q 30 -s 0 -e" 4 "$TWO_CRAMS" "$ALL_VCF"

# -P: proper pairs only
run_case "flags_proper_pairs(-P)" "-c -C -Q 20 -q 30 -s 0 -P" 4 "$TWO_CRAMS" "$ALL_VCF"

# -T 5: trim 5 bp from read ends
run_case "flags_trim_len(-T5)" "-c -C -Q 20 -q 30 -s 0 -T 5" 4 "$TWO_CRAMS" "$ALL_VCF"

# -s 2: require >=2 reads supporting an allele
run_case "flags_min_support(-s2)" "-c -C -Q 20 -q 30 -s 2" 4 "$TWO_CRAMS" "$ALL_VCF"

# -a 1: require >=1 read on each strand
run_case "flags_min_strand(-a1)" "-c -C -Q 20 -q 30 -s 0 -a 1" 4 "$TWO_CRAMS" "$ALL_VCF"

# -p 0.05: drop alleles with AF < 5%
run_case "flags_min_af(-p0.05)" "-c -C -Q 20 -q 30 -s 0 -p 0.05" 4 "$TWO_CRAMS" "$ALL_VCF"

# -v only (no -c): variant-only text output (no VCF format, no ref needed same)
run_case "flags_var_only(-v_no_VCF_fmt)" "-v -Q 20 -q 30 -s 0" 4 "$TWO_CRAMS" "$ALL_VCF"

# -y: preset (-vcC -a2 -s5 -q30 -Q20)
run_case "flags_preset(-y)" "-y" 4 "$TWO_CRAMS" "$ALL_VCF"

# Looser baseQ/mapQ
run_case "flags_loose_quality(-Q10_-q10)" "-c -C -Q 10 -q 10 -s 0" 4 "$TWO_CRAMS" "$ALL_VCF"

# Combined strict filters
run_case "flags_strict_combined(-s2_-a1_-p0.05)" "-c -C -Q 20 -q 30 -s 2 -a 1 -p 0.05" 4 "$TWO_CRAMS" "$ALL_VCF"
echo ""

# ---------------------------------------------------------------------------
echo "=== Case 3: Single-sample input (one CRAM) ==="
for nt in 4 16 300; do
    run_case "single_sample_t=${nt}" "$BASE_FLAGS" "$nt" "$CRAM1" "$ALL_VCF"
done
echo ""

# ---------------------------------------------------------------------------
echo "=== Case 4: Chromosome-restricted VCFs ==="
echo "  chrX+chrY ($N_XY sites) — sex chromosomes, unusual coverage patterns"
for nt in 4 16; do
    run_case "chrXY_t=${nt}" "$BASE_FLAGS" "$nt" "$TWO_CRAMS" "$CHRXY_VCF"
done

echo "  chr1 only ($N_CHR1 sites)"
for nt in 4 16; do
    run_case "chr1_t=${nt}" "$BASE_FLAGS" "$nt" "$TWO_CRAMS" "$CHR1_VCF"
done
echo ""

# ---------------------------------------------------------------------------
echo "=== Case 5: Thread count > site count on small VCFs ==="
echo "  chrXY ($N_XY sites) with t=50 and t=300"
for nt in 50 300; do
    run_case "chrXY_more_threads_than_sites_t=${nt}" "$BASE_FLAGS" "$nt" "$TWO_CRAMS" "$CHRXY_VCF"
done
echo "  chr1 ($N_CHR1 sites) with t=50 and t=300"
for nt in 50 300; do
    run_case "chr1_more_threads_than_sites_t=${nt}" "$BASE_FLAGS" "$nt" "$TWO_CRAMS" "$CHR1_VCF"
done
echo ""

# ---------------------------------------------------------------------------
echo "=== Case 6: -D (overlap deduplication) thread correctness ==="
BASE_FLAGS_D="$BASE_FLAGS -D"

echo "  Thread sweep (t=2,4,8,16,50,300) with -D"
for nt in 2 4 8 16 50 300; do
    run_case "D_t=1_vs_t=${nt}_all_sites" "$BASE_FLAGS_D" "$nt" "$TWO_CRAMS" "$ALL_VCF"
done

echo "  -D combined with other filters (t=1 vs t=4)"
run_case "D_proper_pairs(-D_-P)"           "$BASE_FLAGS_D -P"        4 "$TWO_CRAMS" "$ALL_VCF"
run_case "D_trim_len(-D_-T5)"              "$BASE_FLAGS_D -T 5"      4 "$TWO_CRAMS" "$ALL_VCF"
run_case "D_min_support(-D_-s2)"           "-c -C -Q 20 -q 30 -s 2 -D" 4 "$TWO_CRAMS" "$ALL_VCF"
run_case "D_loose_quality(-D_-Q10_-q10)"   "-c -C -Q 10 -q 10 -s 0 -D" 4 "$TWO_CRAMS" "$ALL_VCF"

echo "  -D single-sample (t=4, t=16)"
for nt in 4 16; do
    run_case "D_single_sample_t=${nt}" "$BASE_FLAGS_D" "$nt" "$CRAM1" "$ALL_VCF"
done

echo "  -D threads > sites (chrXY t=50, t=300)"
for nt in 50 300; do
    run_case "D_chrXY_t=${nt}" "$BASE_FLAGS_D" "$nt" "$TWO_CRAMS" "$CHRXY_VCF"
done
echo ""

# ---------------------------------------------------------------------------
echo "=== SUMMARY ==="
echo "PASS: $PASS"
echo "FAIL: $FAIL"
if [[ $FAIL -gt 0 ]]; then
    echo -e "Failed cases:$FAIL_CASES"
    exit 1
else
    echo "All edge cases passed."
fi
echo ""
echo "Full output: $TESTDIR"

#!/bin/bash
# Test that minipileup2 single-thread and multi-thread produce identical output.
set -euo pipefail

CRAM1="${1:-/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/SMHT007/SMHT007-3AK/illuminaNovaseq_bulkWgs/seq_data/SMAFI9SM15L9.cram}"
CRAM2="${2:-/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/SMHT007/SMHT007-3AK/illuminaNovaseq_bulkWgs/seq_data/SMAFIIRDMCE5.cram}"
VCF="/n/data1/hms/dbmi/park-smaht_dac/DATA/Globus_DataShare/SNV_Indel/DAC/somatic_snvIndel/p25/filtered_calls/SMHT007/SMAFIF98HJA2.vcf.gz"
REF=~/hg38_no_alt.fa
MP2=/n/data1/hms/dbmi/park/corinne/smaht/minipileup/minipileup2
MINPILEUP_ARGS="-c -C -Q 20 -q 30 -s 0"
N_SITES="${3:-100}"
NTHREADS="${4:-4}"

TESTDIR=$(mktemp -d)
echo "Working directory: $TESTDIR"
echo "CRAM1: $(basename $CRAM1)"
echo "CRAM2: $(basename $CRAM2)"
echo "N_SITES: $N_SITES  NTHREADS: $NTHREADS"
echo ""

# Step 1: Subset VCF
echo "=== Step 1: Subset VCF to $N_SITES sites ==="
set +o pipefail
zcat "$VCF" | awk '/^#/{print; next} {print; n++; if(n>='$N_SITES') exit}' \
    | bgzip > "$TESTDIR/sites.vcf.gz"
set -o pipefail
tabix -p vcf "$TESTDIR/sites.vcf.gz"
echo "$(bcftools view -H "$TESTDIR/sites.vcf.gz" | wc -l) variant sites"
echo ""

# Step 2: minipileup2 single-thread
echo "=== Step 2: minipileup2 t=1 ==="
T0=$(date +%s%N)
"$MP2" -f "$REF" $MINPILEUP_ARGS -x "$TESTDIR/sites.vcf.gz" "$CRAM1" "$CRAM2" \
    2>"$TESTDIR/log_t1.txt" \
    | grep -v '^#' > "$TESTDIR/out_t1.txt"
T1=$(date +%s%N)
echo "Time: $(( (T1 - T0) / 1000000 )) ms"
echo "Output lines: $(wc -l < "$TESTDIR/out_t1.txt")"
echo ""

# Step 3: minipileup2 multi-thread
echo "=== Step 3: minipileup2 t=$NTHREADS ==="
T0=$(date +%s%N)
"$MP2" -t "$NTHREADS" -f "$REF" $MINPILEUP_ARGS -x "$TESTDIR/sites.vcf.gz" "$CRAM1" "$CRAM2" \
    2>"$TESTDIR/log_t${NTHREADS}.txt" \
    | grep -v '^#' > "$TESTDIR/out_t${NTHREADS}.txt"
T1=$(date +%s%N)
echo "Time: $(( (T1 - T0) / 1000000 )) ms"
echo "Output lines: $(wc -l < "$TESTDIR/out_t${NTHREADS}.txt")"
echo ""

# Step 4: Normalize and compare
echo "=== Step 4: Normalize + compare ==="

extract_counts() {
    local infile=$1 outfile=$2
    awk '
    {
        if ($9 !~ /ADF/) next
        split($9, fmt, ":")
        adf_idx = adr_idx = 0
        for (i=1; i<=length(fmt); i++) {
            if (fmt[i]=="ADF") adf_idx=i
            if (fmt[i]=="ADR") adr_idx=i
        }
        if (!adf_idx || !adr_idx) next
        n_samples = NF - 9
        for (s=1; s<=n_samples; s++) {
            split($(9+s), val, ":")
            split(val[adf_idx], adf, ",")
            split(val[adr_idx], adr, ",")
            print $1, $2, $4, $5, s, adf[1], adf[2], adr[1], adr[2]
        }
    }' "$infile" | sort > "$outfile"
}

extract_counts "$TESTDIR/out_t1.txt"        "$TESTDIR/norm_t1.txt"
extract_counts "$TESTDIR/out_t${NTHREADS}.txt" "$TESTDIR/norm_t${NTHREADS}.txt"

echo "  t=1 entries:           $(wc -l < "$TESTDIR/norm_t1.txt")"
echo "  t=$NTHREADS entries:           $(wc -l < "$TESTDIR/norm_t${NTHREADS}.txt")"
echo ""

python3 - "$TESTDIR" "$NTHREADS" <<'PYEOF'
import sys, os

def load_norm(path):
    d = {}
    with open(path) as f:
        for line in f:
            p = line.strip().split()
            if len(p) >= 9:
                key = (p[0], p[1], p[2], p[3], p[4])
                val = (int(p[5]), int(p[6]), int(p[7]), int(p[8]))
                d[key] = val
    return d

testdir, nt = sys.argv[1], sys.argv[2]
t1 = load_norm(os.path.join(testdir, 'norm_t1.txt'))
tn = load_norm(os.path.join(testdir, f'norm_t{nt}.txt'))

only_t1 = set(t1) - set(tn)
only_tn = set(tn) - set(t1)
shared   = set(t1) & set(tn)
match    = sum(1 for k in shared if t1[k] == tn[k])

print(f"  Only in t=1:    {len(only_t1)}")
print(f"  Only in t={nt}: {len(only_tn)}")
print(f"  Shared entries: {len(shared)}")
print(f"  Exact match:    {match} / {len(shared)}")

if only_t1:
    print(f"  Entries only in t=1 (first 5):")
    for k in sorted(only_t1)[:5]:
        print(f"    {k[0]} {k[1]} {k[2]} {k[3]} sample={k[4]} {t1[k]}")
if only_tn:
    print(f"  Entries only in t={nt} (first 5):")
    for k in sorted(only_tn)[:5]:
        print(f"    {k[0]} {k[1]} {k[2]} {k[3]} sample={k[4]} {tn[k]}")
if match < len(shared):
    print(f"  Count mismatches (first 5):")
    n = 0
    for k in sorted(shared):
        if t1[k] != tn[k]:
            print(f"    {k[0]} {k[1]} {k[2]} {k[3]} sample={k[4]}  t1:{t1[k]}  t{nt}:{tn[k]}")
            n += 1
            if n >= 5: break

print()
if len(only_t1) == 0 and len(only_tn) == 0 and match == len(shared):
    print(f"PASS: t=1 and t={nt} produce identical output ({match} entries match).")
else:
    print(f"FAIL: t=1 and t={nt} differ (only_t1={len(only_t1)}, only_t{nt}={len(only_tn)}, match={match}/{len(shared)}).")

print(f"\nFull output: {testdir}")
PYEOF

#!/bin/bash
# Compare minipileup (original), samtools mpileup, and minipileup2 on N_SITES.
# All three tools are run with two CRAM/BAM inputs simultaneously.
set -euo pipefail

CRAM1="${1:-/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/SMHT007/SMHT007-3AK/illuminaNovaseq_bulkWgs/seq_data/SMAFI9SM15L9.cram}"
CRAM2="${2:-/n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/SMHT007/SMHT007-3AK/illuminaNovaseq_bulkWgs/seq_data/SMAFIIRDMCE5.cram}"
VCF="/n/data1/hms/dbmi/park-smaht_dac/DATA/Globus_DataShare/SNV_Indel/DAC/somatic_snvIndel/p25/filtered_calls/SMHT007/SMAFIF98HJA2.vcf.gz"
REF=~/hg38_no_alt.fa
MP1=~/software/bin/minipileup
MP2=/n/data1/hms/dbmi/park/corinne/smaht/minipileup/minipileup2
# Shared pileup flags: VCF output, per-strand counts, baseQ>=20, mapQ>=30, no min-allele-depth filter
MINPILEUP_ARGS="-c -C -Q 20 -q 30 -s 0"
N_SITES=10000

TESTDIR=$(mktemp -d)
echo "Working directory: $TESTDIR"
echo "CRAM1: $(basename $CRAM1)"
echo "CRAM2: $(basename $CRAM2)"
echo ""

# --------------------------------------------------------------------------
# Step 1: Subset VCF to N_SITES
# --------------------------------------------------------------------------
echo "=== Step 1: Subset VCF to $N_SITES sites ==="
set +o pipefail
zcat "$VCF" | awk '/^#/{print; next} {print; n++; if(n>='$N_SITES') exit}' \
    | bgzip > "$TESTDIR/sites.vcf.gz"
set -o pipefail
tabix -p vcf "$TESTDIR/sites.vcf.gz"

bcftools query -f '%CHROM:%POS-%END\n' "$TESTDIR/sites.vcf.gz" > "$TESTDIR/regions.txt"
# chr, 1-based POS, REF, ALT (first ALT only) for mpileup parser
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$TESTDIR/sites.vcf.gz" \
    | awk '{split($4,a,","); print $1,$2,$3,a[1]}' > "$TESTDIR/sites_lookup.txt"

echo "$(wc -l < "$TESTDIR/regions.txt") regions written"
echo ""

# --------------------------------------------------------------------------
# Step 2: Original minipileup (per-region BAMs from both CRAMs, -r flag)
# --------------------------------------------------------------------------
echo "=== Step 2: minipileup original (two per-region BAMs) ==="
> "$TESTDIR/out_mp1.txt"
> "$TESTDIR/log_mp1.txt"
T0=$(date +%s%N)
while IFS= read -r region; do
    safe="${region//:/_}"; safe="${safe//-/__}"
    bam1="$TESTDIR/s1_${safe}.bam"
    bam2="$TESTDIR/s2_${safe}.bam"
    samtools view --write-index --reference "$REF" -b -o "$bam1" "$CRAM1" "$region" \
        2>>"$TESTDIR/log_mp1.txt"
    samtools view --write-index --reference "$REF" -b -o "$bam2" "$CRAM2" "$region" \
        2>>"$TESTDIR/log_mp1.txt"
    set +e
    "$MP1" -f "$REF" $MINPILEUP_ARGS -r "$region" "$bam1" "$bam2" \
        2>>"$TESTDIR/log_mp1.txt" \
        | grep -v '^#' >> "$TESTDIR/out_mp1.txt"
    set -e
    rm -f "$bam1" "${bam1}.csi" "$bam2" "${bam2}.csi"
done < "$TESTDIR/regions.txt"
T1=$(date +%s%N)
echo "Time: $(( (T1 - T0) / 1000000 )) ms"
echo "Output lines: $(wc -l < "$TESTDIR/out_mp1.txt")"
echo ""

# --------------------------------------------------------------------------
# Step 3: samtools mpileup (per-region BAMs from both CRAMs)
# --------------------------------------------------------------------------
echo "=== Step 3: samtools mpileup (two per-region BAMs) ==="
> "$TESTDIR/mpileup_raw.txt"
> "$TESTDIR/log_smp.txt"
T0=$(date +%s%N)
while IFS= read -r region; do
    safe="${region//:/_}"; safe="${safe//-/__}"
    bam1="$TESTDIR/s1_${safe}.bam"
    bam2="$TESTDIR/s2_${safe}.bam"
    samtools view --write-index --reference "$REF" -b -o "$bam1" "$CRAM1" "$region" \
        2>>"$TESTDIR/log_smp.txt"
    samtools view --write-index --reference "$REF" -b -o "$bam2" "$CRAM2" "$region" \
        2>>"$TESTDIR/log_smp.txt"
    samtools mpileup -f "$REF" -q 30 -Q 20 -B -x -d 0 -r "$region" "$bam1" "$bam2" \
        2>>"$TESTDIR/log_smp.txt" >> "$TESTDIR/mpileup_raw.txt" || true
    rm -f "$bam1" "${bam1}.csi" "$bam2" "${bam2}.csi"
done < "$TESTDIR/regions.txt"
T1=$(date +%s%N)
echo "Time: $(( (T1 - T0) / 1000000 )) ms"
echo "Raw mpileup lines: $(wc -l < "$TESTDIR/mpileup_raw.txt")"

# Parse samtools mpileup output into VCF-like multi-sample format
python3 - "$TESTDIR" <<'PYEOF'
import sys, re, os

testdir = sys.argv[1]

lookup = {}
with open(os.path.join(testdir, 'sites_lookup.txt')) as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 4:
            continue
        chrom, pos, ref, alt = parts[0], int(parts[1]), parts[2], parts[3]
        lookup[(chrom, pos)] = (ref, alt)

def clean_reads(reads_str):
    s = re.sub(r'\^.', '', reads_str)
    s = s.replace('$', '')
    while True:
        m = re.search(r'[+](\d+)', s)
        if not m:
            break
        n = int(m.group(1))
        s = s[:m.start()] + s[m.start() + len(m.group(0)) + n:]
    while True:
        m = re.search(r'[-](\d+)', s)
        if not m:
            break
        n = int(m.group(1))
        s = s[:m.start()] + s[m.start() + len(m.group(0)) + n:]
    return s

out_lines = []
with open(os.path.join(testdir, 'mpileup_raw.txt')) as fin:
    for line in fin:
        parts = line.rstrip('\n').split('\t')
        # mpileup with N samples: CHROM POS REF (DEPTH READS QUALS)*N
        if len(parts) < 6:
            continue
        chrom, pos_str, ref_base = parts[0], parts[1], parts[2]
        pos = int(pos_str)
        if (chrom, pos) not in lookup:
            continue
        vcf_ref, vcf_alt = lookup[(chrom, pos)]
        if len(vcf_ref) != 1 or len(vcf_alt) != 1:
            continue
        n_samples = (len(parts) - 3) // 3
        sample_fields = []
        any_alt = False
        for s in range(n_samples):
            reads_str = parts[4 + s * 3]
            cleaned = clean_reads(reads_str)
            ref_fwd = cleaned.count('.')
            ref_rev = cleaned.count(',')
            alt_fwd = cleaned.count(vcf_alt.upper())
            alt_rev = cleaned.count(vcf_alt.lower())
            if alt_fwd + alt_rev > 0:
                any_alt = True
            sample_fields.append(f'0/1:{ref_fwd},{alt_fwd}:{ref_rev},{alt_rev}')
        if not any_alt:
            continue
        row = f'{chrom}\t{pos}\t.\t{vcf_ref}\t{vcf_alt}\t.\t.\t.\tGT:ADF:ADR'
        for sf in sample_fields:
            row += f'\t{sf}'
        out_lines.append(row + '\n')

with open(os.path.join(testdir, 'out_smp.txt'), 'w') as fout:
    fout.writelines(out_lines)

print(f'  Parsed {len(out_lines)} variant sites from samtools mpileup')
PYEOF

echo ""

# --------------------------------------------------------------------------
# Step 4: minipileup2 single-thread (full CRAMs + subset VCF, indexed access)
# --------------------------------------------------------------------------
echo "=== Step 4a: minipileup2 single-thread (two CRAMs + VCF indexed) ==="
T0=$(date +%s%N)
"$MP2" -f "$REF" $MINPILEUP_ARGS -x "$TESTDIR/sites.vcf.gz" "$CRAM1" "$CRAM2" \
    2>"$TESTDIR/log_mp2_t1.txt" \
    | grep -v '^#' > "$TESTDIR/out_mp2_t1.txt"
T1=$(date +%s%N)
echo "Time: $(( (T1 - T0) / 1000000 )) ms"
echo "Output lines: $(wc -l < "$TESTDIR/out_mp2_t1.txt")"
echo ""

# --------------------------------------------------------------------------
# Step 4b: minipileup2 multi-thread (same input, -t 4)
# --------------------------------------------------------------------------
echo "=== Step 4b: minipileup2 multi-thread -t 4 (two CRAMs + VCF indexed) ==="
T0=$(date +%s%N)
"$MP2" -t 4 -f "$REF" $MINPILEUP_ARGS -x "$TESTDIR/sites.vcf.gz" "$CRAM1" "$CRAM2" \
    2>"$TESTDIR/log_mp2_t4.txt" \
    | grep -v '^#' > "$TESTDIR/out_mp2_t4.txt"
T1=$(date +%s%N)
echo "Time: $(( (T1 - T0) / 1000000 )) ms"
echo "Output lines: $(wc -l < "$TESTDIR/out_mp2_t4.txt")"
echo ""

# Use single-thread output as the canonical mp2 result for cross-tool comparison
cp "$TESTDIR/out_mp2_t1.txt" "$TESTDIR/out_mp2.txt"

# --------------------------------------------------------------------------
# Step 5: Normalise each tool's output to: chr pos ref alt sample rf af rr ar
# --------------------------------------------------------------------------
echo "=== Step 5: Normalise outputs ==="

# Extract all sample columns from VCF lines with GT:ADF:ADR FORMAT
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
            # chr pos ref alt sample rf af rr ar
            print $1, $2, $4, $5, s, adf[1], adf[2], adr[1], adr[2]
        }
    }' "$infile" | sort > "$outfile"
}

extract_counts "$TESTDIR/out_mp1.txt"      "$TESTDIR/norm_mp1.txt"
extract_counts "$TESTDIR/out_smp.txt"      "$TESTDIR/norm_smp.txt"
extract_counts "$TESTDIR/out_mp2_t1.txt"   "$TESTDIR/norm_mp2_t1.txt"
extract_counts "$TESTDIR/out_mp2_t4.txt"   "$TESTDIR/norm_mp2_t4.txt"
cp "$TESTDIR/norm_mp2_t1.txt" "$TESTDIR/norm_mp2.txt"

echo "  mp1 entries:    $(wc -l < "$TESTDIR/norm_mp1.txt")"
echo "  smp entries:    $(wc -l < "$TESTDIR/norm_smp.txt")"
echo "  mp2 t1 entries: $(wc -l < "$TESTDIR/norm_mp2_t1.txt")"
echo "  mp2 t4 entries: $(wc -l < "$TESTDIR/norm_mp2_t4.txt")"
echo ""

# --------------------------------------------------------------------------
# Step 6: Comparisons
# --------------------------------------------------------------------------
echo "=== Step 6: Comparison ==="

python3 - "$TESTDIR" <<'PYEOF'
import sys, os

def load_norm(path):
    """Load norm file: (chr,pos,ref,alt,sample) -> (rf,af,rr,ar)."""
    d = {}
    pos_d = {}  # (chr,pos,sample) -> list of (ref,alt,val)
    with open(path) as f:
        for line in f:
            p = line.strip().split()
            if len(p) >= 9:
                key = (p[0], p[1], p[2], p[3], p[4])
                val = (int(p[5]), int(p[6]), int(p[7]), int(p[8]))
                d[key] = val
                pos_d.setdefault((p[0], p[1], p[4]), []).append((p[2], p[3], val))
    return d, pos_d

def alleles_match_at_pos(pos_a, pos_b):
    counts_a = sorted(v for _, _, v in pos_a)
    counts_b = sorted(v for _, _, v in pos_b)
    return counts_a == counts_b

def compare(a, a_pos, b, b_pos, name_a, name_b):
    only_a_keys = set(a) - set(b)
    only_b_keys = set(b) - set(a)
    shared = set(a) & set(b)
    match = sum(1 for k in shared if a[k] == b[k])

    n_base_diff = set()
    for k in only_a_keys:
        pos_key = (k[0], k[1], k[4])  # chr, pos, sample
        if pos_key in b_pos:
            if alleles_match_at_pos(a_pos.get(pos_key, []), b_pos.get(pos_key, [])):
                n_base_diff.add(f'{k[0]}:{k[1]}:{k[4]}')
    truly_only_a = {k for k in only_a_keys if f'{k[0]}:{k[1]}:{k[4]}' not in n_base_diff}
    truly_only_b = {k for k in only_b_keys if f'{k[0]}:{k[1]}:{k[4]}' not in n_base_diff}

    print(f"  Only in {name_a}: {len(truly_only_a)}"
          + (f"  (+ {len(only_a_keys)-len(truly_only_a)} N-base allele diffs, counts match)" if len(only_a_keys) > len(truly_only_a) else ""))
    print(f"  Only in {name_b}: {len(truly_only_b)}"
          + (f"  (+ {len(only_b_keys)-len(truly_only_b)} N-base allele diffs, counts match)" if len(only_b_keys) > len(truly_only_b) else ""))
    print(f"  Shared entries: {len(shared)}")
    print(f"  Exact count match: {match} / {len(shared)}")
    if truly_only_a:
        print(f"  Entries only in {name_a} (first 5):")
        for k in sorted(truly_only_a)[:5]:
            print(f"    {k[0]} {k[1]} {k[2]} {k[3]} sample={k[4]} {a[k]}")
    if truly_only_b:
        print(f"  Entries only in {name_b} (first 5):")
        for k in sorted(truly_only_b)[:5]:
            print(f"    {k[0]} {k[1]} {k[2]} {k[3]} sample={k[4]} {b[k]}")
    if match < len(shared):
        print(f"  Count mismatches (first 5):")
        n = 0
        for k in sorted(shared):
            if a[k] != b[k]:
                print(f"    {k[0]} {k[1]} {k[2]} {k[3]} sample={k[4]}  {name_a}:{a[k]}  {name_b}:{b[k]}")
                n += 1
                if n >= 5:
                    break
    return len(truly_only_a), len(truly_only_b), match, len(shared)

testdir = sys.argv[1]
mp1,    mp1p    = load_norm(os.path.join(testdir, 'norm_mp1.txt'))
mp2_t1, mp2_t1p = load_norm(os.path.join(testdir, 'norm_mp2_t1.txt'))
mp2_t4, mp2_t4p = load_norm(os.path.join(testdir, 'norm_mp2_t4.txt'))
smp,    smpp    = load_norm(os.path.join(testdir, 'norm_smp.txt'))

print("--- mp2 t=1 vs mp2 t=4 (threading consistency) ---")
ot1, ot4, mt_match, mt_shared = compare(mp2_t1, mp2_t1p, mp2_t4, mp2_t4p, "mp2_t1", "mp2_t4")
print()
print("--- mp1 vs mp2 (t=1) ---")
o1, o2, match, shared = compare(mp1, mp1p, mp2_t1, mp2_t1p, "mp1", "mp2_t1")
print()
print("--- mp1 vs samtools mpileup ---")
compare(mp1, mp1p, smp, smpp, "mp1", "smp")
print()
print("--- mp2 (t=1) vs samtools mpileup ---")
compare(mp2_t1, mp2_t1p, smp, smpp, "mp2_t1", "smp")
print()
print("Full output files:", testdir)
print()
print("=== SUMMARY ===")
thread_ok = (ot1 == 0 and ot4 == 0 and mt_match == mt_shared)
mp_ok     = (o1  == 0 and o2  == 0 and match    == shared)
if thread_ok:
    print("PASS (threading): mp2 t=1 and mp2 t=4 produce identical counts.")
else:
    print(f"FAIL (threading): mp2 t=1 and mp2 t=4 differ (only_t1={ot1}, only_t4={ot4}, match={mt_match}/{mt_shared}).")
if mp_ok:
    print("PASS (cross-tool): minipileup and minipileup2 produce identical counts.")
else:
    print(f"FAIL (cross-tool): minipileup and minipileup2 differ (only_mp1={o1}, only_mp2={o2}, match={match}/{shared}).")
PYEOF

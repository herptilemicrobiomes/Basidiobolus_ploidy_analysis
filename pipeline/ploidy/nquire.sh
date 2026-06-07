#!/bin/bash -l
#SBATCH -p short -c 8 --mem 24gb --out logs/nquire.%a.log

# nQuire ploidy estimation
# Runs whole-genome and per-scaffold analysis for one BAM per array task.
# Usage: sbatch --array=1-N pipeline/ploidy/nquire.sh

module load nQuire
module load samtools

CPU=${SLURM_CPUS_ON_NODE:-2}

N=${SLURM_ARRAY_TASK_ID:-$1}
if [ -z "$N" ]; then
    echo "ERROR: no array task ID or cmdline argument"
    exit 1
fi

# ---- Paths ----------------------------------------------------------------
ALN=aln
GENEDIR=bed
RESULT=results/nquire
SCRATCH_BASE=working
#${SCRATCH:-/tmp}

mkdir -p "$RESULT"

# ---- Select BAM by array index --------------------------------------------
mapfile -t BAMFILES < <(ls "$ALN"/*.bam 2>/dev/null | sort)
if [ ${#BAMFILES[@]} -eq 0 ]; then
    echo "ERROR: no BAM files found in $ALN/"
    exit 1
fi

BAM=${BAMFILES[$((N-1))]}
if [ -z "$BAM" ] || [ ! -f "$BAM" ]; then
    echo "ERROR: no BAM at index $N (found ${#BAMFILES[@]} files)"
    exit 1
fi

STRAIN=$(basename "$BAM" .bam)
WORKDIR="$SCRATCH_BASE/nquire_${STRAIN}_$$"
mkdir -p "$WORKDIR"

echo "N=$N  BAM=$BAM  STRAIN=$STRAIN"

# ---- Build per-scaffold BED from BAM header -------------------------------
# Format required by nQuire: Chr Start End Name
BED="$WORKDIR/${STRAIN}.scaffolds.bed"
samtools view -H "$BAM" \
    | awk '/^@SQ/ {
        split($2, sn, ":"); split($3, ln, ":");
        print sn[2] "\t0\t" ln[2] "\t" sn[2]
      }' > "$BED"

NSCAFFOLDS=$(wc -l < "$BED")
echo "Scaffolds in BED: $NSCAFFOLDS"

if [ "$NSCAFFOLDS" -eq 0 ]; then
    echo "ERROR: no @SQ lines found in BAM header for $BAM"
    rm -rf "$WORKDIR"
    exit 1
fi

# ==========================================================================
# Whole-genome analysis
# ==========================================================================
echo "--- Whole genome ---"
WG_PREFIX="$WORKDIR/${STRAIN}.wg"

nQuire create -b "$BAM" -o "$WG_PREFIX" -c 5 -q 20
nQuire denoise -o "${WG_PREFIX}.denoised" "${WG_PREFIX}.bin"

# lrdmodel outputs: header line then one data line per input file
WG_TSV="$RESULT/${STRAIN}.nquire.tsv"
{
    # header: replace leading "name" column with "strain"
    nQuire lrdmodel -t "$CPU" "${WG_PREFIX}.denoised.bin" | awk -v s="$STRAIN" '
        NR==1 { sub(/^[^\t]+/, "strain"); print; next }
              { sub(/^[^\t]+/, s);        print }
    '
} > "$WG_TSV"
echo "Whole-genome result -> $WG_TSV"

# ==========================================================================
# Per-scaffold analysis
# ==========================================================================
echo "--- Per scaffold ---"
PERCHROM_TSV="$RESULT/${STRAIN}.perchrom.nquire.tsv"
HEADER_WRITTEN=0

mkdir -p "$WORKDIR/${STRAIN}/perchrom" 
mkdir -p "$WORKDIR/${STRAIN}/perchrom.denoised" 

#> "$PERCHROM_TSV"   # clear/create

nQuire create -b "$BAM" -r "$BED" -o "$WORKDIR/${STRAIN}/perchrom" -c 4 -q 20
for file in $(ls "${WORKDIR}/${STRAIN}/perchrom/*.bin")
do
	PREF=$(basename $file .bin)
	nQuire denoise -o "$WORKDIR/${STRAIN}/perchrom.denoised" $file
	#"$$WORKDIR/${STRAIN}.perchrom.bin"
	nQuire lrdmodel -t "$CPU" "$WORKDIR/${STRAIN}/perchrom.denoised/${PREF}.denoised.bin" 
done > $PERCHROM_TSV

# =====
# Per gene analysis
# ====
#
# echo "--- Per gene ---"
PERGENE_TSV="$RESULT/${STRAIN}.pergene.nquire.tsv"
HEADER_WRITTEN=0
GENE=$GENEDIR/${STRAIN}.bed
mkdir -p "${WORKDIR}/${STRAIN}/pergene"
nquire create -b "$BAM" -o "$WORKDIR/${STRAIN}/pergene" -r $GENE
for file in $(ls "${WORKDIR}/${STRAIN}/pergene/*.bin")
do
	nQuire denoise -o "${WORKDIR}/${STRAIN}/pergene.denoised" $file
done

#while IFS=$'\t' read -r CHROM START END NAME; do
#    SC_PREFIX="$WORKDIR/${STRAIN}.${NAME}"
#    SC_BED="$WORKDIR/${NAME}.bed"
#    printf '%s\t%s\t%s\t%s\n' "$CHROM" "$START" "$END" "$NAME" > "$SC_BED"
#
#    
 #   if [ ! -f "${SC_PREFIX}.bin" ]; then
 #       echo "  WARN: nQuire create produced no output for $NAME, skipping"
 #       continue
 #   fi
#
#    nQuire denoise -o "${SC_PREFIX}.denoised" "${SC_PREFIX}.bin"
#
#    nQuire lrdmodel -t "$CPU" "${SC_PREFIX}.denoised.bin" | awk -v s="$STRAIN" -v sc="$NAME" -v h="$HEADER_WRITTEN" '
#        NR==1 {
#            if (h == 0) { sub(/^[^\t]+/, "strain\tscaffold"); print }
#            next
#        }
#        { sub(/^[^\t]+/, s "\t" sc); print }
#    ' >> "$PERCHROM_TSV"
#    HEADER_WRITTEN=1
#
#done < "$BED"

echo "Per-scaffold result -> $PERCHROM_TSV"

# ---- Cleanup scratch ------------------------------------------------------
#rm -rf "$WORKDIR"
echo "Done: $STRAIN"

#!/bin/bash -l
#SBATCH -p short -c 8 --mem 24gb --out logs/nquire.%a.log -a 1

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
SCAFFOLDBED=scaffold_bed
RESULT=results/nquire2
SCRATCH_BASE=${SCRATCH:-/tmp}
mkdir -p "$RESULT" $SCAFFOLDBED

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
WORKDIR="$SCRATCH_BASE/nquire_${STRAIN}"
mkdir -p "$WORKDIR"

echo "N=$N  BAM=$BAM  STRAIN=$STRAIN"

# ---- Build per-scaffold BED from BAM header -------------------------------
# Format required by nQuire: Chr Start End Name
BED="$SCAFFOLDBED/${STRAIN}.scaffolds.bed"
if [ ! -f $BED ]; then
    samtools view -H "$BAM" \
	| awk '/^@SQ/ {
          split($2, sn, ":"); split($3, ln, ":");
        print sn[2] "\t0\t" ln[2] "\t" sn[2]
      	}' > "$BED"
fi
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
WG_TSV="$RESULT/${STRAIN}.nquire.tsv"
if [ ! -s $WG_TSV ]; then
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
fi

# ==========================================================================
# Per-scaffold analysis
# ==========================================================================
echo "--- Per scaffold ---"
PERCHROM_TSV="$RESULT/${STRAIN}.perchrom.nquire.tsv"
HEADER_WRITTEN=0

mkdir -p "$WORKDIR/perchrom" 
mkdir -p "$WORKDIR/perchrom.denoised" 

#> "$PERCHROM_TSV"   # clear/create
if [ ! -s $PERCHROM_TSV ]; then
    echo "RUNNING nQuire create -b \"$BAM\" -r \"$BED\" -o \"$WORKDIR/perchrom/${STRAIN}\" -c 4 -q 20"
    nQuire create -b "$BAM" -r "$BED" -o "$WORKDIR/perchrom/${STRAIN}" -c 4 -q 20
    for file in $(ls "${WORKDIR}/perchrom/${STRAIN}-"*.bin)
    do
	PREF=$(basename $file .bin)
    # -o takes a prefix; nQuire appends .bin → pass PREF.denoised to get PREF.denoised.bin
	nQuire denoise -o "$WORKDIR/perchrom.denoised/${PREF}" $file
	nQuire lrdmodel -t "$CPU" "$WORKDIR/perchrom.denoised/${PREF}.bin" 1> $WORKDIR/perchrom.denoised/${PREF}.model.out 2> /dev/null
    done
    head -n 1 $(ls $WORKDIR/perchrom.denoised/*.model.out | head -n 1) | perl -p -e 's/file/scaffold/' > $PERCHROM_TSV
    grep -vh 'free' $WORKDIR/perchrom.denoised/*.model.out | perl -p -e 's/^\S+\-(\S+)\.bin/$1/' | sort -t_ -k 2,2n >> $PERCHROM_TSV
fi

echo "Per-scaffold result -> $PERCHROM_TSV"
# =====
# Per gene analysis
# ====
#
# echo "--- Per gene ---"
PERGENE_TSV="$RESULT/${STRAIN}.pergene.nquire.tsv"
HEADER_WRITTEN=0
GENE=$GENEDIR/${STRAIN}.bed
mkdir -p "${WORKDIR}/pergene" "${WORKDIR}/pergene.denoised"
if [ ! -s $PERGENE_TSV ]; then
    nQuire create -b "$BAM" -o "$WORKDIR/pergene/" -r $GENE -c 4 -q 20

    bin_count=$(ls "${WORKDIR}/pergene"/*.bin 2>/dev/null | wc -l)
    echo "[nQuire pergene] $STRAIN: $bin_count gene bins" >&2

    # Phase 1: denoise all bins with retry
    for file in "${WORKDIR}/pergene"/*.bin; do
	[ -f "$file" ] || continue
	PREF=$(basename "$file" .bin)
	attempt=0
	while [ $attempt -lt 3 ]; do
            # -o takes a prefix; nQuire appends .bin → pass PREF.denoised to get PREF.denoised.bin
            nQuire denoise -o "${WORKDIR}/pergene.denoised/${PREF}" "$file" && break
            attempt=$((attempt + 1))
            echo "[nQuire pergene] retry $attempt/3: denoise failed for $PREF" >&2
	done
	[ $attempt -eq 3 ] && echo "[nQuire pergene] WARN: denoise failed for $PREF" >&2
    done

    # Phase 2: lrdmodel all denoised bins, header printed once
    HEADER_WRITTEN=0
    for denoised in "${WORKDIR}/pergene.denoised"/*.bin; do
	[ -f "$denoised" ] || continue
	PREF=$(basename "$denoised" .bin)
	nQuire lrdmodel -t "$CPU" "$WORKDIR/pergene.denoised/${PREF}.bin" 1> $WORKDIR/pergene.denoised/${PREF}.model.out 2> /dev/null
    done
    head -n 1 $(ls $WORKDIR/pergene.denoised/*.model.out | head -n 1) | perl -p -e 's/file/gene/' > $PERGENE_TSV
    grep -vh 'free' $WORKDIR/pergene.denoised/*.model.out | perl -p -e 's/^\S+\-(\S+)\.bin/$1/' | sort -t_ -k 2,2n >> $PERGENE_TSV
fi


echo "Per-gene result -> $PERGENE_TSV"

# ---- Cleanup scratch ------------------------------------------------------
#rm -rf "$WORKDIR"
echo "Done: $STRAIN"

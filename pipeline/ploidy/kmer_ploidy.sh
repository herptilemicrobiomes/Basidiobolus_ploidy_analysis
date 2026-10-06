#!/bin/bash -l
#SBATCH -p epyc -c 16 --mem 64gb -t 1-00:00:00 --out logs/kmer_ploidy.%a.log -a 1-8

# Reference-free ploidy estimation from short-read k-mers.
#   FastK      count k-mers in the paired reads
#   Histex     k-mer count histogram (GenomeScope format)
#   GenomeScope2  fit p = 2, 3, 4 models to the histogram
#   Smudgeplot hetmer pairs -> smudgeplot (AB / AAB / AAAB / AABB ...)
#
# No reference or mapping is involved, so collapsed paralogs, retained
# haplotigs and mapping bias that affect nQuire do not affect this test in the
# same way. Collapsed-paralog k-mer pairs show up as smudges at 2x coverage
# rather than as AAAB at 1x.
#
# One array task per short-read row of samples.csv (rows with a Reverse file).
# Usage:  sbatch pipeline/ploidy/kmer_ploidy.sh          (8 short-read samples)
#         bash pipeline/ploidy/kmer_ploidy.sh 3           (run sample 3 locally)
#         K=31 sbatch pipeline/ploidy/kmer_ploidy.sh      (change k)
#
# Commands follow the Smudgeplot >=0.3 (FastK-based) README; not yet run on
# the cluster. Module names below are a guess; check with `module avail`.
# CBS931.73 short reads are ~1000x; if GenomeScope fails to converge, try
# downsampling those reads to ~100x first (e.g. seqtk sample).

module load fastk
module load genomescope
module load smudgeplot

CPU=${SLURM_CPUS_ON_NODE:-4}
MEM_GB=${FASTK_MEM_GB:-48}
K=${K:-21}

N=${SLURM_ARRAY_TASK_ID:-$1}
if [ -z "$N" ]; then
    echo "ERROR: no array task ID or cmdline argument"
    exit 1
fi

SAMPLES=samples.csv
IN=reads
RESULT=results/kmer_ploidy
SCRATCH_BASE=${SCRATCH:-kmer_scratch}

# short-read rows only: Reverse column is non-empty
ROW=$(tail -n +2 $SAMPLES | awk -F, '$3 != ""' | sed -n ${N}p)
if [ -z "$ROW" ]; then
    echo "ERROR: no short-read sample at index $N"
    exit 1
fi
IFS=, read ID FWD REV <<< "$ROW"
echo "N=$N ID=$ID FWD=$FWD REV=$REV K=$K"

OUT=$RESULT/$ID
WORK=$SCRATCH_BASE/kmer_$ID
mkdir -p $OUT $WORK

# ---- 1. k-mer counting ---------------------------------------------------
TABLE=$WORK/${ID}.k$K
if [ ! -s ${TABLE}.ktab ]; then
    # -t4: keep k-mers seen >= 4 times in the table (smudgeplot input);
    # the .hist still contains the full spectrum.
    FastK -v -t4 -k$K -M$MEM_GB -T$CPU -P$WORK $IN/$FWD $IN/$REV -N$TABLE
fi

# ---- 2. histogram --------------------------------------------------------
HIST=$OUT/${ID}.k$K.hist
if [ ! -s $HIST ]; then
    # wide range so very deep samples (CBS931.73 is ~1000x) are not truncated
    Histex -G -h1:32767 $TABLE > $HIST
fi

# max k-mer coverage for GenomeScope: highest count with a non-zero bin
MAXCOV=$(awk '$2 > 0 {m = $1} END {print m}' $HIST)

# ---- 3. GenomeScope2, p = 2, 3, 4 ----------------------------------------
for P in 2 3 4; do
    GS=$OUT/genomescope_p$P
    if [ ! -s $GS/summary.txt ]; then
        genomescope2 -i $HIST -o $GS -k $K -p $P -n ${ID}_p$P -m $MAXCOV
    fi
done

# one-line-per-model summary: model fit and heterozygosity rows
{
    printf "sample\tp\tline\n"
    for P in 2 3 4; do
        S=$OUT/genomescope_p$P/summary.txt
        [ -s $S ] || continue
        grep -E "Model Fit|Homozygous|Heterozygous|^a[a-z]+ |Haploid Length" $S \
            | sed "s/^/${ID}\t$P\t/"
    done
} > $OUT/genomescope_summary.tsv

# ---- 4. Smudgeplot -------------------------------------------------------
# lower count cutoff L = first valley of the histogram after the error peak
L=$(awk 'NR > 1 && $2 > prev && $1 > 4 {print $1 - 1; exit} {prev = $2}' $HIST)
L=${L:-12}
echo "Smudgeplot lower cutoff L=$L" | tee $OUT/smudgeplot_L.txt

PAIRS=$WORK/${ID}.k$K.pairs
if [ ! -s ${PAIRS}_text.smu ]; then
    smudgeplot.py hetmers -L $L -t $CPU -o $PAIRS --verbose $TABLE
fi
if [ ! -s $OUT/${ID}_smudgeplot.pdf ] && [ ! -s $OUT/${ID}_smudgeplot_log10.pdf ]; then
    smudgeplot.py all -o $OUT/${ID} ${PAIRS}_text.smu
fi

echo "Done: $ID -> $OUT"
# FastK tables are large; remove once results look right
# Fastrm $TABLE

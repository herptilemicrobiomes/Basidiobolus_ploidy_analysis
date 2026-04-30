#!/bin/bash -l
#SBATCH -p short -c 96 --mem 96gb --out logs/bwamem.%a.log

: "${SCRATCH:?SCRATCH environment variable is not set}"

module load samtools

CPU=${SLURM_CPUS_ON_NODE:-2}
MOSTHREADS=$(( CPU < 4 ? CPU : 4 ))

N=${SLURM_ARRAY_TASK_ID:-$1}
if [ -z "$N" ]; then
    echo "cannot find a cmdline option or array/-a option"
    exit 1
fi
echo "N is $N"
SAMPLES=samples.csv
OUT=aln
DB=genome
IN=reads
BEDIN=bed
RESULT=results/mosdepth

mkdir -p $OUT $RESULT
tail -n +2 $SAMPLES | sed -n ${N}p | while IFS=, read ID FWD REV
do

    if [[ $FWD == *pacbio* ]]; then
        BAM=$OUT/${ID}_pb.bam
        PRESET=map-pb
    elif [[ $FWD == *ONT* ]]; then
        BAM=$OUT/${ID}_ont.bam
        PRESET=map-ont
    else
        BAM=$OUT/${ID}.bam
        PRESET=""
    fi
    echo "BAM is $BAM"
    if [ ! -s $BAM ]; then
        if [ -n "$PRESET" ]; then
            module load minimap2 samtools
            TYPE=${PRESET#map-}  # 'pb' or 'ont'
            minimap2 -t $CPU -ax $PRESET $DB/$ID.masked.fasta $IN/$FWD > $SCRATCH/${ID}.sam
            samtools flagstat $SCRATCH/${ID}.sam > $OUT/${ID}_${TYPE}.flagstat.txt
        else
            module load bwa-mem2 samtools
            # -M -k 15 -c 1000 old bwa-mem params don't seem to make much difference
            bwa-mem2 mem -o $SCRATCH/${ID}.sam -t $CPU $DB/${ID}.masked.fasta $IN/$FWD $IN/$REV
            samtools flagstat $SCRATCH/${ID}.sam > $OUT/${ID}.self.flagstat.txt
        fi
        samtools view --threads $MOSTHREADS -OBAM -F 12 -o $SCRATCH/${ID}.bam $SCRATCH/${ID}.sam
        samtools sort -OBAM --threads $CPU -o $BAM $SCRATCH/${ID}.bam
    fi

    if [ ! -f ${BAM}.bai ]; then
        samtools index $BAM
    fi
    module load mosdepth
    export MOSDEPTH_Q0=NO_COVERAGE
    export MOSDEPTH_Q1=LOW_COVERAGE
    export MOSDEPTH_Q2=CALLABLE
    export MOSDEPTH_Q3=HIGH_COVERAGE
    export MOSDEPTH_Q4=VERY_HIGH_COVERAGE
    PREFIX=$RESULT/${ID}${PRESET:+_${PRESET#map-}}
    echo "pref=$PREFIX bam=$BAM bed=$BEDIN/$ID.bed"
    mosdepth -n --quantize 0:1:4:100:200: -t $MOSTHREADS -b $BEDIN/$ID.bed $PREFIX $BAM
done

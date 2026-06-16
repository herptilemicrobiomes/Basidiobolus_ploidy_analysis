#!/usr/bin/env nextflow
/*
 * nquire.nf — DSL2 nQuire ploidy-estimation workflow.
 *
 * Three analyses per sample:
 *   1. Whole-genome   → results/nquire/<base>.nquire.tsv
 *   2. Per-scaffold   → results/nquire/<base>.perchrom.nquire.tsv
 *   3. Per-gene       → results/nquire/<base>.pergene.nquire.tsv
 *
 * Naming follows pipeline/ploidy/nquire.sh:
 *   base   = BAM basename (e.g. STP1710.7_ont)        → used for output filenames
 *   strain = base with trailing _pb/_ont stripped     → used for gene-BED lookup,
 *            scaffold/gene nQuire prefixes, and the strain column value
 *
 * Input: BAM files (with .bai indices) discovered from params.alndir.
 * Optional per-strain gene BED files from params.beddir (<strain>.bed).
 *
 * Usage:
 *   nextflow run pipeline/ploidy/nquire.nf -c nextflow.config -profile nquire -resume
 */

nextflow.enable.dsl = 2

params.alndir = "${launchDir}/aln"
params.beddir = "${launchDir}/bed"
params.outdir = "${launchDir}/results/nquire"
params.n_test = 0           // >0: limit to the first N samples (for testing)
params.min_cov_wg   = 5    // -c flag for whole-genome nQuire create
params.min_cov_sc   = 4    // -c flag for per-scaffold nQuire create
params.min_cov_gene = 4    // -c flag for per-gene nQuire create
params.min_mapq     = 20   // -q flag for all nQuire create calls

// ── Processes ────────────────────────────────────────────────────────────────

/*
 * Extract a scaffold BED from the BAM header (samtools view -H).
 * Output format: Chr\t0\tLen\tChr   (required by nQuire -r)
 */
process BAM_TO_SCAFFOLD_BED {
    tag { base }
    label 'quick'
    publishDir "${launchDir}/scaffold_bed", mode: 'copy', pattern: "*.bed"

    input:
    tuple val(base), val(strain), path(bam), path(bai)

    output:
    tuple val(base), val(strain), path(bam), path(bai), path("${base}.scaffolds.bed")

    script:
    """
    module load samtools
    samtools view -H ${bam} \
        | awk '/^@SQ/ {
            split(\$2, sn, ":"); split(\$3, ln, ":");
            print sn[2] "\\t0\\t" ln[2] "\\t" sn[2]
          }' > ${base}.scaffolds.bed

    if [ ! -s ${base}.scaffolds.bed ]; then
        echo "ERROR: no @SQ lines found in BAM header for ${bam}" >&2
        exit 1
    fi
    """

    stub:
    """
    printf 'scaffold1\\t0\\t1000000\\tscaffold1\\n' > ${base}.scaffolds.bed
    """
}

/*
 * Whole-genome nQuire: create → denoise → lrdmodel.
 * All three steps run in one process to avoid staging intermediate .bin files.
 */
process NQUIRE_WHOLE_GENOME {
    tag { base }
    label 'nquire'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(base), val(strain), path(bam), path(bai)

    output:
    tuple val(base), path("${base}.nquire.tsv")

    script:
    """
    module load nQuire
    nQuire create -b ${bam} -o ${base}.wg -c ${params.min_cov_wg} -q ${params.min_mapq}
    nQuire denoise -o ${base}.wg.denoised ${base}.wg.bin

    # lrdmodel outputs: header line then one data line per input file
    # header: replace leading "name" column with "strain"; data: substitute strain
    nQuire lrdmodel -t ${task.cpus} ${base}.wg.denoised.bin \
        | awk -v s="${strain}" '
            NR==1 { sub(/^[^\\t]+/, "strain"); print; next }
                  { sub(/^[^\\t]+/, s);        print }
        ' > ${base}.nquire.tsv
    """

    stub:
    """
    printf 'strain\\tdip\\trip\\ttet\\tbest\\n' > ${base}.nquire.tsv
    printf '${strain}\\t0\\t0\\t0\\tdip\\n'   >> ${base}.nquire.tsv
    """
}

/*
 * Per-scaffold nQuire: create across all scaffolds → denoise each bin →
 * lrdmodel each bin to individual .model.out → combine into one TSV.
 */
process NQUIRE_PER_SCAFFOLD {
    tag { base }
    label 'nquire'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(base), val(strain), path(bam), path(bai), path(scaffold_bed)

    output:
    tuple val(base), path("${base}.perchrom.nquire.tsv")

    script:
    """
    mkdir -p perchrom perchrom_denoised
    module load nQuire
    nQuire create -b ${bam} -r ${scaffold_bed} \
        -o perchrom/${strain} -c ${params.min_cov_sc} -q ${params.min_mapq}

    for bin in \$(find perchrom -name "${strain}-*.bin"); do
        [ -f "\$bin" ] || continue
        pref=\$(basename "\$bin" .bin)
        # -o takes a prefix; nQuire appends .bin → pass pref.denoised to get pref.denoised.bin
        nQuire denoise -o "perchrom_denoised/\${pref}" "\$bin"
        nQuire lrdmodel -t ${task.cpus} "perchrom_denoised/\${pref}.bin" \
            1> "perchrom_denoised/\${pref}.model.out" 2>/dev/null
    done

    head -n 1 \$(find perchrom_denoised -name "*.model.out" | head -n 1) \
        | perl -p -e 's/file/scaffold/' > ${base}.perchrom.nquire.tsv
    find perchrom_denoised -name "*.model.out" | xargs -r grep -vh 'free' \
        | perl -p -e 's/^\\S+-(\\S+)\\.bin/\$1/' \
        | sort -t_ -k 2,2n >> ${base}.perchrom.nquire.tsv

    touch ${base}.perchrom.nquire.tsv
    """

    stub:
    """
    printf 'scaffold\\tdip\\trip\\ttet\\tbest\\n' > ${base}.perchrom.nquire.tsv
    """
}

/*
 * Per-gene nQuire: create across gene BED → denoise each bin (with retry) →
 * lrdmodel each bin to individual .model.out → combine into one TSV.
 * Skipped automatically when no gene BED is found for a strain.
 */
process NQUIRE_PER_GENE {
    tag { base }
    label 'nquire'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(base), val(strain), path(bam), path(bai), path(gene_bed)

    output:
    tuple val(base), path("${base}.pergene.nquire.tsv")

    script:
    """
    mkdir -p pergene pergene_denoised
    module load nQuire
    nQuire create -b ${bam} -r ${gene_bed} \
        -o pergene/${strain} -c ${params.min_cov_gene} -q ${params.min_mapq}

    bin_count=\$(find pergene -name "${strain}-*.bin" 2>/dev/null | wc -l)
    echo "[nQuire pergene] ${strain}: \${bin_count} gene bins" >&2

    # Phase 1: denoise all bins with retry
    for bin in \$(find pergene -name "*.bin"); do
        [ -f "\$bin" ] || continue
        pref=\$(basename "\$bin" .bin)
        attempt=0
        while [ \$attempt -lt 3 ]; do
            # -o takes a prefix; nQuire appends .bin → pass pref.denoised to get pref.denoised.bin
            nQuire denoise -o "pergene_denoised/\${pref}" "\$bin" && break
            attempt=\$((attempt + 1))
            echo "[nQuire pergene] retry \${attempt}/3: denoise failed for \${pref}" >&2
        done
        if [ \$attempt -eq 3 ]; then
            echo "[nQuire pergene] WARN: denoise failed after 3 attempts for \${pref}, skipping" >&2
        fi
    done

    # Phase 2: lrdmodel all denoised bins to individual .model.out files
    for bin in \$(find pergene_denoised -name "*.bin"); do
        [ -f "\$bin" ] || continue
        pref=\$(basename "\$bin" .bin)
        nQuire lrdmodel -t ${task.cpus} "\$bin" \
            1> "pergene_denoised/\${pref}.model.out" 2>/dev/null
    done

    head -n 1 \$(find pergene_denoised -name "*.model.out" | head -n 1) \
        | perl -p -e 's/file/gene/' > ${base}.pergene.nquire.tsv
    find pergene_denoised -name "*.model.out" | xargs -r grep -vh 'free' \
        | perl -p -e 's/^\\S+-(\\S+)\\.bin/\$1/' \
        | sort -t_ -k 2,2n >> ${base}.pergene.nquire.tsv

    touch ${base}.pergene.nquire.tsv
    """

    stub:
    """
    printf 'gene\\tdip\\trip\\ttet\\tbest\\n' > ${base}.pergene.nquire.tsv
    """
}

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow {

    // Discover BAMs; derive base (full basename) and strain (base minus _pb/_ont).
    // .take(-1) means "all"; params.n_test > 0 limits to the first N for quick tests.
    bam_ch = Channel
        .fromPath("${params.alndir}/*.bam", checkIfExists: true)
        .map { bam ->
            def base   = bam.baseName
            def strain = base.replaceAll(/_(pb|ont)$/, '')
            def bai    = file("${bam}.bai", checkIfExists: false)
            if ( !bai.exists() ) bai = file("${bam.parent}/${base}.bai", checkIfExists: false)
            tuple(base, strain, bam, bai)
        }
        .take( (params.n_test as int) > 0 ? params.n_test as int : -1 )

    // 1. Whole-genome analysis (no scaffold BED needed)
    NQUIRE_WHOLE_GENOME(bam_ch)

    // 2. Per-scaffold: derive scaffold BED from each BAM header
    BAM_TO_SCAFFOLD_BED(bam_ch)
    NQUIRE_PER_SCAFFOLD(BAM_TO_SCAFFOLD_BED.out)

    // 3. Per-gene: join BAMs with the matching gene BED by strain (bed/<strain>.bed).
    //    Inner join → strains with no BED are skipped automatically.
    gene_bed_ch = Channel
        .fromPath("${params.beddir}/*.bed", checkIfExists: false)
        .map { bed -> tuple(bed.baseName, bed) }   // keyed by strain

    bam_with_gene_ch = bam_ch
        .map { base, strain, bam, bai -> tuple(strain, base, bam, bai) }   // key by strain
        .join(gene_bed_ch, by: 0, remainder: false)
        .map { strain, base, bam, bai, bed -> tuple(base, strain, bam, bai, bed) }

    NQUIRE_PER_GENE(bam_with_gene_ch)
}

#!/usr/bin/env nextflow
/*
 * nquire.nf — DSL2 nQuire ploidy-estimation workflow.
 *
 * Three analyses per sample:
 *   1. Whole-genome   → results/nquire/<strain>.nquire.tsv
 *   2. Per-scaffold   → results/nquire/<strain>.perchrom.nquire.tsv
 *   3. Per-gene       → results/nquire/<strain>.pergene.nquire.tsv
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
    tag { strain }
    label 'quick'

    input:
    tuple val(strain), path(bam), path(bai)

    output:
    tuple val(strain), path(bam), path(bai), path("${strain}.scaffolds.bed")

    script:
    """
    module load samtools
    samtools view -H ${bam} \
        | awk '/^@SQ/ {
            split(\$2, sn, ":"); split(\$3, ln, ":");
            print sn[2] "\\t0\\t" ln[2] "\\t" sn[2]
          }' > ${strain}.scaffolds.bed

    if [ ! -s ${strain}.scaffolds.bed ]; then
        echo "ERROR: no @SQ lines found in BAM header for ${bam}" >&2
        exit 1
    fi
    """

    stub:
    """
    printf 'scaffold1\\t0\\t1000000\\tscaffold1\\n' > ${strain}.scaffolds.bed
    """
}

/*
 * Whole-genome nQuire: create → denoise → lrdmodel.
 * All three steps run in one process to avoid staging intermediate .bin files.
 */
process NQUIRE_WHOLE_GENOME {
    tag { strain }
    label 'nquire'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(strain), path(bam), path(bai)

    output:
    tuple val(strain), path("${strain}.nquire.tsv")

    script:
    """
    module load nQuire
    nQuire create -b ${bam} -o ${strain}.wg -c ${params.min_cov_wg} -q ${params.min_mapq}
    nQuire denoise -o ${strain}.wg.denoised ${strain}.wg.bin

    nQuire lrdmodel -t ${task.cpus} ${strain}.wg.denoised.bin \
        | awk -v s="${strain}" '
            NR==1 { sub(/^[^\\t]+/, "strain"); print; next }
                  { sub(/^[^\\t]+/, s);        print }
        ' > ${strain}.nquire.tsv
    """

    stub:
    """
    printf 'strain\\tdip\\trip\\ttet\\tbest\\n' > ${strain}.nquire.tsv
    printf '${strain}\\t0\\t0\\t0\\tdip\\n'   >> ${strain}.nquire.tsv
    """
}

/*
 * Per-scaffold nQuire: create across all scaffolds → denoise each bin →
 * lrdmodel each bin → concatenate into one TSV.
 */
process NQUIRE_PER_SCAFFOLD {
    tag { strain }
    label 'nquire'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(strain), path(bam), path(bai), path(scaffold_bed)

    output:
    tuple val(strain), path("${strain}.perchrom.nquire.tsv")

    script:
    """
    mkdir -p perchrom perchrom_denoised
    module load nQuire
    nQuire create -b ${bam} -r ${scaffold_bed} \
        -o perchrom -c ${params.min_cov_sc} -q ${params.min_mapq}

    header_written=0
    for bin in perchrom/*.bin; do
        [ -f "\$bin" ] || continue
        pref=\$(basename "\$bin" .bin)
        nQuire denoise -o perchrom_denoised "\$bin"
        denoised="perchrom_denoised/\${pref}.denoised.bin"
        [ -f "\$denoised" ] || continue

        if [ \$header_written -eq 0 ]; then
            nQuire lrdmodel -t ${task.cpus} "\$denoised" \
                | awk -v s="${strain}" -v sc="\$pref" '
                    NR==1 { sub(/^[^\\t]+/, "strain\\tscaffold"); print; next }
                          { sub(/^[^\\t]+/, s "\\t" sc);          print }
                ' >> ${strain}.perchrom.nquire.tsv
            header_written=1
        else
            nQuire lrdmodel -t ${task.cpus} "\$denoised" \
                | awk -v s="${strain}" -v sc="\$pref" '
                    NR==1 { next }
                          { sub(/^[^\\t]+/, s "\\t" sc); print }
                ' >> ${strain}.perchrom.nquire.tsv
        fi
    done

    # ensure output exists even if no scaffolds passed coverage filter
    touch ${strain}.perchrom.nquire.tsv
    """

    stub:
    """
    printf 'strain\\tscaffold\\tdip\\trip\\ttet\\tbest\\n' > ${strain}.perchrom.nquire.tsv
    """
}

/*
 * Per-gene nQuire: create across gene BED → denoise each bin →
 * lrdmodel each bin → concatenate into one TSV.
 * Skipped automatically when no gene BED is found for a strain.
 */
process NQUIRE_PER_GENE {
    tag { strain }
    label 'nquire'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(strain), path(bam), path(bai), path(gene_bed)

    output:
    tuple val(strain), path("${strain}.pergene.nquire.tsv")

    script:
    """
    mkdir -p pergene pergene_denoised
    module load nQuire
    nQuire create -b ${bam} -r ${gene_bed} \
        -o pergene -c ${params.min_cov_gene} -q ${params.min_mapq}

    header_written=0
    for bin in pergene/*.bin; do
        [ -f "\$bin" ] || continue
        pref=\$(basename "\$bin" .bin)
        nQuire denoise -o pergene_denoised "\$bin"
        denoised="pergene_denoised/\${pref}.denoised.bin"
        [ -f "\$denoised" ] || continue

        if [ \$header_written -eq 0 ]; then
            nQuire lrdmodel -t ${task.cpus} "\$denoised" \
                | awk -v s="${strain}" -v g="\$pref" '
                    NR==1 { sub(/^[^\\t]+/, "strain\\tgene"); print; next }
                          { sub(/^[^\\t]+/, s "\\t" g);       print }
                ' >> ${strain}.pergene.nquire.tsv
            header_written=1
        else
            nQuire lrdmodel -t ${task.cpus} "\$denoised" \
                | awk -v s="${strain}" -v g="\$pref" '
                    NR==1 { next }
                          { sub(/^[^\\t]+/, s "\\t" g); print }
                ' >> ${strain}.pergene.nquire.tsv
        fi
    done

    touch ${strain}.pergene.nquire.tsv
    """

    stub:
    """
    printf 'strain\\tgene\\tdip\\trip\\ttet\\tbest\\n' > ${strain}.pergene.nquire.tsv
    """
}

// ── Workflow ─────────────────────────────────────────────────────────────────

workflow {

    // Discover BAMs; derive strain name by stripping .bam suffix.
    // .take(-1) means "all"; params.n_test > 0 limits to the first N for quick tests.
    bam_ch = Channel
        .fromPath("${params.alndir}/*.bam", checkIfExists: true)
        .map { bam ->
            def strain = bam.baseName
            def bai    = file("${bam}.bai", checkIfExists: false)
            if ( !bai.exists() ) bai = file("${bam.parent}/${strain}.bai", checkIfExists: false)
            tuple(strain, bam, bai)
        }
        .take( (params.n_test as int) > 0 ? params.n_test as int : -1 )

    // 1. Whole-genome analysis (no scaffold BED needed)
    NQUIRE_WHOLE_GENOME(bam_ch)

    // 2. Per-scaffold: derive scaffold BED from each BAM header
    BAM_TO_SCAFFOLD_BED(bam_ch)
    NQUIRE_PER_SCAFFOLD(BAM_TO_SCAFFOLD_BED.out)

    // 3. Per-gene: join BAMs with the matching gene BED (skip strains with no BED)
    gene_bed_ch = Channel
        .fromPath("${params.beddir}/*.bed", checkIfExists: false)
        .map { bed -> tuple(bed.baseName, bed) }

    bam_with_gene_ch = bam_ch
        .join(gene_bed_ch, by: 0, remainder: false)   // inner join: strains with a BED only

    NQUIRE_PER_GENE(bam_with_gene_ch)
}

Ploidy analysis of Basidiobolus

Separate analysis for Lluvia Vargas paper on Basidiobolus genomes

Jason Stajich, jason.stajich<at>ucr.edu

---

## nQuire Ploidy Estimation (Nextflow)

### What was built

| File | Purpose |
|---|---|
| `pipeline/ploidy/nquire.nf` | DSL2 workflow with 4 processes |
| `nextflow.config` | Root config: SLURM executor, shared defaults, profiles |
| `conf/profile_nquire.config` | Per-run params, resources, trace/report/timeline paths |
| `conf/test_nquire.config` | Stub/test profile (local executor, no modules) |
| `run_nquire.sh` | sbatch launcher for the Nextflow head process |

### Workflow processes

| Process | Queue | CPUs/Mem | Output |
|---|---|---|---|
| `BAM_TO_SCAFFOLD_BED` | `short` | 2 / 4 GB | per-strain scaffold BED from BAM header |
| `NQUIRE_WHOLE_GENOME` | `epyc` | 8 / 24 GB | `results/nquire/<strain>.nquire.tsv` |
| `NQUIRE_PER_SCAFFOLD` | `epyc` | 8 / 24 GB | `results/nquire/<strain>.perchrom.nquire.tsv` |
| `NQUIRE_PER_GENE` | `epyc` | 8 / 24 GB | `results/nquire/<strain>.pergene.nquire.tsv` |

Inputs are discovered automatically from `aln/*.bam`. Gene BED files are matched by strain name from `bed/<strain>.bed`; strains without a matching BED file are skipped for the per-gene analysis.

### How to run

Validate the workflow graph without submitting any jobs (runs in seconds on the login node):

```bash
nextflow run pipeline/ploidy/nquire.nf -c nextflow.config -profile test -stub-run
```

Full run (submits SLURM jobs via the Nextflow head process):

```bash
sbatch run_nquire.sh
```

Pass extra parameters inline, e.g. to limit to 3 samples for a quick test:

```bash
sbatch run_nquire.sh --n_test 3
```

Resume a partially completed run (Nextflow caches completed tasks):

```bash
sbatch run_nquire.sh   # -resume is already included in the launcher
```

### Output

All results land in `results/nquire/`. Nextflow logs (trace, HTML report, timeline) go to `logs/nextflow/`.

---

## Per-gene ploidy plots

`pipeline/plot/plot_pergene_ploidy.R` turns the per-gene nQuire tables into
contig-organized ploidy figures (one PDF per strain plus a combined overview).

### Input expected

| Path | Content |
|---|---|
| `results/nquire2/<strain>.pergene.nquire.tsv.gz` | `gene, free, dip, tri, tet, d_dip, d_tri, d_tet` |
| `bed/<strain>.bed` | `contig, start, end, gene` (gene names match the per-gene table) |
| `genome/<strain>.masked.fasta.fai` | optional; supplies true contig lengths for "by size" ordering. If absent, contig length is taken from the largest gene end coordinate in the BED. |

### Winner / missing-data rules

For each gene the **winning ploidy is the smallest delta-log-likelihood** among
`d_dip` (2n), `d_tri` (3n), `d_tet` (4n). Values of `0`, `nan`, and `-nan` are
treated as missing/unestimable; a gene whose three deltas are all missing gets
no call and is omitted from the plots (it shows up as a gap).

### How to run

```bash
module load R          # R 4.5.x with ggplot2, data.table, scales

# all strains, default = winning ploidy on the y-axis
Rscript pipeline/plot/plot_pergene_ploidy.R

# same plots but with the delta-log-likelihood score on the y-axis
Rscript pipeline/plot/plot_pergene_ploidy.R --mode score

# limit to one strain and/or change the number of faceted contigs
Rscript pipeline/plot/plot_pergene_ploidy.R --strain CBS931.73_pb --top 30
```

Options: `--mode ploidy|score` (default `ploidy`), `--top N` (faceted contigs,
default 50), `--strain NAME` (repeatable; default = all strains found).

### Output (`results/nquire2/plots/`)

| File | Contents |
|---|---|
| `<strain>.pergene_ploidy.<mode>.pdf` | **p1** top-`N` contigs by size, faceted, genes by position (coloured by ploidy call). **p2** genome-wide Manhattan, contigs ordered by size, **alternating black/red per contig** so same-contig genes group together. |
| `ALL_strains.pergene_ploidy.<mode>.pdf` | **p1** stacked composition bar of ploidy calls per strain. **p2** genome-wide Manhattan faceted by strain for cross-strain comparison. |
| `ploidy_call_summary.tsv` | per-strain counts: total / called / missing / n2 / n3 / n4. |

### Notes / findings

- **These assemblies are fragmented — there are no true chromosomes.** Each
  genome has ~1,000–5,900 contigs/scaffolds with at most ~40–160 genes on the
  largest one, so "chromosome-organized" is really "contig-organized". The
  faceted top-`N` page is the clean primary view; in the Manhattan the per-contig
  black/red alternation is only visually resolvable on the larger (left-hand)
  contigs.
- **Empty inputs are skipped automatically.** The long-read per-gene tables
  (`*_ont`, `NRRL2992_pb`, `STP1710.7`) are 0-row and are skipped. `NRRL2992`
  has no `.fai`, so its contig sizes come from the BED, and it has very few
  estimable genes (1,540 / 19,174).
- **Ploidy call summary** (estimable genes; winner = smallest delta-LL):

  | strain | called | missing | 2n | 3n | 4n | dominant |
  |---|---|---|---|---|---|---|
  | CBS931.73_pb | 12,270 | 4,256 | 230 | 1,000 | 11,040 | strongly 4n (90%) |
  | Bran_AGB5 | 11,991 | 6,644 | 2,482 | 1,704 | 7,805 | 4n |
  | CBS931.73 | 9,161 | 7,365 | 2,268 | 2,054 | 4,839 | 4n / mixed |
  | STP1717.1 | 3,384 | 13,264 | 985 | 615 | 1,784 | 4n / mixed |
  | NRRL2992 | 1,540 | 17,634 | 516 | 279 | 745 | mixed (sparse) |
  | UHM260.5136 | 7,055 | 10,064 | 3,167 | 1,149 | 2,739 | 2n / mixed |
  | UHM207.4505 | 6,176 | 13,247 | 3,113 | 783 | 2,280 | 2n |
  | UHM516.7697 | 6,026 | 13,213 | 3,181 | 717 | 2,128 | 2n |
  | UHM520.7734 | 7,969 | 10,635 | 4,752 | 785 | 2,432 | 2n |

  The UHM strains skew diploid; the CBS/Bran strains skew tetraploid. The same
  CBS931.73 isolate looks more strongly tetraploid in its PacBio assembly
  (`_pb`) than in the short-read assembly, which is worth keeping in mind when
  comparing calls across assembly types.

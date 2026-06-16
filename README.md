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
- **Ploidy call summary** (estimable genes; winner = smallest delta-LL;
  percentages are of *called* genes). Source: `results/nquire2/plots/ploidy_call_summary.tsv`.

  | strain | total | called | missing | 2n | 3n | 4n | dominant |
  |---|---|---|---|---|---|---|---|
  | CBS931.73_pb | 16,526 | 12,270 | 4,256 | 232 (2%) | 997 (8%) | 11,041 (90%) | strongly 4n |
  | Bran_AGB5 | 18,635 | 11,988 | 6,647 | 2,483 (21%) | 1,699 (14%) | 7,806 (65%) | 4n |
  | CBS931.73 | 16,526 | 9,161 | 7,365 | 2,265 (25%) | 2,069 (23%) | 4,827 (53%) | 4n / mixed |
  | STP1717.1 | 16,648 | 3,385 | 13,263 | 986 (29%) | 616 (18%) | 1,783 (53%) | 4n / mixed |
  | NRRL2992 | 19,174 | 1,540 | 17,634 | 516 (34%) | 279 (18%) | 745 (48%) | mixed (sparse) |
  | UHM260.5136 | 17,119 | 7,056 | 10,063 | 3,160 (45%) | 1,155 (16%) | 2,741 (39%) | 2n / mixed |
  | UHM207.4505 | 19,423 | 6,176 | 13,247 | 3,113 (50%) | 782 (13%) | 2,281 (37%) | 2n |
  | UHM516.7697 | 19,239 | 6,026 | 13,213 | 3,183 (53%) | 716 (12%) | 2,127 (35%) | 2n |
  | UHM520.7734 | 18,604 | 7,970 | 10,634 | 4,750 (60%) | 789 (10%) | 2,431 (31%) | 2n |
  | STP1710.7 | 27,000 | 7,345 | 19,655 | 4,618 (63%) | 327 (4%) | 2,400 (33%) | 2n |

  The UHM strains and STP1710.7 skew diploid; the CBS/Bran strains skew
  tetraploid. The same CBS931.73 isolate looks more strongly tetraploid in its
  PacBio assembly (`_pb`, 90% 4n) than in the short-read assembly (53% 4n),
  which is worth keeping in mind when comparing calls across assembly types.
  Note also that "called" is often a minority of "total" genes — many genes are
  unestimable (all three deltas missing) and are excluded from the percentages,
  so strains with low `called` counts (e.g. NRRL2992, STP1717.1) carry the most
  uncertainty.

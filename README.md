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

# Interpreting nQuire ploidy calls in *Basidiobolus*: findings and recommendations

Status as of 2026-10-06. Based on `results/nquire2/` (whole-genome, per-scaffold
and per-gene nQuire runs) and `results/mosdepth/` (per-gene read depth).

## Question

Many genes, and some whole genomes, give nQuire calls whose best-fitting
model is tetraploid (allele balance near 25:75). Is a tetraploid genome the
only way to get 25:75 allele frequencies?

**Short answer: no.** Several other mechanisms give 25:75, and for most
strains the per-gene 4n calls are better explained by read-depth artefacts
than by genome-wide tetraploidy. CBS931.73 is the one strain where a
genome-wide 4n-like signal holds up, but a mixed culture would look the same
to nQuire.

## Ways to get 25:75 besides tetraploidy

| Cause | How it makes 25:75 | How to recognise it |
| --- | --- | --- |
| Collapsed paralogs or repeats | Several diverged copies map to one reference locus; one variant copy out of four reads gives 25% | Read depth ~2× or more |
| Assembly keeps both haplotypes as separate contigs | Reads split between the two copies, plus some cross-mapping | About half depth, often across whole scaffolds |
| Mixed or non-clonal culture | A 1:3 mix of two haploid genotypes puts every difference between them at 25% | Same frequency genome-wide; peaks at 0.25 and 0.75 with no 0.5 peak; changes after single-spore re-isolation |
| Partial aneuploidy or segmental duplication | Trisomy gives 33:67; extra copies of one homolog give 25:75 | Long runs of same-call genes along scaffolds, with matching depth steps |
| Low depth and few SNPs | Random sampling noise in small per-gene samples | Mostly in low-depth genes |
| How nQuire's models work | The 4n model has peaks at 0.25, 0.5 and 0.75, so any broad or noisy spread fits it better than the 2n model | 4n wins wherever the data are noisy (e.g. long-read errors) |
| True autotetraploid (AAAB, AABB) | — | Normal depth, genome-wide, peaks at 0.25, 0.5 and 0.75 together, no spatial clustering |

## Findings so far

### 1. Short- vs long-read concordance (`scripts/nquire_platform_concordance.py`)

Outputs: `results/nquire2/concordance/`

- Long-read runs (ONT, PacBio reads mapped to the same reference) call many
  more genes 4n than the short-read run of the same strain (e.g. CBS931.73
  53% → 90%; UHM207.4505 37% → 66%). Every ONT whole-genome run is called 4n,
  including strains whose short-read whole-genome call is 2n. Read error
  spreading allele frequencies is the likely cause, so allele-frequency ploidy
  inference should use short reads (or HiFi).
- Per-gene calls agree poorly between platforms (Cohen's κ 0.06–0.21; 0.03–0.29
  for confident calls only). Individual per-gene calls should not be read as
  gene-level ploidy.

### 2. Read depth and spatial clustering of per-gene calls (`scripts/nquire_depth_spatial.py`)

Outputs: `results/nquire2/depth_spatial/`

- **Per-gene 4n calls track depth.** Fraction of confident (margin ≥ 5)
  per-gene calls that are 4n, by gene depth relative to the sample median
  (short reads):

  | Sample | <0.5× | 0.5–0.75× | 0.75–1.25× | 1.25–1.5× | 1.5–2× | >2× |
  | --- | --- | --- | --- | --- | --- | --- |
  | Bran_AGB5 | 1.00 | 0.99 | 0.71 | 0.42 | 0.36 | 0.59 |
  | CBS931.73 | 0.49 | 0.71 | 0.57 | 0.71 | 0.64 | 0.47 |
  | STP1710.7 | 0.90 | 0.73 | 0.16 | 0.21 | 0.12 | 0.79 |
  | STP1717.1 | 0.78 | 0.47 | 0.42 | 0.52 | 0.56 | 0.52 |
  | UHM207.4505 | 0.94 | 0.88 | 0.37 | 0.15 | 0.22 | 0.55 |
  | UHM260.5136 | 0.92 | 0.79 | 0.29 | 0.27 | 0.45 | 0.55 |
  | UHM516.7697 | 0.88 | 0.66 | 0.18 | 0.24 | 0.31 | 0.77 |
  | UHM520.7734 | 0.69 | 0.53 | 0.18 | 0.55 | 0.48 | 0.82 |

  Low-depth and high-depth genes are mostly called 4n; genes at single-copy
  depth mostly are not, except in CBS931.73 (and partly STP1717.1), where the
  4n fraction is flat across depth.

- **A group of half-depth genes.** Each strain has about 2,000–3,300 genes at
  ~0.6× depth in *both* short- and long-read data, called 4n on both platforms
  (62–80% of short-read calls). Many sit on whole scaffolds at half depth
  (e.g. ~300 of 682 scaffolds with ≥10 genes in UHM516.7697; ~60% of
  half-depth genes in CBS931.73, STP1710.7 and UHM516.7697 are on such
  scaffolds). The leading hypothesis is that the assemblies keep both
  haplotypes (or diverged paralogs) as separate contigs, with reads split
  between them and some cross-mapping. **Not yet confirmed.** If true, part
  of the high duplicated-gene content in these assemblies may also be
  haplotype duplication rather than real paralogs.
- **Weak spatial clustering.** Same-call adjacency for short-read calls is
  only slightly above permutation (e.g. 0.63 vs 0.61 in UHM520.7734); the
  longest run of 4n genes is ≤10 genes except in CBS931.73 (20). No evidence
  for large aneuploid or duplicated segments.

### 3. Per-strain reading (short reads)

| Strain | Whole-genome call (relfit) | Reading |
| --- | --- | --- |
| CBS931.73 | 4n (0.22, good fit) | Genome-wide 4n-like signal at every depth. Candidate tetraploid **or** mixed culture; needs the 0.5-peak test and k-mer analysis |
| UHM520.7734 | 2n (0.20, good fit) | Diploid; per-gene 4n calls are depth artefacts |
| UHM516.7697 | 2n (0.60) | Probably diploid |
| UHM207.4505 | 2n (0.86, poor fit) | Probably diploid; whole-genome fit poor |
| UHM260.5136 | 4n (0.89, poor fit, small margin) | Unresolved; per-gene pattern looks diploid |
| STP1710.7 | no estimate | Per-gene pattern looks diploid; rerun whole-genome nQuire |
| STP1717.1 | no estimate | Unresolved (42% 4n even at single-copy depth); rerun whole-genome nQuire |
| Bran_AGB5 | 4n (0.92, poor fit) | Unusual: many scattered half-depth genes, almost all called 4n. Check that the reference matches this strain, and check culture purity |
| NRRL2992 | 2n (0.21) | Few genes called (1,540). `samples.csv` lists only PacBio reads for NRRL2992, so confirm what the un-suffixed BAM contains |

"relfit" = best fixed-model ΔLL / free-model LL; near 0 means a fixed ploidy
fits about as well as the free mixture, near 1 means none fits.

## Recommendations (most useful first)

1. **Look at the allele-frequency spectra directly.** Run `nQuire histo` and
   `nQuire estmodel` on the existing whole-genome denoised `.bin` files, split
   by depth class (half, single-copy, ≥2×). A 0.5 peak alongside 0.25/0.75
   points to tetraploidy; no 0.5 peak points to a 1:3 mixture. This decides
   CBS931.73.
2. **Reference-free k-mer analysis (GenomeScope2 + Smudgeplot)** on the short
   reads: `pipeline/ploidy/kmer_ploidy.sh` (SLURM array, 8 short-read
   samples). Compare GenomeScope2 fits at p = 2, 3, 4 and see where the main
   smudge falls (AB vs AAB vs AAAB/AABB). Collapsed-paralog pairs appear as
   smudges at 2× coverage rather than as AAAB at 1×.
3. **Check the assemblies for retained haplotigs**: Merqury spectra-cn plot,
   purge_dups coverage histogram, BUSCO duplicated fraction. Confirms or rules
   out the half-depth explanation.
4. **Re-run nQuire on a cleaner gene set**: short reads only, single-copy
   genes at 0.75–1.25× depth (or BUSCO single-copy), higher minimum depth
   (`-c 20`), pooled per scaffold rather than per gene so each fit has enough
   sites. Also rerun the failed whole-genome runs for STP1710.7 and STP1717.1.
5. **Read-backed phasing with long reads** (e.g. WhatsHap polyphase) in
   regions with dense SNPs: a true tetraploid shows up to four haplotypes;
   collapsed paralogs give pairs.
6. **Lab checks**: single-spore re-isolation of CBS931.73 and re-sequencing
   (tests mixed culture); nuclear DNA content by flow or image cytometry.

## Files

| Path | Contents |
| --- | --- |
| `scripts/nquire_platform_concordance.py` | short- vs long-read call concordance |
| `scripts/nquire_depth_spatial.py` | per-gene call vs depth, spatial clustering |
| `pipeline/ploidy/kmer_ploidy.sh` | FastK → GenomeScope2 (p = 2,3,4) + Smudgeplot, SLURM array |
| `results/nquire2/concordance/` | concordance outputs |
| `results/nquire2/depth_spatial/` | depth/spatial outputs, incl. `depth_by_call.pdf` |

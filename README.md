# zAMPExplorer

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-universe Status](https://metagenlab.r-universe.dev/badges/zAMPExplorer)](https://metagenlab.r-universe.dev)

**zAMPExplorer** is a Shiny application for reproducible downstream exploration and statistical analysis of 16S rRNA amplicon microbiome data stored as a [`phyloseq`](https://bioconductor.org/packages/phyloseq/) object. It is designed to work directly with the phyloseq output produced by the [zAMP pipeline](https://zamp.readthedocs.io/en/latest/), but any compatible phyloseq object can be uploaded.

## What changed in 0.2.0

Version 0.2.0 reorganizes the application into independent Shiny modules and introduces one publication plotting/export layer shared by the interface and downloaded figures. This removes duplicated plotting code and makes screen and saved figures consistent.

Main improvements:

- modular application architecture instead of one monolithic server file;
- publication-oriented figure theme and colour palettes;
- vector PDF/SVG and 600-dpi PNG export from the same underlying plot object used in the app;
- optional automatic saving of figures and result tables to a results directory;
- restored heatmap and MaAsLin2 differential-abundance modules;
- read-depth filtering propagated consistently to all downstream modules;
- updated beta-diversity workflow with PERMANOVA and dispersion testing;
- publication-style DMM driver plots and RDA/dbRDA plots;
- reproducible Docker build that installs the checked-out source rather than a previously published release;
- automated R package checks on pull requests.

## Analyses

zAMPExplorer currently provides:

1. **Dataset overview** — metadata, abundance/taxonomy table and basic dataset summaries.
2. **Read QC** — sequencing-depth distributions, group comparisons, read-depth filtering and rarefaction curves.
3. **Taxa overview** — prevalence, mean abundance, core taxa and shared/unique taxa summaries.
4. **Composition** — filtered taxonomic relative-abundance barplots.
5. **Heatmap** — publication heatmaps using relative abundance, log abundance or CLR values.
6. **Alpha diversity** — Observed, Chao1, ACE, Shannon, Simpson, inverse Simpson and Fisher diversity with pairwise Wilcoxon tests and BH correction.
7. **Beta diversity** — PCoA, PCA and NMDS with Bray-Curtis, Jaccard or Aitchison distance where appropriate; PERMANOVA and beta-dispersion tests.
8. **Differential abundance** — multivariable association testing with MaAsLin2.
9. **Community typing** — Dirichlet multinomial mixture (DMM) model selection, assignments and community-type drivers.
10. **RDA / dbRDA** — constrained ordination with permutation tests and leading taxon arrows.

## Input

Upload an `.rds` file containing a `phyloseq` object with, at minimum:

- `otu_table()` / ASV abundance table;
- `tax_table()`;
- `sample_data()`.

A phylogenetic tree and reference sequences are supported but are not required for the current analysis modules.

The app removes zero-depth samples and zero-abundance taxa after upload. Further read-depth filtering can be applied in **Read QC** and is then propagated to all downstream tabs.

## Installation

### Recommended: install the released package

```r
install.packages(
  "zAMPExplorer",
  repos = c(
    "https://metagenlab.r-universe.dev",
    "https://cloud.r-project.org"
  )
)

library(zAMPExplorer)
zAMPExplorer_app()
```

### Development version from GitHub

```bash
git clone https://github.com/metagenlab/zAMPExplorer.git
cd zAMPExplorer
```

Then in R:

```r
install.packages("pak")
pak::local_install(".")

library(zAMPExplorer)
zAMPExplorer_app()
```

For active development:

```r
install.packages("devtools")
devtools::load_all()
zAMPExplorer_app()
```

If needed, install dependencies with:

```r
source("install_dependencies.R")
```

## Docker

Build from the repository root:

```bash
docker build -t zampexplorer:latest .
```

Run it:

```bash
docker run --rm \
  -p 3838:3838 \
  -v "$(pwd)/results:/results" \
  zampexplorer:latest
```

Open `http://localhost:3838` and set the results folder to `/results` if automatically generated outputs should persist on the host.

## Publication-quality exports

Interactive views are generated from publication plot objects wherever possible. Download controls allow:

- **PDF** — preferred vector format for manuscripts;
- **SVG** — editable vector graphics;
- **PNG** — 600 dpi raster output.

Width and height are specified in inches. When automatic saving is enabled on the upload tab, figures are additionally saved as PDFs in the selected results folder. Heatmaps use `ComplexHeatmap` and have their own PDF/SVG/600-dpi PNG exporter.

## MaAsLin2 notes

The differential-abundance tab exposes fixed and random effects, normalization, transformation, analysis method, abundance/prevalence filtering, reference levels and the q-value threshold. For categorical variables, reference levels use the MaAsLin2 format, for example:

```text
group,Control;sex,Female
```

Choose covariates, transformations and reference levels according to the experimental design rather than treating app defaults as a universal prescription.

## Development structure

```text
zAMPExplorer/
├── R/
│   ├── app_config.R
│   ├── app_ui.R
│   ├── app_server.R
│   ├── run_app.R
│   ├── phyloseq_utils.R
│   ├── plot_utils.R
│   ├── mod_overview.R
│   ├── mod_qc.R
│   ├── mod_taxa.R
│   ├── mod_composition.R
│   ├── mod_heatmap.R
│   ├── mod_alpha.R
│   ├── mod_beta.R
│   ├── mod_maaslin2.R
│   ├── mod_dmm.R
│   └── mod_rda.R
├── inst/app/www/zamp.css
├── tests/testthat/
├── .github/workflows/R-CMD-check.yaml
├── DESCRIPTION
├── NAMESPACE
├── Dockerfile
└── env.yml
```

Each module creates its analysis result once and reuses the same reactive plot for interactive display, download and automatic saving.

## Testing

From the repository root:

```r
install.packages("devtools")
devtools::test()
devtools::check()
```

Pull requests also run an automated R CMD check on GitHub Actions.

## Citation and acknowledgements

zAMPExplorer integrates established R/Bioconductor packages including `phyloseq`, `vegan`, `Maaslin2`, `DirichletMultinomial`, `ComplexHeatmap`, `ggplot2` and `plotly`. Please cite the relevant methods/packages when analyses produced by zAMPExplorer are used in publications.

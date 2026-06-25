# Breast cancer proteomics data (TCGA, RPPA)

Reverse-phase protein array (RPPA) measurements of 163 proteins across
346 breast cancer tumor samples from The Cancer Genome Atlas (TCGA),
along with the PAM50 molecular subtype of each sample. Used as the
running example in Tous & Chiquet (2026).

## Usage

``` r
brca_rppa
```

## Format

A list with 3 elements:

- expr:

  a 346 x 163 numeric matrix of protein expression levels (normalized
  log-ratios), samples in rows (named by TCGA sample identifier),
  proteins in columns.

- covariates:

  a data frame with one row per sample, in the same order as the rows of
  \`expr\`: \`sampleId\` (matches \`rownames(expr)\`) and 11 clinical
  variables, factors unless noted otherwise – \`CN_CLUSTER\`
  (copy-number cluster, 5 levels), \`CONVERTED_STAGE\` (tumor stage),
  \`ER_STATUS\` (estrogen receptor status), \`FRACTION_GENOME_ALTERED\`
  (numeric, in \\0, 1\\), \`HER2_STATUS\`, \`METASTASIS_CODED\`,
  \`METHYLATION_CLUSTER\` (5 levels), \`MUTATION_COUNT\` (numeric),
  \`PAM50_SUBTYPE\` (5 levels: Basal-like, HER2-enriched, Luminal A,
  Luminal B, Normal-like), \`RPPA_CLUSTER\` (7 levels) and \`AGE\`
  (numeric, in years). A few of these have missing values for a handful
  of samples (\`FRACTION_GENOME_ALTERED\`, \`HER2_STATUS\`,
  \`METASTASIS_CODED\`).

- gene_annotation:

  a data frame with one row per protein, in the same order as the
  columns of \`expr\`: \`protein\` (matches \`colnames(expr)\`),
  \`entrezGeneId\`, \`hugoGeneSymbol\` and \`type\` (\`"protein"\` or
  \`"phosphoprotein"\`).

## Source

TCGA breast cancer cohort (Cancer Genome Atlas Network, 2012,
[doi:10.1038/nature11412](https://doi.org/10.1038/nature11412) ),
downloaded via cBioPortal (\`brca_tcga_pub\` study,
<https://www.cbioportal.org>).

## Examples

``` r
expr <- brca_rppa$expr
X <- model.matrix(~ 0 + PAM50_SUBTYPE, data = brca_rppa$covariates)
nb_data <- NormalBlockData$new(expr, X)
```

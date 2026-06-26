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
  \`entrezGeneId\`, \`hugoGeneSymbol\`, \`type\` (\`"protein-coding"\`
  or \`"phosphoprotein"\`) and \`go_bp_term\` (a Gene Ontology
  Biological Process classification, one term per protein, or
  \`"unknown"\` if none could be found – see Details).

## Source

TCGA breast cancer cohort (Cancer Genome Atlas Network, 2012,
[doi:10.1038/nature11412](https://doi.org/10.1038/nature11412) ),
downloaded via cBioPortal (\`brca_tcga_pub\` study,
<https://www.cbioportal.org>). \`go_bp_term\` (see Details) was added
afterwards from \`org.Hs.eg.db\`.

## Details

\`go_bp_term\` is derived from \`org.Hs.eg.db\`, queried by
\`entrezGeneId\` for \`"protein-coding"\` rows. Phospho-site antibodies
(\`type == "phosphoprotein"\`, e.g. \`"EGFR_PY1068"\`) carry a
placeholder, negative \`entrezGeneId\` (there is no separate gene for a
specific phosphorylation site, only the underlying gene) and are instead
queried by gene symbol (the part of \`protein\` before the first
\`"\_"\`). Each protein typically has dozens of GO Biological Process
terms; \`go_bp_term\` keeps, for each protein, the one term shared by
the most \*other\* proteins in this dataset – giving a handful of
non-trivial groups rather than one near-singleton group per protein,
more useful for comparing against an inferred clustering. See
\`data-raw/brca_rppa_go_annotation.R\` for the full extraction code.

## Examples

``` r
expr <- brca_rppa$expr
X <- model.matrix(~ 0 + PAM50_SUBTYPE, data = brca_rppa$covariates)
nb_data <- NormalBlockData$new(expr, X)
table(brca_rppa$gene_annotation$go_bp_term)
#> 
#>                                                DNA damage response 
#>                                                                 12 
#>                              G1/S transition of mitotic cell cycle 
#>                                                                  5 
#>                                                  apoptotic process 
#>                                                                 19 
#>                                                      axon guidance 
#>                                                                  1 
#>                                                      cell adhesion 
#>                                                                  3 
#>                                               cell differentiation 
#>                                                                  2 
#>                                                      cell division 
#>                                                                  1 
#>                                    cell-cell junction organization 
#>                                                                  1 
#>                              cellular response to insulin stimulus 
#>                                                                  1 
#>                                  intracellular signal transduction 
#>                                                                  1 
#>                                            lipid metabolic process 
#>                                                                  2 
#>                           negative regulation of apoptotic process 
#>                                                                 21 
#>               negative regulation of cell population proliferation 
#>                                                                  1 
#>                             negative regulation of gene expression 
#>                                                                  2 
#>          negative regulation of transcription by RNA polymerase II 
#>                                                                  1 
#> phosphatidylinositol 3-kinase/protein kinase B signal transduction 
#>                                                                  1 
#>                 positive regulation of DNA-templated transcription 
#>                                                                  1 
#>               positive regulation of cell population proliferation 
#>                                                                  5 
#>                             positive regulation of gene expression 
#>                                                                  1 
#>          positive regulation of transcription by RNA polymerase II 
#>                                                                 15 
#>                                            protein phosphorylation 
#>                                                                  1 
#>                          regulation of DNA-templated transcription 
#>                                                                  2 
#>                                    regulation of apoptotic process 
#>                                                                  2 
#>                   regulation of transcription by RNA polymerase II 
#>                                                                  1 
#>                                          regulation of translation 
#>                                                                  1 
#>                                       response to oxidative stress 
#>                                                                  1 
#>                                    response to xenobiotic stimulus 
#>                                                                  3 
#>                                                signal transduction 
#>                                                                 56 
```

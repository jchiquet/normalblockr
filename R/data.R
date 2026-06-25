#' Breast cancer proteomics data (TCGA, RPPA)
#'
#' Reverse-phase protein array (RPPA) measurements of 163 proteins across 346
#' breast cancer tumor samples from The Cancer Genome Atlas (TCGA), along with
#' the PAM50 molecular subtype of each sample. Used as the running example in
#' Tous & Chiquet (2026) and in
#' `inst/NORMAL_BLOCK_PAPER_ILLUSTRATIONS/1_breast_cancer_proteomics`.
#'
#' @format A list with 3 elements:
#' \describe{
#'   \item{expr}{a 346 x 163 numeric matrix of protein expression levels
#'   (normalized log-ratios), samples in rows (named by TCGA sample
#'   identifier), proteins in columns.}
#'   \item{covariates}{a data frame with one row per sample, in the same
#'   order as the rows of `expr`: `sampleId` (matches `rownames(expr)`) and 11
#'   clinical variables, factors unless noted otherwise --
#'   `CN_CLUSTER` (copy-number cluster, 5 levels), `CONVERTED_STAGE` (tumor
#'   stage), `ER_STATUS` (estrogen receptor status), `FRACTION_GENOME_ALTERED`
#'   (numeric, in \[0, 1\]), `HER2_STATUS`, `METASTASIS_CODED`,
#'   `METHYLATION_CLUSTER` (5 levels), `MUTATION_COUNT` (numeric),
#'   `PAM50_SUBTYPE` (5 levels: Basal-like, HER2-enriched, Luminal A, Luminal
#'   B, Normal-like), `RPPA_CLUSTER` (7 levels) and `AGE` (numeric, in
#'   years). A few of these have missing values for a handful of samples
#'   (`FRACTION_GENOME_ALTERED`, `HER2_STATUS`, `METASTASIS_CODED`).}
#'   \item{gene_annotation}{a data frame with one row per protein, in the
#'   same order as the columns of `expr`: `protein` (matches `colnames(expr)`),
#'   `entrezGeneId`, `hugoGeneSymbol` and `type` (`"protein"` or
#'   `"phosphoprotein"`).}
#' }
#' @source TCGA breast cancer cohort (Cancer Genome Atlas Network, 2012,
#' \doi{10.1038/nature11412}), downloaded via cBioPortal (`brca_tcga_pub`
#' study, \url{https://www.cbioportal.org}).
#' @examples
#' expr <- brca_rppa$expr
#' X <- model.matrix(~ 0 + PAM50_SUBTYPE, data = brca_rppa$covariates)
#' nb_data <- NormalBlockData$new(scale(expr), X)
"brca_rppa"

#' French stream fish community data (ONEMA / OFB electrofishing surveys)
#'
#' Fish biomass per species, aggregated by sampling station, together with
#' environmental covariates for each station. Derived from the fish community
#' monitoring of stream sections across metropolitan France (1995-2018) by
#' the French Office for Biodiversity (formerly ONEMA, "Office National de
#' l'Eau et des Milieux Aquatiques"), using electrofishing.
#'
#' The source data is organized around fishing *operations* (`opcod`):
#' several electrofishing operations can be carried out at the same
#' *station* over time. `onema` aggregates this to one row per station:
#' every operation is first mapped to its station (`fishing_protocol.csv`),
#' then biomass (in grams) is summed over every operation recorded at a
#' given station, separately for each species (i.e. `biomass[s, ]` is the
#' total biomass of each species ever caught at station `s`, not a single
#' operation's catch). Stations without environmental data are dropped.
#' Species observed at fewer than 10 stations are also dropped: at that
#' level of rarity, a zero-inflated fit's iterative initialization can fit a
#' species' handful of non-zero observations exactly, driving its estimated
#' noise precision to infinity (these species also carry essentially no
#' information for clustering anyway).
#'
#' @format A list with 2 elements:
#' \describe{
#'   \item{biomass}{a 399 x 46 numeric matrix of total biomass (grams) per
#'   species (3-letter species code, columns) and station (row names).}
#'   \item{covariates}{a data frame with one row per station, in the same
#'   order as the rows of `biomass`: `station` (matches `rownames(biomass)`)
#'   and 14 environmental variables -- `slope`, `alt` (altitude), `d_source`
#'   (distance to source), `strahler` (Strahler stream order),
#'   `width_river_mean`/`width_river_cv`, `avg_depth_station_mean`/
#'   `avg_depth_station_cv`, `DBO_med`/`DBO_cv` (biological oxygen demand),
#'   `flow_med`/`flow_cv` and `temperature_med`/`temperature_cv` (`_med` =
#'   median, `_cv` = coefficient of variation, over the monitoring period).}
#' }
#' @source Danet, A., Mouchet, M., Bonnaffe, W., Thebault, E., Fontaine, C.
#' (2021) "Species richness and food-web structure jointly drive community
#' biomass and its temporal stability in fish communities", data set,
#' Zenodo, \doi{10.5281/zenodo.5095656}.
#' @examples
#' Y <- log(1 + onema$biomass)
#' X <- model.matrix(~ 1, data = onema$covariates)
#' nb_data <- NormalBlockData$new(Y, X)
"onema"

#' Breast cancer proteomics data (TCGA, RPPA)
#'
#' Reverse-phase protein array (RPPA) measurements of 163 proteins across 346
#' breast cancer tumor samples from The Cancer Genome Atlas (TCGA), along with
#' the PAM50 molecular subtype of each sample. Used as the running example in
#' Tous & Chiquet (2026).
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
#'   `entrezGeneId`, `hugoGeneSymbol`, `type` (`"protein-coding"` or
#'   `"phosphoprotein"`) and `go_bp_term` (a Gene Ontology Biological Process
#'   classification, one term per protein, or `"unknown"` if none could be
#'   found -- see Details).}
#' }
#' @details `go_bp_term` is derived from `org.Hs.eg.db`, queried by
#' `entrezGeneId` for `"protein-coding"` rows. Phospho-site antibodies
#' (`type == "phosphoprotein"`, e.g. `"EGFR_PY1068"`) carry a placeholder,
#' negative `entrezGeneId` (there is no separate gene for a specific
#' phosphorylation site, only the underlying gene) and are instead queried by
#' gene symbol (the part of `protein` before the first `"_"`). Each protein
#' typically has dozens of GO Biological Process terms; `go_bp_term` keeps,
#' for each protein, the one term shared by the most *other* proteins in this
#' dataset -- giving a handful of non-trivial groups rather than one
#' near-singleton group per protein, more useful for comparing against an
#' inferred clustering. See `data-raw/brca_rppa_go_annotation.R` for the full
#' extraction code.
#' @source TCGA breast cancer cohort (Cancer Genome Atlas Network, 2012,
#' \doi{10.1038/nature11412}), downloaded via cBioPortal (`brca_tcga_pub`
#' study, \url{https://www.cbioportal.org}). `go_bp_term` (see Details) was
#' added afterwards from `org.Hs.eg.db`.
#' @examples
#' expr <- brca_rppa$expr
#' X <- model.matrix(~ 0 + PAM50_SUBTYPE, data = brca_rppa$covariates)
#' nb_data <- NormalBlockData$new(expr, X)
#' table(brca_rppa$gene_annotation$go_bp_term)
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
#' out <- normal_block(nb_data, 2:15, control = NB_control(clustering_init = "ward2"))
"onema"

#' University webpages text data (CMU "4 Universities" / WebKB)
#'
#' Term-frequency data derived from the student personal webpages of the CMU
#' "World Wide Knowledge Base" (WebKB) project's "4 Universities" dataset
#' (Cornell, Texas, Washington, Wisconsin computer science departments,
#' collected January 1997). Used as one of the running examples in Tous &
#' Chiquet (2026).
#'
#' Each of the 504 pages (rows) is a document; each of the 1867 columns is a
#' term retained after lower-casing, stripping URLs/HTML tags/control
#' characters/punctuation/numbers, removing English stopwords and dropping
#' terms occurring in fewer than 2 documents. `entropies` ranks every term by
#' how evenly it is spread across documents (Shannon entropy of each term's
#' normalized document distribution, in \[0, 1\], 1 = perfectly uniform);
#' `terms` is the 100 highest-entropy terms, i.e. the terms that are the most
#' informative for distinguishing documents from one another rather than
#' just reflecting a few documents' idiosyncratic vocabulary -- the
#' transformation used by Tan et al. (2015).
#'
#' @format A list with 3 elements:
#' \describe{
#'   \item{frequencies}{a 504 x 1867 numeric matrix of term frequencies (row
#'   sums to 1): for each document (row, named by its original file path) and
#'   term (column), the fraction of that document's (post-preprocessing) word
#'   count made up of that term.}
#'   \item{entropies}{a named numeric vector of length 1867 (one value per
#'   column of `frequencies`, same names/order), the normalized Shannon
#'   entropy of each term's distribution across documents.}
#'   \item{terms}{a character vector of the 100 column names of
#'   `frequencies` with the highest `entropies`, ordered by decreasing
#'   entropy.}
#' }
#' @source CMU Text Learning Group, "World Wide Knowledge Base (Web->KB)
#' project", \url{http://www.cs.cmu.edu/~webkb/}; "4 Universities" subset,
#' \url{https://www.cs.cmu.edu/afs/cs/project/theo-20/www/data/}. Tan, P.-N.,
#' Steinbach, M., Kumar, V. (2015) "Introduction to Data Mining" (transform
#' used to derive `entropies`/`terms`).
#' @examples
#' Y <- log(1 + university$frequencies[, university$terms])
#' nb_data <- NormalBlockData$new(Y, X = matrix(1, nrow(Y), 1))
#' out <- normal_block(nb_data, 2:15)
"university"

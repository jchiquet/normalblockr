## Adds a Gene Ontology (Biological Process) classification to
## brca_rppa$gene_annotation, as a `go_bp_term` column -- one human-readable
## GO BP term per protein, or "unknown" if none could be found. Kept here as
## a trace of how that column was derived; re-running this script regenerates
## data/brca_rppa.rda from the brca_rppa object already shipped with the
## package (it does not redo the original TCGA/cBioPortal data preparation).
##
## This extraction was coded with Claude Code, originally as part of
## inst/CSDA_analyses/analysis_breast_cancer_proteomics.qmd; promoted to a
## stored column for convenience (mainly pedagogical: comparing the
## Normal-Block clustering against a known biological classification without
## redoing the org.Hs.eg.db lookup in every vignette/script that needs it).

library(normalblockr)
library(org.Hs.eg.db)

data(brca_rppa)
ga         <- brca_rppa$gene_annotation
is_phospho <- ga$type == "phosphoprotein"

# Phospho-site antibodies (e.g. "EGFR_PY1068") carry a placeholder, negative
# entrezGeneId (see gene_annotation$type): there is no separate gene for "EGFR
# phosphorylated at Y1068", only EGFR itself, so these are looked up by gene
# symbol (the part of the protein name before the first "_") instead.
ga$gene_key <- ifelse(is_phospho, sub("_.*", "", ga$protein), ga$entrezGeneId)

key_map_entrez <- dplyr::select(dplyr::filter(ga, !is_phospho), protein, gene_key)
key_map_symbol <- dplyr::select(dplyr::filter(ga, is_phospho),  protein, gene_key)

# dplyr::select()/rename()/filter() are called with an explicit namespace
# throughout this script: AnnotationDbi (loaded by org.Hs.eg.db) exports S4
# generics of the same names that otherwise silently shadow their dplyr
# counterparts once org.Hs.eg.db is attached after dplyr/tidyverse.
go_entrez <- AnnotationDbi::select(
  org.Hs.eg.db, keys = key_map_entrez$gene_key,
  columns = c("GO", "ONTOLOGY"), keytype = "ENTREZID"
)
go_entrez <- dplyr::rename(go_entrez, gene_key = ENTREZID)
go_entrez <- dplyr::inner_join(go_entrez, key_map_entrez, by = "gene_key", relationship = "many-to-many")

go_symbol <- AnnotationDbi::select(
  org.Hs.eg.db, keys = key_map_symbol$gene_key,
  columns = c("GO", "ONTOLOGY"), keytype = "SYMBOL"
)
go_symbol <- dplyr::rename(go_symbol, gene_key = SYMBOL)
go_symbol <- dplyr::inner_join(go_symbol, key_map_symbol, by = "gene_key", relationship = "many-to-many")

# Several phospho-sites of the same gene (e.g. SRC_PY527/SRC_PY416) share one
# gene_key, so each of that key's GO rows must join to *every* matching
# protein -- a plain match() would silently keep only the first and lose the
# others. Restricting to Biological Process (BP) keeps the resulting groups
# comparable to one another (mixing BP/CC/MF terms would not).
go_annotations <- dplyr::bind_rows(go_entrez, go_symbol)
go_annotations <- dplyr::filter(go_annotations, !is.na(GO), ONTOLOGY == "BP")
go_annotations <- dplyr::distinct(go_annotations, protein, GO) # dedup evidence-code rows

# Each protein typically has dozens of BP terms; keeping the first one
# returned is essentially arbitrary and fragments proteins into many
# near-singleton groups. Instead, for each protein, keep the GO term shared
# by the *most* other proteins in this dataset -- giving fewer, non-trivial
# groups that are more meaningful to compare against the Normal-Block
# clustering.
term_freq <- dplyr::count(go_annotations, GO, name = "n_proteins")
go_best   <- dplyr::left_join(go_annotations, term_freq, by = "GO")
go_best   <- dplyr::group_by(go_best, protein)
go_best   <- dplyr::slice_max(go_best, n_proteins, n = 1, with_ties = FALSE)
go_best   <- dplyr::ungroup(go_best)

# Final membership vector, in the same order as gene_annotation's rows
# (== colnames(brca_rppa$expr)): every protein gets its most-shared BP term,
# or "unknown" if none was found. Matched by name rather than relying on row
# order, in case gene_annotation is ever reordered upstream.
go_bp_term <- setNames(rep("unknown", nrow(ga)), ga$protein)
go_bp_term[match(go_best$protein, names(go_bp_term))] <- AnnotationDbi::Term(go_best$GO)
stopifnot(all(names(go_bp_term) == brca_rppa$gene_annotation$protein))

brca_rppa$gene_annotation$go_bp_term <- unname(go_bp_term)

save(brca_rppa, file = "data/brca_rppa.rda", compress = "xz", version = 3)

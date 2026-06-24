# computes ARI between two clusterings; test-only helper (hence the
# dependency on aricode, which only needs to be in Suggests)
matching_group_scores <- function(groups1, groups2) {
  aricode::ARI(groups1, groups2)
}

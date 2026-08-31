This is a resubmission. It addresses the feedback from the previous round:

* Explained the "EM" acronym in the DESCRIPTION text.
* Replaced `\dontrun{}` with runnable example code (the wrapped example ran
  in well under 5s, so it was unwrapped entirely rather than switched to
  `\donttest{}`).
* Replaced a few `cat()`/`print()` calls that could not be suppressed with
  `message()`, or gated them behind `verbose` (in
  `NormalBlockVarCollectionSparsity`, `NormalBlockVarCollectionClustersSparsity`
  and `SelectionNClusters`); calls inside `print`/`summary` methods were left
  as is.

## Tested environments

* tested locally on Ubuntu Linux 24.04, R-release, GCC

* tested remotely with github-action

- Linux ubuntu 24.04, R-release (github-action)
- Linux ubuntu 24.04, R-oldrel-1 (github-action)
- Linux ubuntu 24.04, R-devel (github-action)
- Windows Server 2022, R-release, 64 bit (github-action)
- macOS latest, R-release (github-action)

* tested remotely with win-builder (R-oldrelease, R-release, R-devel)
  1 NOTE for new submission
  + False positive for mispelled words in DESCRIPTION (authors' names)

## Local R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission -- expected for a first release.

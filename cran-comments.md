## Submission

normalblockr 1.0.0 is a major release. It accompanies a submission to the
Journal of Statistical Software describing the package, whose replication
code requires this version.

Main changes since 0.3.0 (see NEWS.md):

* model-selection criteria aligned with the literature (one BIC penalty per
  group of parameters; ICL completing the partition only);
* new penalty-free selection tools: held-out log-likelihood and K-fold
  cross-validation (`normal_block_cv()`); StARS for both model families;
* defaults unified across the two model families (initialization, rescaling
  of the data, local search after fitting a collection);
* extractors (`coef()`, `sigma()`, `fitted()`, `predict()`) consistently on
  the original scale of the data, with names.

## Test environments

* local: Ubuntu 24.04, R-release, GCC (`R CMD check --as-cran`)
* win-builder: R-release (4.6.2), R-oldrelease (4.5.3), R-devel (2026-11-17 r90926): Status OK
* GitHub Actions: Ubuntu 26.04 (R-release, R-devel, R-oldrel-1), macOS 26
  arm64 (R-release), Windows Server 2025 (R-release)

## R CMD check results

0 errors | 0 warnings | 0 notes attributable to the package.

The local `--as-cran` run reports notes that are artifacts of this machine
rather than of the package: the non-portable `-mno-omit-leaf-frame-pointer`
that Debian/Ubuntu's own `r-base` puts in `Makeconf` (the package's
`src/Makevars` sets only `CXX_STD`, `-DARMA_WARN_LEVEL=1` and the BLAS/LAPACK
libs), HTML Tidy not being installed here, and -- intermittently, when GitHub
rate-limits us -- a 429 on the README's last-commit badge URL, which resolves
normally otherwise.

## Previous submission

The 0.3.0 ERROR on r-devel-linux-x86_64-fedora-gcc (`could not find function
"initRefFields"`, through `sbm::estimateSimpleSBM()` and blockmodels'
Reference Classes) came from an r-devel regression in `methods`, since
resolved: all flavors are OK for 0.3.0.

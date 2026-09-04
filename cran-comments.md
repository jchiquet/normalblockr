This is a feature release. The main addition is a second, complementary model
family (mean-block: variables clustered by their regression response to the
covariates rather than by their covariance), with its zero-inflated
counterparts and collections.

The package's own graphical-lasso solver replaces the `glassoFast`
dependency, which has therefore been dropped from Imports. The
`RhpcBLASctl` suggestion has been dropped too, along with the control
argument that used it.

## Tested environments

* locally on Ubuntu Linux 24.04, R-release, GCC (`R CMD check --as-cran`)

* remotely with github-actions:
  - Linux ubuntu 24.04, R-release
  - Linux ubuntu 24.04, R-devel
  - Linux ubuntu 24.04, R-oldrel-1
  - macOS, R-release
  - Windows Server, R-release

## R CMD check results

0 errors | 0 warnings | 3 notes

* "Days since last update: 1" -- please see the note below.

* "Compilation used the following non-portable flag(s):
  '-mno-omit-leaf-frame-pointer'". This flag comes from the local R
  installation's own `CXXFLAGS` (Ubuntu's r-base build), not from the
  package: `src/Makevars` sets only `CXX_STD`, `PKG_CPPFLAGS` and
  `PKG_LIBS`.

* "Skipping checking HTML validation: no command 'tidy' found" -- local
  toolchain only.

The URL flagged as possibly invalid
(`https://github.com/jchiquet/normalblockr/commits/master`, status 429) is
GitHub rate-limiting the checker; the URL resolves normally in a browser.

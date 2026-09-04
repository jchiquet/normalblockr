This release fixes the ERROR currently shown on
`r-devel-linux-x86_64-fedora-gcc`, where re-building
`breast-cancer-proteomics.Rmd` aborts with

    *** caught segfault ***
    address 0x580, cause 'memory not mapped'

It is submitted sooner than the usual interval for that reason; please let
us know if you would rather we wait.

## The fix

The crash is in the vignette's graphical-lasso section. Until now that step
called back into R -- `glassoFast::glassoFast()` -- from inside the C++
(V)EM loop, once per M-step, i.e. thousands of R re-entries from within a
single `.Call()` for one penalty path. One memory-safety bug on that
boundary had already been found and fixed in an earlier version (an R object
cached in a C++ `static`, outliving R's protection discipline), and the
crash signature matched what our own CI showed intermittently on Linux
across several R versions.

This release removes that boundary rather than working around it: the
graphical lasso is now implemented in the package's own C++
(`src/graphical_lasso.h`), and `glassoFast` is no longer a dependency. There
is no longer any call from C++ back into R inside the recursion.

We should be straightforward about the limits of this: we were never able to
reproduce the crash locally (it appeared on roughly one CI run in ten,
across R versions and platforms, which is the expected profile of heap
corruption whose manifestation depends on memory layout), so we cannot
demonstrate the fix directly on a failing case. The argument is that the
suspected mechanism is gone, not that we caught it in the act. As
supporting evidence, the vignette's full pipeline has been run 40 times
under `MALLOC_CHECK_=3` on the new code without incident.

The new solver was validated against `glassoFast` while both were installed:
477 random problems agree to a worst relative difference of 2.2e-13, and on
40 covariance/penalty pairs captured from a real penalty path it also takes
the identical number of inner coordinate-descent passes.

## Other changes

The release also adds a second model family (variables clustered by their
regression response to the covariates rather than by their covariance). See
NEWS.md.

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

* "Days since last update" -- see above; this submission answers the check
  ERROR on r-devel-linux-x86_64-fedora-gcc.

* "Compilation used the following non-portable flag(s):
  '-mno-omit-leaf-frame-pointer'". This flag comes from the local R
  installation's own `CXXFLAGS` (Ubuntu's r-base build), not from the
  package: `src/Makevars` sets only `CXX_STD`, `PKG_CPPFLAGS` and
  `PKG_LIBS`.

* "Skipping checking HTML validation: no command 'tidy' found" -- local
  toolchain only.

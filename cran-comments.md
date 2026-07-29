This is a new package.

## Tested environments

* tested locally on Ubuntu Linux 24.04, R-release, GCC

* tested remotely with github-action

- Linux ubuntu 24.04, R-release (github-action)
- Linux ubuntu 24.04, R-oldrel-1 (github-action)
- Linux ubuntu 24.04, R-devel (github-action)
- Windows Server 2022, R-release, 64 bit (github-action)
- macOS latest, R-release (github-action)

* tested remotely with win-builder (R-oldrelease, R-release, R-devel) -- submitted, results pending

## Local R CMD check results

0 errors | 0 warnings | 3 notes

* checking CRAN incoming feasibility ... NOTE
  New submission -- expected for a first release.

* checking compilation flags used ... NOTE
  Compilation used the non-portable flag '-mno-omit-leaf-frame-pointer'.
  This flag comes from the local toolchain's default R Makeconf, not from
  this package's own src/Makevars (which sets none of the compilation
  flags), so it is not under this package's control.

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found.
  Local environment limitation (HTML Tidy not installed); not expected on
  CRAN's own check machines.

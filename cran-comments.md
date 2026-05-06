## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local Windows 11, R 4.5.3
* GitHub Actions: windows-latest (R release)
* GitHub Actions: macos-latest (R release)
* GitHub Actions: ubuntu-latest (R devel, release, oldrel-1, oldrel-2)

## Notes

* This is a new submission.
* The package uses compiled C++ code (Rcpp/RcppParallel) with vendored
  xtensor headers. The vendored source is documented in
  `src/vendor/VENDORED_XTENSOR.md`.
* SystemRequirements: GNU make, C++17.

# xsdm methodology — knowledge base

This directory captures the engineering recipes that turned `xsdm` from an
R-only package into a hybrid R/C++ pipeline that scales to continental
rasters and runs the optimizer without R-callback overhead. It is written
as **reusable reference material** — every recipe should transplant cleanly
into a new R package that needs the same capability.

| # | Document | What you get out of it |
|---|----------|------------------------|
| 1 | [01-xtensor-cpp-arrays.md](01-xtensor-cpp-arrays.md)        | How to add `xtensor` to an Rcpp package and write zero-copy C++ kernels that accept raw R pointers and return R-friendly buffers. |
| 2 | [02-log_prob_detect_cpp.md](02-log_prob_detect_cpp.md)      | How to collapse a multi-function R call chain into a single xtensor + RcppParallel C++ kernel, with parity tests against the R reference. |
| 3 | [03-terra-blockwise-io.md](03-terra-blockwise-io.md)        | How to stream a list of `terra::SpatRaster` time series block-by-block through the C++ kernel using raw `double*` and stay memory-bounded. |
| 4 | [04-ucminfcpp-xptr.md](04-ucminfcpp-xptr.md)                | How to wrap a pure-C++ objective in `Rcpp::XPtr<ucminf::ObjFun>`, eliminate R-callback overhead inside the optimizer, and harden the closure against fatal errors. |
| 5 | [05-optimization-convergence-and-equivalence-class.md](05-optimization-convergence-and-equivalence-class.md) | What "good convergence" actually means when the model has a discrete equivalence class, why loose `grtol`/`xtol` defaults silently sabotage multi-start fits, and how to canonicalize parameter vectors so multi-start results become directly comparable. |

Each document follows the same structure:

1. **Why** — what problem the technique solves and what shape of code it
   replaces.
2. **Prerequisites** — package metadata, compiler flags, system deps.
3. **Step-by-step recipe** — numbered steps with the exact files to add or
   modify and short inline code snippets.
4. **Tests / parity checks** — how to verify the new path matches the
   reference at numerical tolerance.
5. **Pitfalls** — concrete bugs we hit during development that future
   reimplementations should avoid.

The reference implementation lives in this repository; every recipe links
to the canonical files and line ranges it is distilled from.

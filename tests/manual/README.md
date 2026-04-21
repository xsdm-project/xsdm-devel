# Manual visual tests for `interpret_parameters()`

These scripts are **not** part of automated `R CMD check`. They require human
judgment to verify that plots look correct.

## Running the tests

Open R in the package root directory (or set working directory to the package
root) and run:

```r
source("tests/manual/test_interpret_parameters_visually.R")
```

All tests will run sequentially. Each plot will appear on the default graphics
device (e.g., RStudio Plots pane). Examine each plot carefully.

## What to verify
### Legacy single-panel (no env_dat/occ)
- Univariate curves should be smooth, peak at mu, and go down to negative values.
- Bivariate contours should be ellipses (or rotated ellipses) centred at mu.
- Changing plot_lims should only change axis ranges, not the shape.

### Two-panel (presence vs non‑detection)
- Left panel (presences): points should be overlaid on the growth function.
- Right panel (non‑detections): same curve, different points.
- For breadth = 1 (default), all points should lie inside the axis limits with
some breathing room (margin).
- For breadth = 0, the axes collapse to a tiny region around mu; you should
see essentially no data points (only a dot at the centre). This is expected.

### auto_plot_lims_ smoke tests
- These run automatically and print status messages to the console.

- They do not require visual inspection; they check that limits cover the
data and that the pinprick is tiny. Look for “OK” messages.

### Error cases
- The script deliberately triggers errors – these should produce error messages
but not stop the entire script. Verify that the error messages are
sensible and that the function fails gracefully.

## Optional: Save plots to PDF
If you want to review the plots later (or share them), uncomment the PDF block
at the top of the script:

```r
save_pdf <- TRUE
if (save_pdf) {
  pdf_file <- tempfile(pattern = "interpret_parameters_", fileext = ".pdf")
  pdf(pdf_file, width = 10, height = 6)
  on.exit(dev.off(), add = TRUE)
  message("Saving plots to: ", pdf_file)
}
```

The PDF will be created in a temporary directory.

## Troubleshooting
- If you get could not find function "auto_plot_lims_", make sure you are
running the script from the package root so that xsdm is loaded and its
internal functions are accessible via xsdm:::.

- If the example data objects (example_1_env_array, etc.) are missing, load
them first with data(example_1_env_array, package = "xsdm") or replace
with your own data.

## When to run these tests
- After modifying interpret_parameters() or auto_plot_lims_().
- Before submitting a package update to CRAN (to catch any unintended visual
regressions).
- When adding new features (e.g., new ... arguments, different colour palettes).

## Contact
- If a plot looks obviously wrong (e.g., missing points, wrong axis scaling),
please file an issue with a screenshot.

### Why this is useful

- **Onboarding** – New contributors can quickly understand what the script does.
- **Reproducibility** – Instructions avoid guesswork (e.g., “Do I need to set a
seed?” No, because the jitter is just for visual reference).
- **Documentation of expected behaviour** – Helps catch regressions (“The left
panel used to show points inside the curve – now they’re outside.”)
- **Separation of concerns** – Clearly distinguishes manual tests from automated
unit tests.



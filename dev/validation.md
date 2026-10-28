# Validation record

2026-10-28, macOS aarch64, R 4.6.1. This is local development validation.

## Completed

- `air format R tests`.
- `devtools::document()` regenerated NAMESPACE and Rd documentation.
- `devtools::test()`: **97 passed, 0 failed, 0 warnings, 0 skipped**.
- README.Rmd rendered with `rmarkdown::github_document(html_preview = FALSE)`
  after `devtools::load_all()`; all evaluated examples completed.
- `R CMD build glyfp`: **OK**, including vignette construction.
- `R CMD check glyfp_0.0.0.9000.tar.gz`: **Status: OK**. Installation, examples,
  installed-package tests, vignette rebuild and PDF manual passed.
- Both workflow YAML files parsed. GitHub's `on` key is YAML 1.1 boolean TRUE in
  R's yaml parser; this is expected and not a GitHub Actions syntax failure.
- SVG rendered with librsvg; the final 424 x 489 transparent PNG was visually
  inspected. Earlier ImageMagick SVG rendering had missing strokes/text and was
  replaced by the verified librsvg export. SVG remains the editable source.
- Independent reviewer found no functional issues. It independently enumerated
  connected vertex sets for a three-residue chain, a three-residue branch, and
  a five-residue N-glycan core; dictionary/counts agreed. The N-glycan example
  had 11 distinct features and 13 occurrences. Branch-order spelling was
  invariant. This small-tree oracle does not establish large-tree performance.
- Added an explicit default-size regression: mono-, di- and trisaccharides
  produce counts (1,0,0), (2,1,0), and (3,2,1) in one shared dictionary.
- `dev/validate-sweetnet.py`: actual official pretrained legacy SweetNet CPU
  inference passed (4 x 128 finite values). Repeated calls and duplicate rows
  were exact; input reordering was exact in the tested batches. Batch size
  1 versus 4 had maximum absolute difference 4.0531158447265625e-06, within
  rtol/atol 1e-5. Strict checkpoint loading and SHA-256 validation succeeded.
  Exact package versions, checkpoint/vocabulary hashes and readout are recorded
  in methods.md. This is a feasibility probe outside the exported R API.

Tested dependency versions: glyrepr 1.2.0, glymotif 1.0.0, glyparse 0.8.1,
igraph 2.3.4, testthat 3.3.2, roxygen2 7.3.3.

## Evidence outside the package source

All paths below are relative to the parent task directory:

- `build.log`: final package build output.
- `check.log`: final check console output.
- `glyfp.Rcheck/00check.log`: final check summary, Status: OK.
- `glyfp.Rcheck/tests/testthat.Rout`: installed-package result, PASS 97.
- `glyfp_0.0.0.9000.tar.gz`: checked source package.
- `check-as-cran.log`: separate unsuccessful online incoming-check attempt.

## Not completed / environment limits

- `R CMD check --as-cran` stopped during CRAN incoming feasibility because the
  sandbox could not connect to the CRAN PACKAGES/PACKAGES.in endpoints (and
  Bioconductor indices). It did not finish and is not reported as passing.
  The normal full R CMD check completed afterwards; repository-index connection
  messages did not cause a final NOTE, WARNING or ERROR.
- `devtools::build_readme()` could not open a pak metadata-cache lock outside the
  writable workspace. No global cache was modified. The README was successfully
  rendered with already-installed dependencies via rmarkdown. Disabling its
  optional HTML preview avoided fetching the external lifecycle badge; the
  generated Markdown retains the badge URL.
- `covr` is not installed, so local numerical coverage was not computed.
  CI contains a coverage job and artifact upload, but no remote CI has run.
- R torch 0.17.0 is installed but `torch::torch_is_installed()` is FALSE. Its
  runtime was not installed and native R numerical parity was not checked.
  Python CPU SweetNet inference succeeded; current GIFFLAR checkpoints and MIT
  licensing were verified, but GIFFLAR inference was not run. The R-facing
  input-normalization/unknown-token contract for an embedding method remains
  integration work; no random/untrained embedding is exposed.

No remote repository, PR, release, publication or deployment was created. A local
Git repository was initialized without commits or a remote. Existing ecosystem
repositories were inspected read-only.

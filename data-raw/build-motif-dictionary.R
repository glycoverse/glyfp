# Run from the glyfp repository root:
# Rscript data-raw/build-motif-dictionary.R /path/to/glydb/data/glydb_data.rda
# Source: glydb commit cafa287eac546a9899977a9e8c777af77cfd2ac3,
# GlyGen reviewed release 2.12.1. The source artifact is hash pinned below.
stopifnot(
  requireNamespace("glyrepr"),
  requireNamespace("glymotif"),
  requireNamespace("digest")
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Provide the pinned glydb_data.rda source file.")
}
source_file <- args[[1L]]
source_sha256 <- "754e00cf6c19307a774e126b621e8103c34be868ecdf3e0d11cf1a8e3fe53796"
stopifnot(identical(
  digest::digest(file = source_file, algo = "sha256"),
  source_sha256
))
source("R/structure-checks.R", local = TRUE)
data_env <- new.env(parent = emptyenv())
load(source_file, envir = data_env)
data <- data_env$glydb_data
stopifnot(nrow(data) == 49728L)
graphs <- as.list(data$glycan_structure)
issues <- lapply(graphs, function(g) {
  if (is.null(g)) "missing_structure" else structure_issues(g)
})
keep <- lengths(issues) == 0L
structures <- unique(data$glycan_structure[keep])
structures <- structures[order(
  enc2utf8(as.character(structures)),
  method = "radix"
)]
cat(
  "Selected",
  sum(keep),
  "accessions and",
  length(structures),
  "unique complete structures\n"
)
features <- character()
# Bounded batches avoid retaining all extracted graph objects simultaneously.
for (start in seq.int(1L, length(structures), by = 250L)) {
  end <- min(start + 249L, length(structures))
  extracted <- glymotif::extract_motif(structures[start:end], max_size = 7L)
  features <- sort(
    unique(c(features, enc2utf8(as.character(extracted)))),
    method = "radix"
  )
  cat("Processed", end, "structures;", length(features), "motifs\n")
}
# Only complete, concrete source structures enter extraction; extracted motifs
# inherit their reducing-end anomer from the original root or attachment bond.
motifs <- glyrepr::as_glycan_structure(features)
stopifnot(all(vapply(
  as.list(motifs),
  function(g) !length(structure_issues(g)),
  TRUE
)))
sizes <- vapply(as.list(motifs), igraph::vcount, 1)
resource <- data.frame(size = as.integer(sizes), motif = features)
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
path <- "inst/extdata/motif-dictionary.tsv"
write.table(
  resource,
  path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  fileEncoding = "UTF-8"
)
metadata <- list(
  id = "glydb-2.12.1-complete-v1",
  max_size = 7L,
  source_package_version = "0.7.2.9000",
  source_commit = "cafa287eac546a9899977a9e8c777af77cfd2ac3",
  source_release = "GlyGen 2.12.1",
  source_sha256 = source_sha256,
  source_url = "https://data.glygen.org/ln2data/releases/data/v-2.12.1/reviewed/",
  source_rows = nrow(data),
  included_accessions = sum(keep),
  unique_structures = length(structures),
  excluded_accessions = sum(!keep),
  exclusion_counts = as.list(table(unlist(issues))),
  sizes = as.list(table(sizes)),
  feature_count = length(features),
  dictionary_sha256 = digest::digest(file = path, algo = "sha256"),
  included_accessions_sha256 = digest::digest(
    paste(sort(data$glytoucan_ac[keep], method = "radix"), collapse = "\n"),
    algo = "sha256",
    serialize = FALSE
  ),
  r_version = as.character(getRversion()),
  package_versions = vapply(
    c("glyrepr", "glymotif", "glyparse"),
    function(x) as.character(utils::packageVersion(x)),
    ""
  )
)
writeLines(
  sub("[[:space:]]+$", "", capture.output(dput(metadata))),
  "inst/extdata/motif-dictionary-metadata.R"
)
cat("Dictionary SHA256:", metadata$dictionary_sha256, "\n")

# Built-in motif dictionary

`build-motif-dictionary.R` builds `glydb-2.12.1-complete-v1` from the real
glydb source artifact, without a runtime glydb dependency or network access.
Run from the package root:

```sh
Rscript data-raw/build-motif-dictionary.R /path/to/glydb/data/glydb_data.rda
```

Source: glydb 0.7.2.9000, commit `cafa287eac546a9899977a9e8c777af77cfd2ac3`.
The script verifies SHA-256
`754e00cf6c19307a774e126b621e8103c34be868ecdf3e0d11cf1a8e3fe53796`.
This snapshot contains 49,728 parsed structures from GlyGen reviewed 2.12.1,
downloaded by glydb on 2026-10-13. Its upstream preparation started with 59,994
WURCS records and excluded 10,266 parsing failures before this dictionary step.
The dictionary inherits glydb/glyparse representation coverage; it does not
recover chemistry the upstream parser could not represent.

Selection uses actual glyrepr graph fields, not the presence of a question mark
in an input string or glydb's coarse structure-level label. Only concrete
residues, fully specified bonds and reducing ends, and positioned substituents
are included. Floating components/substituents and alditols are excluded.
There are 8,416 included accessions, 8,211 unique structures and 41,312 excluded
accessions. Reasons overlap: uncertain reducing end 31,654; uncertain linkage
23,836; floating parts 6,218; unsupported alditol 4,443; nonconcrete residue 2,681;
uncertain substituent 637; floating substituent 27.

The generator sorts unique canonical structures by UTF-8 bytes, enumerates all
connected motifs up to size seven with glymotif, deduplicates and sorts motifs
by the same rule. Reducing-end information is retained, including when a
fragment inherits its root anomer from an internal glycosidic bond. Counts by
exact size 1–7 are 134, 906, 2,133, 3,110, 4,355, 5,849 and 7,637 (24,124 total).
The default sizes 1–3 select 3,173 dimensions. TSV SHA-256:
`141b7d72541cd6331c385da04afba7ebd0962ec3b3b6fbb2715ca0c6cdd3445a`.
`inst/extdata/motif-dictionary-metadata.R` also records included-accession hash,
source identity, exclusion counts and generator dependency versions. Rebuild
with those versions for byte-for-byte comparison. Changes require a new
dictionary identity and recomputation of compared fingerprints.

GlyGen licenses database sets under CC BY 4.0; glydb's package code is MIT.
See `inst/COPYRIGHTS` for attribution, links and the description of changes.
Source: https://data.glygen.org/ln2data/releases/data/v-2.12.1/reviewed/
License: https://www.glygen.org/license/

The 2012 paper *Glycan Fingerprints: Calculating Diversity in Glycan Libraries*
(DOI 10.1021/cb300003z) motivates connected 1–7 residue subtrees. Its method hashes
and folds binary strings; it does not publish this explicit dictionary. The
fixed glydb dictionary is a separate engineering design. No dataset-dependent
occupancy correction or modified Tanimoto score is implemented here.

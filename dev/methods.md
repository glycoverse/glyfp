# Method review and first-version decisions

Reviewed 2026-10-28. This is a source-backed implementation review, not a claim
that published model accuracy has been reproduced.

## Representative methods

| Family | Source / implementation | Feature and counting unit | Alignment and limitations |
|---|---|---|---|
| Binary subtree fingerprints | Rademacher & Paulson (2012), DOI [10.1021/cb300003z](https://doi.org/10.1021/cb300003z), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC3356451/) | Presence of connected residue fragments, including branches, up to seven residues; the published method hashes fragments to bits and folds the bit string. | Fixed hashing defines coordinates but allows collisions. Residue and linkage distinctions matter; not an atomic/3D representation. Original scripts are described as available on request, not downloaded or reproduced here. |
| Explicit binary dictionary | Cao et al., GlyMDB (2020), DOI [10.1093/bioinformatics/btz934](https://doi.org/10.1093/bioinformatics/btz934), [full text](https://europepmc.org/articles/PMC7178435) | Enumerates connected substructures and uses a unique-fragment dictionary. A glycan-by-fragment binary matrix supports supervised motif discovery. | An explicit shared dictionary avoids hash collisions. GlyMDB is a web tool, not an R fingerprint package. Its classifier and recursive feature elimination are outside glyfp scope. |
| Substructure counts | [glymotif count_motifs source](https://github.com/glycoverse/glymotif/blob/main/R/count-motif.R), [extract_motif source](https://github.com/glycoverse/glymotif/blob/main/R/extract-motif.R) | Unique matched residue sets, allowing overlaps. Swapping identical branches within one occurrence does not add a count. | Uses the same explicit dictionary as binary fingerprints. Strict mode supports directional generic/ambiguous motif matching; it is not exact-label equality. This is a count extension using an existing implementation, not a claim that the two papers used counts. |
| Learned graph embedding | Joeres & Bojar (2024), [Higher-Order Message Passing for Glycan Representation Learning](https://arxiv.org/html/2409.13467v1), [GIFFLAR](https://github.com/BojarLab/GIFFLAR) | Higher-order message passing over atoms, bonds and monosaccharides, followed by graph pooling. Features are learned latent coordinates, not fragment counts. | Fixed architecture, preprocessing, checkpoint, readout and inference mode are required for comparable vectors. Atomic graphs require sufficiently specified input; the paper filters structures that cannot be translated through GlyLES. |
| Lightweight learned graph embedding | Burkholz, Quackenbush & Bojar (2021), [SweetNet](https://doi.org/10.1016/j.celrep.2021.109251), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9208909/), [glycowork implementation](https://github.com/BojarLab/glycowork) | Monosaccharide and linkage tokens form graph nodes, followed by three graph convolutions and graph pooling. The pretrained species model exposes a 128-coordinate representation. | Token chemistry and connectivity are represented, not atom coordinates/3D conformations. A fixed token vocabulary, architecture, weight file and readout define coordinates; different SweetNet generations are not interchangeable. |

The 2012 paper's example configuration uses 1024 bits folded once. glyfp does
**not** reproduce that hashed/folded representation and does not use its name for
an incompatible implementation. GlyMDB demonstrates the explicit dictionary
alternative: the reported CFG 5.1 example has 541 glycans and 14,973 fragments.
Those numbers describe that source dataset, not a built-in glyfp dictionary.

## Implemented methods

`substructure_binary` and `substructure_count` clearly identify both the unit and
encoding. Binary is the default because presence fingerprints are the starting
point in both cited glycan applications. Both methods use exactly the same counts;
binary applies `count > 0`. `max_size = 3` matches glymotif's public default and
covers the smallest branching motifs while keeping exploratory fitting modest.
It includes **every connected fragment size from one through three**, not only
exact-size trisaccharides. Mono- and disaccharides therefore retain features.
The first API accepts 1–7 residues, including the heptasaccharide setting in the
2012 work. This upper bound is a practical first-version scope, not a theorem
about adequate glycan coverage. Large branching trees can still be expensive.

Production extraction and search are delegated to glymotif's exported APIs:

1. Normalize characters with `glyparse::auto_parse(on_failure = "error")`, matching
   glymotif's input policy; preserve `glyrepr_structure` vectors and names.
2. Remove missing entries and deduplicate before `glymotif::extract_motif()`.
3. Store unique canonical IUPAC strings sorted using radix ordering. No arbitrary
   graph IDs, hashes, random seeds, or training order become coordinates.
4. Rebuild motif objects from these strings and call `glymotif::count_motifs()`
   with substructure alignment, linkages included, strict substituents and strict
   matching mode. Deduplicate targets before matching and expand afterwards.

Canonical keys preserve the extracted motif half-linkage. A beta-linked residue
and an alpha-linked residue can therefore yield different single-residue motifs.
Generic residues in motifs match compatible specific residues; unknown motif
positions are constraints that can match several concrete positions. Columns may
therefore overlap semantically. Unknown-feature policies refer to **new canonical
extracted motifs**, even if a broader known motif matches the same structure.
There is no claim that all motif columns are disjoint or all counts sum to a
glycan's number of fragments.

Substituents and linkage labels are retained, but omitted upstream chemistry
cannot be recovered. Atomic stereochemistry beyond those labels, aglycones and
3D conformation are not represented. Floating metadata and alditols are rejected
in this first version: allowing extraction to discard such information would be
misleading. The local glymotif source did not expose an alditol-sensitive matching
contract, so this is deliberately tested as an unsupported boundary.

## Batch, missingness and persistence contract

Top-level calls fit on the supplied batch and return only a vector list. Equal
vector lengths from separate calls do not imply comparable coordinates. For
comparisons, combine all glycans once or fit once and reuse one encoder. For model
evaluation, fit on training data only.

The encoder stores `features`, `method`, `params`, schema version and exact
versions of glyrepr/glymotif/glyparse. No graph pointers are serialized. Versions
must match at transform time; to upgrade, refit and recompute all batches together.
This conservative requirement avoids silently accepting upstream changes to
parsing, canonicalization or matching. Dictionary order and method are never
modified during transform. Do not edit encoder internals by hand.

`NA` entries do not add features and yield all-NA vectors. Empty/all-NA fitting
creates an empty dictionary; missingness cannot be encoded in length-zero vectors,
so callers must retain the original mask. Unknown features are ignored, warned
about once per batch, or rejected according to the stored policy. They never
append coordinates. No automatic growing dictionary is provided.

## Embeddings: R torch versus a Python bridge

The original GIFFLAR [requirements](https://github.com/BojarLab/GIFFLAR/blob/main/requirements.txt)
include torch, torch_geometric >= 2.6.1, RDKit, GlyLES, Lightning, and other Python
packages. The inspected [DownstreamGGIN implementation](https://github.com/BojarLab/GIFFLAR/blob/main/gifflar/model/downstream.py)
returns `graph_embed` separately from predictions. Its output size is configured
by the model, so a universal embedding dimension should not be invented.

| Route | What works in principle | Work still needed |
|---|---|---|
| Native R torch | CPU tensors and neural modules using LibTorch; no Python required once the architecture and preprocessing are ported. | Reproduce heterogeneous graph construction, PyG operations, pooling, tensor names/dtypes, and checkpoint mapping. Verify numerical parity against the original implementation. |
| Exported scripted model | `torch::jit_load()` loads serialized scripted modules/functions. | Produce and validate a compatible export; PyG graph containers/custom operations and chemical preprocessing are not automatically covered. No GIFFLAR export was verified. |
| Optional Python bridge | Keep the authors' PyTorch/PyG/GlyLES/RDKit pipeline and expose graph embeddings through a small adapter. | Pin Python/dependency versions, use a verified checkpoint, validate CPU execution and batch invariance, and keep this optional rather than burdening binary/count users. |

R torch's official [Python-to-R guide](https://torch.mlverse.org/docs/articles/python-to-r)
shows rebuilding the architecture with matching parameter names before loading a
Python state dictionary. Its [installation guide](https://torch.mlverse.org/docs/articles/installation.html)
documents CPU builds and additional LibTorch/LibLantern installation; the
[`jit_load` reference](https://torch.mlverse.org/docs/reference/jit_load) documents
scripted artifacts. These are viable mechanisms, not proof that a GIFFLAR Python
checkpoint is directly loadable in R.

R torch 0.17.0 is installed locally, but `torch::torch_is_installed()` returns
FALSE. This only means the R runtime is not installed; it does not prevent
Python CPU inference or establish that the model is unavailable.

### Official checkpoints and licenses

The standalone GIFFLAR repository at commit
`7a006ff2dac6449a647532d752f519f89e8a7183` contains no root LICENSE or checked-in
checkpoint, and its README only documents setup. That observation must not be
generalized to the newer official glycowork distribution:

- [glycowork](https://github.com/BojarLab/glycowork), reviewed at commit
  `199efdde6bc137dbbcd900f84f927fdbbca52eae` (version 1.11.0), is MIT licensed.
  Its model loader retrieves weights from the author's Hugging Face repository.
- [DBojar/glycowork_models](https://huggingface.co/DBojar/glycowork_models)
  is public, ungated, and declares `license: mit` in its model card. Reviewed
  repository revision: `52187c054b9c988ddb2a5a5af9ee64b51004fd1e`.
- The [file tree](https://huggingface.co/DBojar/glycowork_models/tree/main)
  provides legacy SweetNet weights (3.29 MB), SweetNet v2 weights (3.82 MB), and
  GIFFLAR species weights (7.78 MB). These are small checkpoints, though the
  Python runtime and chemical dependencies are substantially larger.
- Current glycowork's GIFFLAR loader explicitly loads its trained checkpoint.
  The older local glycowork 1.9.0 instead warns and returns random GIFFLAR weights
  when `trained=True`: this path must never be used as a pretrained embedding.
- SweetNet v2 adds residual/norm graph layers and mean+max pooling; legacy
  SweetNet uses mean pooling. Matching tensor shapes alone is not enough to
  interchange model generations or vocabularies.

No model code or weight file is bundled in the R package. The paper's license
is not used as a substitute for the separately declared code/weight licenses.

### CPU feasibility probe and next integration boundary

`dev/validate-sweetnet.py` uses the authors' actual legacy SweetNet architecture
and official pretrained checkpoint, verifies its SHA-256, loads the complete
state dictionary strictly with `weights_only=True`, and runs in CPU evaluation
mode under `torch.inference_mode()`. It performs no training or downloads.
The probe is separate from the exported R API; binary/count users need no Python.

Checkpoint: `glycowork_sweetnet_species.pt`, pinned Hugging Face revision
`0a3981df888532e52b4c2d910bfb56b5a45ad955`, SHA-256
`4ee19fc0593db0562c16585b556050694144816754750ebfe1876ed21ba0f4be`.

The CPU probe **passed** on four inputs (three distinct glycans, including a
branched glycan and a monosaccharide). Output shape was 4 x 128, all finite;
duplicate rows and repeated calls were bit-identical. Reversing input order and
using batches of two gave zero difference after restoring order. Batches of one
versus four differed by at most `4.0531158447265625e-06`, passing `rtol=1e-5`,
`atol=1e-5`. This is a small functional probe, not accuracy reproduction or a
guarantee of bitwise identity across platforms and BLAS implementations.

Tested environment: Python 3.11.16, glycowork 1.9.0, torch 2.8.0,
torch-geometric 2.6.1, GlyLES 1.2.2, RDKit 2026.9.1, NumPy 2.4.6, pandas 3.0.5,
pydot 1.4.2, antlr4-python3-runtime 4.13.2. The sorted compact JSON token mapping
has SHA-256 `838281860e3bbddc0de91ac635bf1ba746d79bf62891e5b86ca08eff537a974e`.
The readout is `bn2(lin2(...))`, before the second activation/dropout/classifier.
PyTorch was already available; missing PyG/chemical imports were installed into
an isolated task directory. No existing Python environment was modified.
Although SweetNet itself uses residue/linkage graphs, this glycowork release
eagerly imports GIFFLAR's RDKit/GlyLES support even for SweetNet, so installing
only PyTorch/PyG was insufficient. The 3.29 MB weight was downloaded separately
into the task directory, not into package source.

For a production optional `sweetnet` method, the encoder must fix the code and
preprocessing versions, checkpoint hash, token-vocabulary hash, readout and
dimension. It would expose numerical latent-coordinate names instead of a motif
dictionary, while preserving the same list-only top-level contract. Before adding
that method, validate conversion from all supported glyrepr inputs to the model's
IUPAC dialect, including reducing-end half-linkages, unsupported modifications,
missing values and unknown tokens. A general `embedding` name would obscure the
model/version distinction. GIFFLAR remains a second candidate after molecular
graph conversion and CPU inference are checked with its current checkpoint.

One concrete normalization gap was reproduced: glyrepr canonicalizes
`Gal(b1-4)GlcNAc` to `Gal(b1-4)GlcNAc(?1-`. glycowork 1.9.0 graph parsing treats
the final `?1-` as an extra node absent from the pretrained vocabulary; explicit
reducing-end `b1-` is also absent. Passing canonical strings through unchanged
would therefore alter the model input. Silently deleting these constraints
would discard chemistry that the structural fingerprints retain. This needs an
explicit model-specific supported-input policy and text/object parity tests,
so the feasibility probe is not advertised as a completed general R method.

An optional Python bridge is the shortest evidence-backed integration route;
native R torch is reasonable after numerical parity with the original model.
Neither a Python state dictionary nor a `.pt` extension proves that
`torch::jit_load()` can read the file as a scripted model. No scripted GIFFLAR
artifact or native R parity was verified here.

## Ecosystem references and delivery scope

Local glyenzy: version 0.10.0.9000, commit `c276d5b`; inspected README.Rmd,
DESCRIPTION, tests, R-CMD-check and coverage workflows. Its README documents
R-universe/GitHub development installation, without a CRAN badge. The public
[R-universe entry](https://glycoverse.r-universe.dev/glyenzy) likewise provides
R-universe installation and no CRAN link. The attempted CRAN package page could
not be retrieved, so this review does not use that failure alone as proof of
absence. No CRAN/release workflow was copied into glyfp.

Local glyrepr checkout: 1.1.0.9000 at `496e0b0`; actual installed/tested glyrepr:
1.2.0. Installed glymotif: 1.0.0; glyparse: 0.8.1. Source inspection and runtime
checks are distinguished deliberately. The task did not edit any of these
repositories or the independent glyparse release worktree.

The development skeleton includes generated documentation, executable README and
vignette, testthat edition 3, a pkgdown index, and check/coverage workflows. CI is
configured only: no remote repository was created, no workflow was run remotely,
and no release, publication or PR was made. The logo follows the inspected
Glycoverse hexagon geometry and is delivered as editable SVG plus transparent PNG.

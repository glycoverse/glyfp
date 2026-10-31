# glyfp 0.0.0.9000

## Breaking changes

* `glycan_fingerprint()` is removed, so existing calls must migrate to `fp_motif_binary()`, `fp_motif_count()`, `fp_sweetnet()`, `fp_glycangt()` or `fp_gifflar()` according to the former `method`. Calls without `method` become `fp_motif_binary()`. Outputs and defaults are unchanged; `max_size` and `unknown` belong only to motif functions, and `cache_dir` and `download` belong only to SweetNet. For GlycanGT and GIFFLAR, supply the required local `checkpoint` and remove `download = FALSE`.

## Fingerprints

* `fp_circular()` encodes residue-centered WL neighborhoods into fixed hashed binary or count vectors, retaining represented residue chemistry, directed linkage labels and reducing-end state without fitting a dictionary.
* `fp_residue_pair()` encodes endpoint WL environments and glycosidic-edge graph distances into fixed hashed vectors. Distinct unordered pairs count once; single residues yield all-zero vectors. Both new functions document collision behavior and representation limits.

* `fp_motif_count()` and `fp_motif_binary()` accept `motif_dict = "dynamic"` to build a sorted dictionary from the current batch, including motifs outside glydb. Independent batches may define different coordinates; the default `motif_dict = "fixed"` preserves the existing stable dictionary and outputs.

* `fp_motif_binary()` and `fp_motif_count()` encode glyrepr structure vectors and glyparse-supported text as lists of binary or count fingerprints of connected residue substructures. Each method has its own API and exposes only relevant options.
* Substructure methods use the fixed `glydb-2.12.1-complete-v1` dictionary, derived from 8,211 complete structures. The default 3,173 coordinates remain aligned across batches. Unseen motifs can be ignored, warned about or rejected without changing dimensions.
* All methods require complete structures, including the reducing-end alpha/beta state. Motif fingerprints retain that state; pretrained embeddings explicitly project it out because their models do not represent it. Unknown structure fields are rejected before inference.

* `fp_sweetnet()` optionally computes 128-dimensional pretrained CPU embeddings with verified weights and explicit input restrictions. Inference runs entirely in R with the optional torch runtime; no Python is required. The optional torch runtime is never installed automatically.
* `fp_sweetnet()` now defaults to SweetNet v2, with legacy weights available through `checkpoint = "sweetnet-v1"`. Outputs from different model versions must not be combined.
* Named checkpoints download from pinned official sources on first use, verify SHA-256, and reuse a per-user cache. `download = FALSE` supports offline use, and `checkpoint` still accepts verified local files.

* `fp_glycangt()` computes 768-dimensional pretrained graph Transformer embeddings entirely in R. It requires a manually downloaded official local checkpoint and never downloads weights. A fixed model, vocabulary and deterministic per-glycan orthogonal node identifiers keep its coordinates stable across batches.

* `fp_gifflar()` adds 384-dimensional pretrained molecular graph embeddings using R torch and 27 fixed residue templates. It requires a verified local checkpoint, preserves atom, bond, ring and directed linkage information, and does not require Python or a chemistry toolkit.

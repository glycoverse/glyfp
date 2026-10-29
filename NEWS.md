# glyfp 0.0.0.9000

* `glycan_fingerprint()` encodes glyrepr structure vectors and glyparse-supported text as lists of binary or count fingerprints of connected residue substructures.
* `fit_fingerprint()` and `transform_fingerprint()` provide reusable encoders with fixed feature dictionaries and explicit missing-value and unknown-feature policies.

* `method = "sweetnet"` optionally computes 128-dimensional pretrained CPU embeddings with verified weights, a fixed model encoder, and explicit input restrictions. Inference runs entirely in R with the optional torch runtime; no Python is required. The optional torch runtime is never installed automatically.
* `method = "sweetnet"` now defaults to SweetNet v2, with legacy weights available through `checkpoint = "sweetnet-v1"`. Encoders record their model identity; outputs from different model versions must not be combined.
* Named checkpoints download from pinned official sources on first use, verify SHA-256, and reuse a per-user cache. `download = FALSE` supports offline use, and `checkpoint` still accepts verified local files.

* `method = "glycangt"` computes 768-dimensional pretrained graph Transformer embeddings entirely in R. It requires a manually downloaded official local checkpoint and never downloads weights. Encoders fix the model identity, vocabulary and deterministic per-glycan orthogonal node identifiers across batches.

* `method = "gifflar"` adds 384-dimensional pretrained molecular graph embeddings using R torch and 27 fixed residue templates. It requires a verified local checkpoint, preserves atom, bond, ring and directed linkage information, and does not require Python or a chemistry toolkit.

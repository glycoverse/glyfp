# glyfp 0.0.0.9000

* `glycan_fingerprint()` encodes glyrepr structure vectors and glyparse-supported text as lists of binary or count fingerprints of connected residue substructures.
* `glycan_fingerprint()` provides the single public interface for all fingerprint methods, returning a list of vectors.
* Substructure methods use the fixed `glydb-2.12.1-complete-v1` dictionary, derived from 8,211 complete structures. The default 3,173 coordinates remain aligned across batches. Unseen motifs can be ignored, warned about or rejected without changing dimensions.
* All methods require complete structures, including the reducing-end alpha/beta state. Motif fingerprints retain that state; pretrained embeddings explicitly project it out because their models do not represent it. Unknown structure fields are rejected before inference.

* `method = "sweetnet"` optionally computes 128-dimensional pretrained CPU embeddings with verified weights and explicit input restrictions. Inference runs entirely in R with the optional torch runtime; no Python is required. The optional torch runtime is never installed automatically.
* `method = "sweetnet"` now defaults to SweetNet v2, with legacy weights available through `checkpoint = "sweetnet-v1"`. Outputs from different model versions must not be combined.
* Named checkpoints download from pinned official sources on first use, verify SHA-256, and reuse a per-user cache. `download = FALSE` supports offline use, and `checkpoint` still accepts verified local files.

* `method = "glycangt"` computes 768-dimensional pretrained graph Transformer embeddings entirely in R. It requires a manually downloaded official local checkpoint and never downloads weights. A fixed model, vocabulary and deterministic per-glycan orthogonal node identifiers keep its coordinates stable across batches.

* `method = "gifflar"` adds 384-dimensional pretrained molecular graph embeddings using R torch and 27 fixed residue templates. It requires a verified local checkpoint, preserves atom, bond, ring and directed linkage information, and does not require Python or a chemistry toolkit.

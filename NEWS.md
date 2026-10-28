# glyfp 0.0.0.9000

* `glycan_fingerprint()` encodes glyrepr structure vectors and glyparse-supported text as lists of binary or count fingerprints of connected residue substructures.
* `fit_fingerprint()` and `transform_fingerprint()` provide reusable encoders with fixed feature dictionaries and explicit missing-value and unknown-feature policies.

* `method = "sweetnet"` optionally computes 128-dimensional pretrained CPU embeddings with verified weights, a fixed model encoder, and explicit input restrictions. Inference runs entirely in R with the optional torch runtime; no Python is required. Dependencies and model weights are never installed or downloaded automatically.

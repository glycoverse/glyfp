"""CPU feasibility probe, not an exported glyfp method.

Usage: python dev/validate-sweetnet.py /path/glycowork_sweetnet_species.pt
Requires the pinned environment recorded in methods.md. No network requests,
training, or random-weight fallback. Writes a compact validation record to stdout.
"""

import hashlib
import importlib.metadata as metadata
import json
from pathlib import Path
import sys

import numpy as np
import torch

if metadata.version("glycowork") != "1.9.0":
    raise RuntimeError("This probe requires glycowork 1.9.0 (legacy SweetNet).")

from glycowork.glycan_data.loader import lib
from glycowork.ml.models import SweetNet
import glycowork.ml.inference as inference

checkpoint = Path(sys.argv[1])
digest = hashlib.sha256(checkpoint.read_bytes()).hexdigest()
expected = "4ee19fc0593db0562c16585b556050694144816754750ebfe1876ed21ba0f4be"
if digest != expected:
    raise RuntimeError("Checkpoint hash differs from the verified official artifact.")

torch.set_num_threads(1)
torch.use_deterministic_algorithms(True)
inference.device = "cpu"
model = SweetNet(len(lib), num_classes=1075, hidden_dim=128)
state = torch.load(checkpoint, map_location="cpu", weights_only=True)
model.load_state_dict(state, strict=True)
model = model.to("cpu").eval()
glycans = ["Gal(b1-4)GlcNAc", "Man(a1-3)[Man(a1-6)]Man", "GlcNAc", "Gal(b1-4)GlcNAc"]


def embed(values, batch_size):
    with torch.inference_mode():
        return inference.glycans_to_emb(
            values, model, libr=lib, batch_size=batch_size, rep=True
        ).to_numpy()


whole = embed(glycans, 4)
single = embed(glycans, 1)
repeated = embed(glycans, 4)
reversed_rows = embed(list(reversed(glycans)), 2)[::-1]
assert whole.shape == (4, 128)
assert np.isfinite(whole).all()
assert np.array_equal(whole[0], whole[3])
assert np.array_equal(whole, repeated)
np.testing.assert_allclose(whole, single, rtol=1e-5, atol=1e-5)
np.testing.assert_allclose(whole, reversed_rows, rtol=1e-5, atol=1e-5)
assert not np.allclose(whole[0], whole[1])
vocabulary_json = json.dumps(lib, sort_keys=True, separators=(",", ":"))
packages = ["glycowork", "torch", "torch-geometric", "glyles", "rdkit", "numpy", "pandas"]
print(json.dumps({
    "checkpoint_sha256": digest,
    "vocabulary_sha256": hashlib.sha256(vocabulary_json.encode()).hexdigest(),
    "device": "cpu", "shape": list(whole.shape),
    "readout": "bn2(lin2(...)), before activation/dropout/classifier",
    "repeated_exact": True, "duplicate_exact": True,
    "max_batch_difference": float(np.max(np.abs(whole - single))),
    "max_order_difference": float(np.max(np.abs(whole - reversed_rows))),
    "finite": True, "versions": {p: metadata.version(p) for p in packages}
}, indent=2))

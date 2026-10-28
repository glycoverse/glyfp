"""Optional, offline CPU bridge for the verified legacy SweetNet checkpoint."""
import hashlib
import importlib
import importlib.metadata as package_metadata
import json
import os
from pathlib import Path
import sys
import tempfile

CHECKPOINT_SHA256 = "4ee19fc0593db0562c16585b556050694144816754750ebfe1876ed21ba0f4be"
VOCABULARY_SHA256 = "838281860e3bbddc0de91ac635bf1ba746d79bf62891e5b86ca08eff537a974e"


class UnknownTokenError(ValueError):
    pass


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(request):
    if package_metadata.version("glycowork") != "1.9.0":
        raise ValueError("Install glycowork==1.9.0; newer SweetNet architectures are incompatible.")
    if sha256(request["checkpoint"]) != CHECKPOINT_SHA256:
        raise ValueError("Checkpoint SHA-256 does not match the supported pretrained SweetNet artifact.")
    # Do not create matplotlib caches in the user's home directory.
    with tempfile.TemporaryDirectory(prefix="glyfp-mpl-") as cache:
        os.environ.setdefault("MPLCONFIGDIR", cache)
        import torch
        import networkx as nx
        from glycowork.glycan_data.loader import lib
        from glycowork.ml.models import SweetNet
        import glycowork.ml.inference as inference
        from glycowork.motif.graph import glycan_to_nxGraph

        vocabulary = json.dumps(lib, sort_keys=True, separators=(",", ":"))
        vocabulary_hash = hashlib.sha256(vocabulary.encode()).hexdigest()
        if vocabulary_hash != VOCABULARY_SHA256:
            raise ValueError("Pretrained token vocabulary has changed.")
        packages = ["glycowork", "torch", "torch-geometric", "glyles", "rdkit",
                    "numpy", "pandas", "networkx", "scipy", "scikit-learn"]
        modules = ["glycowork.ml.models", "glycowork.ml.processing",
                   "glycowork.ml.inference", "glycowork.motif.graph",
                   "glycowork.motif.processing", "glycowork.glycan_data.loader"]
        metadata = {
            "model": "sweetnet-glycowork-1.9.0-species",
            "checkpoint_sha256": CHECKPOINT_SHA256,
            "vocabulary_sha256": vocabulary_hash,
            "preprocessing": "glyfp-sweetnet-v1:unknown-inferred-root-only;exact-internal-links;directed-graph-isomorphism",
            "readout": "bn2(lin2):before-activation-dropout-classifier",
            "dimensions": 128, "device": "cpu", "batch_size": 1,
            "versions": {name: package_metadata.version(name) for name in packages},
            "code_sha256": {
                name: sha256(importlib.import_module(name).__file__) for name in modules
            },
        }
        metadata["code_sha256"]["glyfp.bridge"] = sha256(__file__)
        metadata["code_sha256"]["glyfp.preprocessing"] = hashlib.sha256(
            request["preprocessing_source"].encode()
        ).hexdigest()
        expected = request.get("expected_metadata")
        if expected is not None and metadata != expected:
            raise ValueError("Encoder model, preprocessing, code, or dependency metadata differs; restore the fitted environment.")

        for item in request["glycans"]:
            graph = glycan_to_nxGraph(item["iupac"])
            labels = [attrs["string_labels"] for _, attrs in graph.nodes(data=True)]
            unknown = sorted(set(labels).difference(lib))
            if unknown:
                raise UnknownTokenError("Unsupported pretrained vocabulary token(s): " + ", ".join(unknown))
            original = nx.DiGraph()
            original.add_nodes_from((i, {"string_labels": label})
                                    for i, label in enumerate(item["labels"]))
            original.add_edges_from(item["edges"])
            if not nx.is_isomorphic(
                original, graph,
                node_match=nx.algorithms.isomorphism.categorical_node_match("string_labels", None)
            ):
                raise ValueError("Model parser does not preserve the glyrepr residue/linkage graph.")

        torch.set_num_threads(1)
        torch.use_deterministic_algorithms(True)
        inference.device = "cpu"
        model = SweetNet(len(lib), num_classes=1075, hidden_dim=128)
        state = torch.load(request["checkpoint"], map_location="cpu", weights_only=True)
        model.load_state_dict(state, strict=True)
        model = model.to("cpu").eval()
        embeddings = []
        if not request["fit"]:
            # Fixed single-graph inference removes caller batch composition
            # from floating-point accumulation and BatchNorm behavior.
            with torch.inference_mode():
                for item in request["glycans"]:
                    values = inference.glycans_to_emb(
                        [item["iupac"]], model, libr=lib, batch_size=1, rep=True
                    ).to_numpy()
                    if values.shape != (1, 128):
                        raise ValueError("Unexpected pretrained readout shape.")
                    embeddings.append(values[0].tolist())
        return {"ok": True, "metadata": metadata, "embeddings": embeddings}


if __name__ == "__main__":
    try:
        request = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
        response = run(request)
    except Exception as error:
        response = {"ok": False, "error": str(error),
                    "error_type": "unknown" if isinstance(error, UnknownTokenError) else "backend"}
    Path(sys.argv[2]).write_text(json.dumps(response, allow_nan=False), encoding="utf-8")
    sys.exit(0 if response["ok"] else 1)

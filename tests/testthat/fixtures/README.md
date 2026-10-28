`sweetnet-reference.csv` records 19 full 128-dimensional CPU outputs from
unmodified glycowork 1.9.0 SweetNet + glycans_to_emb(batch_size=1, rep=True),
PyTorch 2.8.0 and torch-geometric 2.6.1. The model was in eval mode with
no gradients and one PyTorch thread. No timing measurements are recorded.

Official legacy species checkpoint SHA-256:
4ee19fc0593db0562c16585b556050694144816754750ebfe1876ed21ba0f4be

Each directed residue/linkage graph was checked for labeled isomorphism against
glycowork's parser before collecting its output. Cases include single residues,
chains, branches, fucosylation, sulfation, Neu5Ac, Neu5Gc and a five-residue core.
The fixture was generated independently of the R implementation and is checked
with a maximum absolute coordinate error of 1e-5.

`sweetnet-v2-reference.csv` uses the same 19 glycans and all 128 coordinates.
It comes from the unmodified SweetNet class, tokenizer and standard-structure
parser definitions at BojarLab/glycowork commit
199efdde6bc137dbbcd900f84f927fdbbca52eae. Definitions were extracted as Python AST
nodes without rewriting their bodies. The current parser's directed labeled
graphs were checked against glyrepr before inference. The model ran on CPU in
eval mode with no gradients and one PyTorch thread.

Official v2 checkpoint SHA-256:
bc3d901e4284554c1a245b8e6b7e291a5631c8932ca9c27eb8668b6d4424a146

Pinned source SHA-256 values:
- models.py: bdda4bdbe64b59b3701409aa61ca1557bbc27b1fd91ba2bb2fbb651ceb600fb1
- motif/processing.py: 73e1821a5670066fd4c6ebb1a861ffe474866291879f152e6ff11e11360bc843
- motif/graph.py: 7e07db2e11b2f8fcdfeee658318cdc6dff073507ff8ed920b43046c18f18bf3e

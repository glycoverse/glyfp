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

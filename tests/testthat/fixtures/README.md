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

`glycangt-untrained.rds` is an **untrained** arithmetic regression fixture, not
pretrained glycan embeddings. It contains ordinary R numeric arrays for a small
2-layer, 16-dimensional, 4-head, FFN-24 model with ORF width 64, initialized once
with PyTorch seed 314159. It also contains five graph inputs, their complete
16-coordinate outputs, and independent ORF references for 1, 3, 8 and 65 nodes.

Reference classes and functions were extracted as unchanged AST definitions from
matsui-lab/GlycanGT commit 96611518c971deb89215ca163deaf9de3a59fa32. Unused imports
for Performer were omitted; that branch is disabled and was never executed.
Executed tokenizer, QR, attention, feedforward, normalization and graph readout
bodies were unchanged. PyTorch 2.8.0 ran on CPU in eval mode without gradients.
The public graph-to-triples code, using glycowork 1.9.0's parser, independently
matched all five canonical graphs; this does not verify the gated inference
distribution or every upstream parser version.

Each reference call resets the CPU generator immediately before inference to
the first seven hexadecimal digits of SHA-256(`1729\n<canonical glyrepr string>`).
Thus its Gaussian input matches the R policy while preserving the original
one-graph ORF construction. The remaining model configuration matches the public
large preset except the explicitly reduced dimensions, layer/head counts and
8-token test vocabularies. No downloaded trained weights are in this fixture.

Selected original source SHA-256 values:
- tokenizer.py: ae027f4747d00cf3a319607e88b3f8d1c7e1eb5f8b07619fe89db886359e64cd
- orf.py: 4398994b619fe383b49202082e1bf4c615d723257fadfa53eed797078e0f4bfa
- tokengt_graph_encoder.py: 18ab384e1ace0dbb4f6489ccb36d6c6f10ed79db01552d603449ac14da5617a3
- tokengt_graph_encoder_layer.py: ce54a9d46a7a93cb01136976243fd9ee8bd74359c031b92c9b5ff6a689da1a3c
- multihead_attention.py: 872c4851bc72252a4e1a25468c09d55712fd075186a30fa343d81f969d3b4604

`glycangt-reference.csv` contains independent CPU outputs from the **official
pretrained large mask-35 checkpoint**, not the small untrained fixture above.
It covers 12 multi-residue glycans and all 768 coordinates, including branched,
sialylated, fucosylated and sulfated structures and an eight-residue chain.
Official weight SHA-256:
d490d18795ee6d4746b9b4abd0a5c733d8f7e4261998d050c1b741b65e530a54

The reference uses the unchanged public GlycanGT definitions at commit
96611518c971deb89215ca163deaf9de3a59fa32, the full large configuration and all
202 official state-dict entries loaded strictly. Public graph-to-triples code
with the authors' specified glycowork 1.6.2 parser independently verifies node
and directed edge identities. The reference model receives the upstream edge
order. Evaluation uses PyTorch 2.8.0, float32, CPU, no gradients and batch size1,
with graph-derived seeds recorded alongside each canonical input. Both original
vocabulary JSON blobs match the pinned Hugging Face repository metadata exactly.
This tests numerical equivalence under the recorded ORF policy; it does not
reproduce stochastic larger minibatches or claim downstream benchmark accuracy.

## GIFFLAR molecular graphs and networks

`gifflar-reference.rds` contains 86 molecular graphs and 384-dimensional outputs
from the official MIT-licensed glycowork implementation at commit
`199efdde6bc137dbbcd900f84f927fdbbca52eae`. All 353 tensors of the official
`glycowork_gifflar_species.pt` state dictionary loaded strictly. The checkpoint
SHA-256 is `0dab7c3b61aad95dab37d4974da4b57b3fd985128dc09ad5da846e82adb9519d`,
from Hugging Face revision `52187c054b9c988ddb2a5a5af9ee64b51004fd1e`.
Only graph arrays and model outputs are bundled, not pretrained weights.

The upstream `smiles.py`, `GIFFLARTransform`, collator and model bodies generated
these graphs and outputs. Memoization was disabled for the reference export;
installed glycowork 1.9 supplied the initial directed residue graph. The pinned
upstream vocabulary was independently downloaded and compared: all 2451 token
indices exactly equal the shared TSV resource. Every supported residue was
checked alone and as both alpha/beta donor to Gal, plus branched and modified
examples. Native tests compare SMILES, node features, all six relation multisets
and directed linkage attributes, then real-model outputs at tolerance 1e-5.

`gifflar-untrained.rds` contains ordinary arrays from a small independent upstream
model: 8 features, 2 layers, 3 classification outputs, initialization seed 4287.
Its 24-coordinate pooled outputs exercise native network operations without a
pretrained checkpoint. It is deliberately UNTRAINED, not a pretrained substitute.
Development-only reference/export scripts are retained outside the package in
`glyfp-development-notes/gifflar/`; Python is never used by the R package.

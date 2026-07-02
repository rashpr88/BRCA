
**A parsimonious positive-unlabelled learning framework for breast-cancer-related gene prioritisation**, integrating ProtT5-XL protein language model embeddings, PPI network propagation (random walk with restart), and TCGA-BRCA somatic mutation frequency.

This repository accompanies the manuscript *"A parsimonious positive-unlabelled framework reveals a network topology ceiling in breast cancer gene prioritization"*.

---

## Overview

SeqNetMut evaluates which molecular data modalities genuinely improve breast-cancer gene prioritisation under a positive-unlabelled (PU) learning protocol, and whether conclusions generalise across PPI network configurations (STRING, CPDB, and their union/intersection). The pipeline runs in eight sequential stages: feature engineering, meta-learner selection, feature ablation, multi-network evaluation, holdout benchmarking, genome-wide candidate prediction, subnetwork/hub analysis, and biological validation.

---

## Repository structure

The repository is organised by pipeline stage rather than by file type. Each folder corresponds to one stage (or a closely related pair of stages) in the manuscript's Section 3 workflow.

| Folder / file | Stage(s) | Contents |
|---|---|---|
| `Preprocessing/` | 1 | Shared preprocessing: embedding normalisation (z-score → L2 unit-hypersphere), graph construction utilities, and statistical-test helpers used consistently across all downstream stages. |
| `Protein Seq_emb/` | 1 | ProtT5-XL (1024d) sequence embedding generation from FASTA input via `Rostlab/prot_t5_xl_half_uniref50-enc`, with long-sequence chunking (512 aa windows, 256 aa overlap). |
| `GearNet_Protein_Stru_emb/` | 1 | PDB/AlphaFold structure retrieval and GearNet structural embedding generation (512d) from 3D coordinates, for the structure-embedding ablation arm. |
| `GO-BP/` | 1, 8 | QuickGO annotation retrieval per UniProt ID and GO Biological Process feature preparation; also used downstream for enrichment analysis of final candidates. |
| `getAnnotations_GO-BP(onto2vec).pl` | 1 | Onto2Vec-style GO annotation pipeline: extracts GO classes, filters non-experimental evidence codes (excludes IEA/ND), and propagates annotations up the ontology hierarchy (true-path rule) to build the ancestor-augmented GO feature set. |
| `Networks_build/` | 1, 4 | Parses STRING v11 and ConsensusPathDB v35, and constructs the four PPI network configurations (STRING, CPDB, STRING∩CPDB, STRING∪CPDB) used for cross-network robustness evaluation. |
| `Meta_leaner/` | 2 | Meta-learner selection: compares candidate classifiers (e.g., logistic regression vs. tree-based alternatives) to fix the downstream model architecture before ablation. |
| `Feature_Network_ablation/` | 3, 4 | Systematic feature ablation (20 configurations, 20×10-fold repeated undersampling cross-validation) and multi-network robustness checks, identifying seq + rwr + mut as the optimal combination. |
| `Benchmark/` | 5 | Holdout benchmarking against classical topology baselines (RWR, LocalNbr, Node2Vec+XGBoost, feature-only RF) and GNN baselines (EMOGI-style GCN, GAT) on the **fixed holdout split** (`seed=2024`, 17 HP test genes) shared with the main pipeline. |
| `MP_Hubs/` | 6, 7 | Genome-wide candidate scoring, modified-precision cutoff selection against TCGA-BRCA DEGs, Leiden community detection, and hub cartography (z-score / participation coefficient classification). |
| `DEGs_validation/` | 6, 8 | Independent validation of predicted candidates against TCGA-BRCA differentially expressed genes (STAR FPKM, UCSC Xena; |log₂FC| > 0.58 threshold). |


## Data availability

Raw and intermediate data files (PPI network `.gml` files, `.h5` embedding stores, TCGA-BRCA expression/mutation tables, the 88-gene hard-positive seed list) are not included in this repository.

Data sources: STRING v11, ConsensusPathDB v35, TCGA-BRCA (UCSC Xena STAR FPKM), COSMIC, DepMap, PAM50.

---

## Requirements

- Python 3.10
- `torch`, `torch_geometric`, `scikit-learn`, `xgboost`, `networkx`, `leidenalg`, `python-igraph`, `scipy`, `pandas`, `h5py`, `transformers`, `kneed`
- Perl 5 with `List::MoreUtils` (for the ontology preprocessing scripts)
- A CUDA-capable GPU is recommended for `Sequence.ipynb` (ProtT5-XL) and `gear_net.py` (GearNet)

Install core Python dependencies:
```bash
pip install torch torch_geometric scikit-learn xgboost networkx leidenalg python-igraph scipy pandas h5py transformers kneed
```

---

## Usage notes

Recommended run order follows the numbered stages above:

```
Preprocessing/  →  Protein Seq_emb/  →  GearNet_Protein_Stru_emb/  →  GO-BP/  →  getAnnotations_GO-BP(onto2vec).pl
      ↓
Networks_build/  →  Meta_leaner/  →  Feature_Network_ablation/
      ↓
Benchmark/
      ↓
MP_Hubs/  →  DEGs_validation/  →  GO-BP/ (enrichment)
```




# SeqNetMut

**A parsimonious positive-unlabelled learning framework for breast-cancer-related gene prioritisation**, integrating ProtT5-XL protein language model embeddings, PPI network propagation (random walk with restart), and TCGA-BRCA somatic mutation frequency.

This repository accompanies the manuscript *"A parsimonious positive-unlabelled framework reveals a network topology ceiling in breast cancer gene prioritization"*.

---

## Overview

SeqNetMut asks which molecular data modalities genuinely improve breast-cancer gene prioritisation under a positive-unlabelled (PU) learning protocol with only 88 curated seed genes, and whether those conclusions hold across PPI network configurations (STRING, CPDB, and their intersection/union).

The central result is negative and deliberate: after systematic ablation, **only three features survive** — maximum ProtT5-XL sequence similarity to a training seed (`max_seq`), random-walk-with-restart proximity to the seed set (`rwr`), and TCGA-BRCA somatic mutation frequency (`mut`) — and on a locked holdout this parsimonious logistic-regression model is statistically indistinguishable from classical topology baselines. Structural embeddings, GO embeddings, two-step Markov diffusion, log-degree, CNA frequency and expression z-scores add nothing detectable. Graph neural networks (EMOGI-style GCN, GAT) do not exceed it either. Performance in this regime is bounded by PPI topology rather than by model capacity.

The pipeline runs in eight sequential stages: feature engineering, meta-learner selection, feature ablation, multi-network evaluation, holdout benchmarking, genome-wide candidate prediction, subnetwork/hub analysis, and biological validation.

### Locked configuration

| Component | Setting |
|---|---|
| Feature set | `max_seq` + `rwr` + `mut` (3 of 9 candidate features) |
| Meta-learner | Logistic regression, `class_weight="balanced"`, `C=1.0` |
| PPI network | STRING ∪ CPDB (17,442 nodes, 291,372 edges) |
| Seed genes | 88 curated breast-cancer-related genes (71 train / 17 holdout) |
| Holdout split | `seed=2024`; 17 positives + 855 background (n = 872, prevalence 0.0195) |
| Cross-validation | 20 repeats × 10-fold stratified, negative:positive ratio 10:1 |

### Headline results

| Metric | Value |
|---|---|
| Cross-validated AUPR (STRING ∪ CPDB) | 0.736 ± 0.028 |
| Holdout AUPR | 0.553 (95% bootstrap CI 0.309–0.792) |
| Holdout AUROC / F<sub>max</sub> | 0.883 / 0.625 |
| Prioritised candidates | 100 (from 17,354 non-seed genes) |
| Candidate subnetwork | 188 nodes, 3,159 edges, Q = 0.449, 3 modules |
| Degree-matched seed permutation | 99 permutations; CV AUPR 0.271 ± 0.021 under permutation |

---

## Repository structure

The repository is organised by pipeline stage rather than by file type. Each folder corresponds to one stage (or a closely related pair of stages) in the manuscript's Section 3 workflow.

| Folder / file | Stage(s) | Contents |
|---|---|---|
| `Preprocessing/` | 1 | Shared preprocessing: embedding normalisation (global z-score → L2 unit-hypersphere), graph construction utilities, and statistical-test helpers used consistently across all downstream stages. |
| `Protein Seq_emb/` | 1 | ProtT5-XL (1024d) sequence embedding generation from FASTA input via `Rostlab/prot_t5_xl_half_uniref50-enc`, with long-sequence chunking (512 aa windows, 256 aa overlap) and mean pooling across chunks. |
| `GearNet_Protein_Stru_emb/` | 1 | PDB/AlphaFold2 structure retrieval and GearNet-Edge structural embedding generation (512d) from Cα contact geometry, for the structure-embedding ablation arm. |
| `GO-BP/` | 1, 8 | QuickGO annotation retrieval per UniProt ID and GO Biological Process feature preparation; also used downstream for DAVID enrichment analysis of final candidates and modules. |
| `getAnnotations_GO-BP(onto2vec).pl` | 1 | Onto2Vec-style GO annotation pipeline (`getClasses.pl` → `getAnnotations.pl` → `AddAncestors.pl`): extracts GO classes from ontology axioms, filters non-experimental evidence codes (excludes IEA/ND), and propagates annotations up the ontology hierarchy (true-path rule) to build the ancestor-augmented GO feature set (200d). |
| `Networks_build/` | 1, 4 | Parses STRING v11 and ConsensusPathDB v35, and constructs the four PPI network configurations (STRING, CPDB, STRING ∩ CPDB, STRING ∪ CPDB) used for cross-network robustness evaluation. |
| `Meta_leaner/` | 2 | Meta-learner selection on the fixed six-feature vector: logistic regression vs. two MLP architectures vs. XGBoost, compared by De Long test on paired out-of-fold probability vectors. |
| `Feature_Network_ablation/` | 3, 4 | Systematic feature ablation (eight primary configurations, each under 20 × 10-fold repeated undersampling cross-validation) plus three supplementary diffusion-variant configurations on the winning network, and multi-network robustness checks. Identifies `seq + rwr + mut` as the locked combination. |
| `Benchmark/` | 5 | Holdout benchmarking against classical topology baselines (RWR, Majority Vote, Hishigaki, NIAPU, two-step Markov) and GNN baselines (EMOGI-style GCN, GAT) on the **fixed holdout split** (`seed=2024`, 17 HP test genes) shared with the main pipeline. `benchmark.py` is the canonical split loader — every stage imports `load_data()` from it so the holdout is defined once. |
| `MP_Hubs/` | 6, 7 | Genome-wide candidate scoring, modified-precision cutoff selection against TCGA-BRCA DEGs, Leiden community detection, and hub cartography (within-module z-score / participation coefficient classification). |
| `DEGs_validation/` | 6, 8 | Independent validation of predicted candidates against TCGA-BRCA differentially expressed genes (STAR FPKM, UCSC Xena; \|log₂FC\| > 0.58, FDR < 0.05). |
| `Permutation/` | 7 | Degree-matched seed permutation control (99 permutations) producing `permutation_results.csv`. |
| `scRNA_validation/` | 8 | Single-cell characterisation of the candidate programme across malignant, immune, and stromal compartments and across molecular subtypes (R). |

---

## Pipeline detail

### Stage 1 — Feature engineering

Nine candidate features are computed per gene. All embedding-derived features are computed **per fold** against training seeds only, so no label information crosses the fold boundary.

| Feature | Description | Source dimension |
|---|---|---|
| `max_seq` | Max cosine similarity to any training seed in ProtT5-XL space | 1024 |
| `max_struct` | Max cosine similarity to any training seed in GearNet-Edge space | 512 |
| `max_go` | Max cosine similarity to any training seed in Onto2Vec GO-BP space | 200 |
| `rwr` | Random walk with restart seeded on training seeds | — |
| `markov2` | Two-step Markov diffusion from training seeds | — |
| `log_deg` | log1p(PPI degree) | — |
| `mut` | TCGA-BRCA somatic mutation frequency | — |
| `cna` | TCGA-BRCA copy-number alteration frequency | — |
| `expr` | TCGA-BRCA mean mRNA expression z-score | — |

Embedding matrices are globally z-score standardised across all *N* proteins, then L2-normalised to the unit hypersphere. Omics features are globally z-scored only (they are scalars, not vectors). Cosine similarity uses leave-self-out masking so a seed gene never matches itself.

**Sequence embeddings.** `Rostlab/prot_t5_xl_half_uniref50-enc`, half precision, encoder only. Sequences longer than 512 residues are chunked into 512-aa windows with 256-aa overlap and mean-pooled.

**Structural embeddings.** torchdrug `models.GearNet` with `hidden_dims=[512]×6`, `num_relation=7`, `edge_input_dim=59`, `num_angle_bin=8`, `short_cut=True`, `readout="sum"`, pretrained GearNet-Edge weights. Graph construction: `AlphaCarbonNode` + `SequentialEdge(max_distance=2)` + `SpatialEdge(radius=10.0 Å, min_distance=5)` + `KNNEdge(k=10, min_distance=5)`, `edge_feature="gearnet"`. Structures from RCSB mmCIF where available, otherwise AlphaFold2 models.

**GO embeddings.** QuickGO annotations per UniProt ID restricted to Biological Process, excluding IEA and ND evidence codes, with true-path ancestor propagation before vectorisation.

**PPI networks.** STRING v11 at the project confidence threshold; ConsensusPathDB v35 restricted to interactions supported by ≥ 2 source databases, with complexes of > 10 members dropped before pairwise expansion (they would otherwise generate thousands of spurious clique edges).

| Network | Nodes | Edges |
|---|---|---|
| STRING v11 | 14,043 | 174,045 |
| CPDB v35 | 15,565 | 141,630 |
| STRING ∩ CPDB | 12,166 | 24,303 |
| STRING ∪ CPDB | 17,442 | 291,372 |

**TCGA-BRCA features.** Somatic mutation frequency, CNA frequency, and mean mRNA expression z-scores are retrieved from the `brca_tcga` study via the cBioPortal public API. Note that these are *not* the same data used for differential expression in Stage 6 — see the Data section.

### Stage 2 — Meta-learner selection

Four classifiers are compared on the fixed six-feature vector (`max_seq`, `max_struct`, `max_go`, `rwr`, `markov2`, `log_deg`) with fixed hyperparameters (no grid search):

- **LR** — `class_weight="balanced"`, `max_iter=2000`, `C=1.0`
- **MLP_16** — `(16,)`, `alpha=1e-2`, early stopping, `validation_fraction=0.15`
- **MLP_32_16** — `(32, 16)`, `alpha=1e-3`, early stopping
- **XGB** — 300 trees, `max_depth=3`, `lr=0.06`, `subsample=0.8`, `colsample_bytree=0.8`, `min_child_weight=5`, `reg_lambda=5.0`, `scale_pos_weight=neg/pos`

Comparison uses the **De Long test on paired out-of-fold predicted-probability vectors** (n ≈ 780), not Wilcoxon on 10 per-fold AUPRs — the latter cannot resolve ΔAUC < 0.05 at n = 10. Logistic regression is retained. The holdout is never touched at this stage.

### Stages 3–4 — Feature ablation and cross-network evaluation

Eight primary feature configurations are evaluated on each of the four networks under 20 × 10-fold repeated undersampling cross-validation (`BASE_SEED=42`, `N_REPEATS=20`, `N_FOLDS=10`, `NEG_RATIO=10`):

```
FULL           max_seq, max_struct, max_go, rwr, markov2, log_deg, mut, cna, expr
seq+rwr+mut ★  max_seq, rwr, mut
seq+rwr        max_seq, rwr
rwr+omics      rwr, mut, cna, expr
seq+rwr+omics  max_seq, rwr, mut, cna, expr
topo_only      rwr, markov2, log_deg
bio_only       max_seq, max_struct, max_go
omics_only     mut, cna, expr
```

Three supplementary configurations are then run on the winning network to test whether two-step Markov diffusion can replace or augment RWR:

```
markov2+seq+mut      max_seq, markov2, mut
seq+rwr+m2+mut       max_seq, rwr, markov2, mut
rwr+m2+seq+omics     max_seq, rwr, markov2, mut, cna, expr
```

Each repeat draws a fresh negative subsample (seed = `BASE_SEED + rep`), so between-repeat SD reflects sensitivity to which negatives were drawn. Configurations are compared to the reference by paired *t*-test and Wilcoxon signed-rank across the 20 repeat-level AUPRs. `seq + rwr + mut` on STRING ∪ CPDB is the locked winner (CV AUPR 0.736 ± 0.028); the union network beats STRING in 16/20 repeats (*t*-test P = 0.031, Wilcoxon P = 0.028) and CPDB and STRING ∩ CPDB in 20/20 (both P < 0.001).

RWR parameters throughout: restart probability α = 0.15, 50 power iterations on the row-normalised transition matrix.

### Stage 5 — Holdout benchmark

The holdout is fixed once in `benchmark.py::load_data()` and imported everywhere else:

```
HP  test : RandomState(2024).permutation(hp_idx)[:int(0.20 × 88)]  → 17 genes
U   test : RandomState(2025).permutation(unlabeled_idx)[:17 × 50]  → 850 genes
NEG test : RandomState(2026).permutation(neg_idx)[:20%]
```

SeqNetMut is compared against RWR, Majority Vote, Hishigaki, NIAPU and two-step Markov (classical topology methods), and against EMOGI-style GCN and GAT trained on **matched inputs and matched training pools** (`GNN_REPS=5`, `GNN_EPOCHS=400`, `lr=5e-4`, `weight_decay=1e-3`, `dropout=0.5`, early stopping on a 10% validation split with patience 50). AUPR uncertainty is quantified by nonparametric bootstrap; pairwise ΔAUPR is tested by paired bootstrap on common resample indices with two-sided empirical P values. AUROC comparisons use De Long. AUPR comparisons do **not** — De Long is valid for AUROC only.

### Stage 6 — Genome-wide prediction and modified-precision cutoff

All 17,354 non-seed genes are scored by the locked model fitted on 71 training seeds plus 710 sampled background genes. Candidate-list size is set by the modified-precision framework against a transcriptomic consistency reference:

- DEGs from UCSC Xena TCGA-BRCA **STAR FPKM**, 1,097 primary tumours vs. 113 solid-tissue normals, Welch's *t*-test with Benjamini–Hochberg FDR, retaining \|log₂FC\| > 0.58 and FDR < 0.05 → 4,291 DEGs.
- Seed genes are removed from the reference set before computing MP, so the cutoff is not circular.
- MP(N) = |top-N candidates ∩ independent DEGs| / N, evaluated for N ∈ [50, 2000].

MP peaks at N = 61 (MP = 0.475, 29/61 supported). The operating cutoff is extended to **N = 100** (MP = 0.470, 47/100 supported) to give the downstream subnetwork sufficient density for community detection. Differential expression is treated as a *consistency reference*, not independent validation, because it remains biologically entangled with the evidence used to curate the seeds.

### Stage 7 — Subnetwork, modules, hubs, and permutation control

The 100 candidates plus all 88 seeds induce a 188-node, 3,159-edge subgraph of STRING ∪ CPDB. Isolated nodes are retained for completeness but excluded from community detection.

- **Leiden** — `RBConfigurationVertexPartition`, minimum module size 6, random seed 42. Three principal modules, Q = 0.449: Module 1 (75 proteins, transcriptional regulation), Module 2 (50, cell division), Module 3 (59, DNA repair). Four isolated seed genes remain unassigned.
- **Hub cartography** — within-module degree z-score with z ≥ 1.5 defining hubs, participation coefficient PC > 0.5 separating inter-modular connector hubs from intra-modular provincial hubs. 14 hubs: 4 connector, 10 provincial; 8 are seeds, 6 are predicted candidates.
- **Degree-matched seed permutation** — the 71 training seeds are replaced by degree-matched random genes over **99 permutations**. Genes are binned into 20 log-degree quantiles and replacements drawn without replacement from the same bin, excluding all 88 curated seeds. Outputs land in `permutation_results.csv` (`perm`, `cv_aupr`, `modularity`, `n_modules`, `overlap_with_real`). CV AUPR collapses to 0.271 ± 0.021 and no permutation reaches the observed value. A mean of 32.7 ± 3.3 of the 100 candidates are still recovered under permutation (max 40), concentrated in the transcriptional-regulation module — these are generic high-degree hubs and are best read as calibration controls rather than discoveries.

### Stage 8 — Biological validation

- **GO-BP enrichment** — DAVID 2022, FDR-adjusted P < 0.05, on the full candidate set and on each module separately. Enrichment is also repeated on the **candidate-only** gene lists (seeds excluded) to confirm the programmes are not an artefact of the seed scaffold.
- **Single-cell** — an independent breast-cancer atlas is used to localise the candidate programme across malignant, immune, and stromal compartments and across molecular subtypes (analysis in R).

---

## Data availability

Raw and intermediate data files (PPI network `.gml` files, `.h5` embedding stores, TCGA-BRCA tables, the 88-gene seed list) are not included in this repository.

| Data | Source | Used for |
|---|---|---|
| STRING v11 | https://string-db.org | PPI network construction |
| ConsensusPathDB v35 | http://cpdb.molgen.mpg.de | PPI network construction |
| TCGA-BRCA mutation / CNA / expression z-scores (`brca_tcga`) | cBioPortal public API, https://www.cbioportal.org | `mut`, `cna`, `expr` model features |
| TCGA-BRCA RNA-seq STAR FPKM | UCSC Xena, https://xenabrowser.net | Differential expression for the modified-precision reference **only** |
| Protein sequences | UniProt | ProtT5-XL embeddings |
| Structures | RCSB PDB / AlphaFold DB | GearNet-Edge embeddings |
| GO annotations | EBI QuickGO | Onto2Vec GO-BP embeddings |

> The two TCGA-BRCA sources are deliberately distinct. The model's `expr` feature uses cBioPortal precomputed expression z-scores; the STAR FPKM matrix from Xena is used *only* for the differential-expression analysis behind the modified-precision cutoff.

Archived code and processed data: **Zenodo DOI — to be inserted on acceptance.**

---

## Requirements

- Python 3.10
- `torch`, `torch_geometric`, `torchdrug`, `scikit-learn`, `xgboost`, `networkx`, `leidenalg`, `python-igraph`, `scipy`, `pandas`, `statsmodels`, `h5py`, `transformers`, `biopython`, `requests`, `kneed`, `matplotlib`
- Perl 5 with `List::MoreUtils` (for the ontology preprocessing scripts)
- R with `Seurat`, `pheatmap`, `ggrepel` (single-cell stage)
- A CUDA-capable GPU is recommended for `Sequence.ipynb` (ProtT5-XL) and `gear_net.py` (GearNet-Edge)

```bash
pip install torch torch_geometric scikit-learn xgboost networkx leidenalg python-igraph \
            scipy pandas statsmodels h5py transformers biopython requests kneed matplotlib
```

`torchdrug` is required only for the GearNet-Edge structural embeddings and is best installed separately, since it pins older PyTorch builds.

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
MP_Hubs/  →  DEGs_validation/  →  Permutation/  →  GO-BP/ (enrichment)  →  scRNA_validation/
```

Notebooks were developed in Google Colab with Google Drive mounted at `/content/drive/MyDrive/MMSE BC/`. Paths at the top of each notebook must be repointed to a local directory before running elsewhere. Historic filenames in the code refer to the project's earlier working names (`MMSE-BC`, `PRIMA-BC`); these denote the same method now reported as **SeqNetMut**.

### Reproducibility

All randomness is seeded. `seed=2024` fixes the train/holdout split, `seed=42` fixes Leiden community detection and the base CV seed, and repeated-CV negative subsampling uses `42 + repeat`. Re-running any stage on the same inputs reproduces the reported numbers exactly. Because the holdout split is defined once in `benchmark.py`, no stage can silently redefine it.

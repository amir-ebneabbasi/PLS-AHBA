🧠 PLS Pipeline for Imaging Transcriptomics
A complete workflow for imaging‑transcriptomics analysis, integrating AHBA gene expression, brain phenotypes, Partial Least Squares (PLS), and pathway enrichment (FGSEA).

📦 Data Requirements
1. AHBA Gene Expression
Rows: Brain regions

Columns: Genes

2. Phenotype Data
Rows: Brain regions

Must include: phenotype column (pheno_name)

3. Spins (Spatial Null Model)
Spatial permutation indices

Each column = one permutation

4. Pathways
Gene sets in GMX (tab‑separated)

🔧 Pipeline Overview
1. AHBA Processing (abagen)
Generates regional gene expression matrices using the abagen toolbox.

Script: run_abagen.sh

Output: ROI × gene expression matrix (input to PLS)

2. Spin Permutations (Spatial Null Model)
Generates spatially constrained permutations to control for spatial autocorrelation.

Script: rotate_parcellation.R

Inputs: ROI coordinates

Methods: "hungarian", "vasa"

Output: ROI × nrot permutation index matrix

🧠 3. PLS + FGSEA Pipeline
Step 1 — Data Preprocessing
Removes genes with excessive missingness

Ensures ROI alignment across datasets

Merges phenotype + gene expression

Step 2 — PLS Regression
Uses plsdepot::plsreg1

Extracts explained variance (R²)

Aligns component signs for interpretability

Step 3 — Spin Permutation Test
Applies spatial permutations

Computes p‑values for explained variance

Step 4 — Bootstrapping
Estimates variability of gene weights

Computes corrected weights (Z‑like scores)

Step 5 — FGSEA
Uses PLS Component 1

Performs pathway enrichment analysis

📤 Outputs
1. Results_PLS_Spin.csv
Variance explained (R²)

Spin‑test p‑values

2. Results_PLS_cWeights.csv
Bootstrapped gene weights

Used for biological interpretation

3. Results_FGSEA.csv
Enriched pathways

Includes:

Normalized Enrichment Score (NES)

Leading‑edge genes

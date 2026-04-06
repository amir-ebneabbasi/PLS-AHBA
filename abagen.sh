#!/bin/bash

# Define paths
ABAGEN_BIN="/path/to/abagen"   # e.g. ~/miniconda3/envs/abagen/bin/abagen
ATLAS="/path/to/HCP-MMP1_on_MNI152_ICBM2009a_nlin.nii.gz" # glasser 
ATLAS_INFO="/path/to/modified/hcp_mmp1.0.csv" # cols should be id, label, hemisphere and structure
OUTPUT="hcp_expression.csv"
DATA_DIR="$HOME/abagen-data"

# Run abagen
$ABAGEN_BIN \
    $ATLAS \
    --atlas_info $ATLAS_INFO \
    --donors all \
    --data_dir $DATA_DIR \
    --n_proc 6 \
    --ibf_threshold 0.5 \
    --probe_selection diff_stability \
    --lr_mirror None \
    --sim_threshold None \
    --missing centroids \
    --tol 2 \
    --sample_norm srs \
    --gene_norm srs \
    --region_agg donors \
    --agg_metric mean \
    --output-file $OUTPUT

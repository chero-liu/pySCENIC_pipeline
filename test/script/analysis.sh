
#!/bin/bash
source /gpfs/oe-scrna/liuchenglong/anaconda3/bin/activate /gpfs/oe-scrna/liuchenglong/envs/pyscenic

database=/gpfs/oe-scrna/liuchenglong/RaD/pyscenic/data
input=/gpfs/oe-scrna/liuchenglong/RaD/pyscenic/pySCENIC_pipeline/test1/pyscenic_test_all_modules_True_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather/data
outdir=/gpfs/oe-scrna/liuchenglong/RaD/pyscenic/pySCENIC_pipeline/test1/pyscenic_test_all_modules_True_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather/data
species=human
num_workers=20
method=grnboost2

# Step 1: GRN inference
pyscenic grn \
    --num_workers $num_workers \
    --output $outdir/expr_mat.adjacencies.tsv \
    --method $method \
    $input/for_scenic.loom \
    $database/$species/allTFs.txt

# Step 2: Context enrichment
pyscenic ctx \
    $outdir/expr_mat.adjacencies.tsv \
    $database/$species/500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname $database/$species/motifs-v10nr_clust-nr-m0.001-o0.0.tbl \
    --expression_mtx_fname $input/for_scenic.loom \
    --mode "custom_multiprocessing" \
    --output $outdir/regulons.csv \
    --num_workers $num_workers \
    --rank_threshold 5000  \
    --auc_threshold 0.05 \
    --nes_threshold 3.0 \
    --min_orthologous_identity 0.0 \
    --max_similarity_fdr 0.001 \
    --chunk_size 100 \
    --thresholds 0.75 0.90 \
    --top_n_targets 5 10 50 \
    --min_genes 20  --all_modules

# Step 3: AUCell scoring
pyscenic aucell \
    --num_workers $num_workers \
    --output $outdir/sce_SCENIC.loom \
    $input/for_scenic.loom \
    $outdir/regulons.csv

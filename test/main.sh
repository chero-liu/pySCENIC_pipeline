#!/bin/bash
source /gpfs/oe-scrna/liuchenglong/anaconda3/bin/activate /gpfs/oe-scrna/liuchenglong/envs/pyscenic

python /gpfs/oe-scrna/liuchenglong/RaD/pyscenic/pySCENIC_pipeline/main.py \
    --input ./data_ob_v3.rds \
    --outdir ./test \
    --species human \
    --groupby group_new_celltype \
    --all_modules 'True'

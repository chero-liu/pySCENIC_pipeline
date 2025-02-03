#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pySCENIC Pipeline - Regulatory Network Analysis for Single-cell Data
---------------------------------------------------------------
Author: Chenglong Liu
Contact: chenglong.liu@oebiotech.com
Created: 2025-02
"""

import argparse
import os
import sys
from tools.utils import Datapreparation, Analysis, Visualize

def main():
    parser = argparse.ArgumentParser(description='pySCENIC pipeline script.')
    
    # Data Preparation arguments
    parser.add_argument('--input', required=True, help='Input RDS file path')
    parser.add_argument('--outdir', required=True, help='Output directory')
    parser.add_argument('--species', required=True, choices=['human', 'mouse'], help='Species for analysis')
    parser.add_argument('--subsampleid', default="no", help='Subsample ID filter')
    parser.add_argument('--subgroup', default="no", help='Subgroup filter')
    parser.add_argument('--subcluster', default="no", help='Subcluster filter')
    parser.add_argument('--subnew_celltype', default="no", help='New celltype filter')
    parser.add_argument('--predicate', default="no", help='Predicate filter')
    parser.add_argument('--shell_run', default=True, help='Run in shell mode')
    
    # Analysis arguments
    parser.add_argument('--method', default='grnboost2', choices=['grnboost2', 'genie3'], help='GRN inference method')
    parser.add_argument('--database', default='/gpfs/oe-scrna/liuchenglong/RaD/pyscenic/data', help='Database directory path')
    parser.add_argument('--tfs', default='allTFs.txt', help='TF list filename')
    parser.add_argument('--genes_vs_motifs', default='auto',
                       help='Gene-motif rankings filename')
    parser.add_argument('--motifs_tbl', default='motifs-v10nr_clust-nr-m0.001-o0.0.tbl', help='Motif table filename')
    parser.add_argument('--rank_threshold', type=int, default=5000, help='Gene rank threshold (default: 5000)')
    parser.add_argument('--auc_threshold', type=float, default=0.05, help='AUC calculation threshold (default: 0.05)')
    parser.add_argument('--nes_threshold', type=float, default=3.0, help='NES significance threshold (default: 3.0)')
    parser.add_argument('--min_orthologous_identity', type=float, default=0.0, help='Minimum orthologous identity (default: 0.0)')
    parser.add_argument('--max_similarity_fdr', type=float, default=0.001, help='Maximum similarity FDR (default: 0.001)')
    parser.add_argument('--chunk_size', type=int, default=100, help='Chunk size for module partitioning (default: 100)')
    parser.add_argument('--thresholds', type=str, default='0.75 0.90', help='Correlation coefficient thresholds (default: 0.75 0.90)')
    parser.add_argument('--top_n_targets', type=str, default='5 10 50', help='Top N target counts (default: 5 10 50)')
    parser.add_argument('--min_genes', type=int, default=20, help='Minimum number of genes (default: 20)')
    parser.add_argument('--all_modules', default='False', help='Include both positive and negative regulatory modules (default: False)')

    # Visualization arguments
    parser.add_argument('--binmethod', default='aucell', help='Binning method for visualization')
    parser.add_argument('--groupby', default='clusters', help='Grouping column for visualization')
    parser.add_argument('--threshold', type=float, default=0.0, help='AUC threshold')
    parser.add_argument('--topGenes', type=int, default=4, help='Number of top genes to show')
    parser.add_argument('--extended', default='FALSE', help='Extended visualization')
    parser.add_argument('--nclust', type=int, default=4, help='Number of clusters')
    parser.add_argument('--script_path', default=f'{os.path.dirname(__file__)}/script', help='R script path for data preparation')
    parser.add_argument('--utils_path', default=f'{os.path.dirname(__file__)}/tools/utils.r', help='R utils path for visualization')

    args = parser.parse_args()

    args.outdir = os.path.abspath(args.outdir)

    # Data Preparation
    dataprep = Datapreparation(
        input=args.input,
        outdir=f"{args.outdir}/data",
        species=args.species,
        subsampleid=args.subsampleid,
        subgroup=args.subgroup,
        subcluster=args.subcluster,
        subnew_celltype=args.subnew_celltype,
        predicate=args.predicate,
        shell_run=args.shell_run,
        script_path=f"{args.script_path}/save_seurat_counts.r"
    )
    dataprep.run()

    # Analysis
    analysis = Analysis(
        input=dataprep.outdir,
        outdir=dataprep.outdir,
        species=args.species,
        subsampleid=args.subsampleid,
        subgroup=args.subgroup,
        subcluster=args.subcluster,
        subnew_celltype=args.subnew_celltype,
        predicate=args.predicate,
        method=args.method,
        database=args.database,
        tfs=args.tfs,
        genes_vs_motifs=args.genes_vs_motifs,
        motifs_tbl=args.motifs_tbl,
        rank_threshold=args.rank_threshold,
        auc_threshold=args.auc_threshold,
        nes_threshold=args.nes_threshold,
        min_orthologous_identity=args.min_orthologous_identity,
        max_similarity_fdr=args.max_similarity_fdr,
        chunk_size=args.chunk_size,
        thresholds=args.thresholds,
        top_n_targets=args.top_n_targets,
        min_genes=args.min_genes,
        all_modules=args.all_modules,
        shell_run=args.shell_run
    )
    analysis.run()

    # Visualization
    visualize = Visualize(
        input=os.path.join(analysis.outdir, "sce_SCENIC.loom"),
        outdir=os.path.join(args.outdir, "result"),
        species=args.species,
        shell_run=args.shell_run,
        rds_filepath=args.input,
        binmethod=args.binmethod,
        groupby=args.groupby,
        threshold=args.threshold,
        topGenes=args.topGenes,
        extended=args.extended,
        nclust=args.nclust,
        subsampleid=args.subsampleid,
        subgroup=args.subgroup,
        subcluster=args.subcluster,
        subnew_celltype=args.subnew_celltype,
        predicate=args.predicate,
        script_path=f"{args.script_path}/visualize.r",
        utils_path=args.utils_path
    )
    visualize.run()

if __name__ == '__main__':
    main()

import os
import sys
import subprocess
import abc
import numpy as np
import scanpy as sc
import loompy as lp

VALID_SPECIES = {'human', 'mouse'}

class Step:
    """
    Step class
    """

    def __init__(self, input, outdir, step, species, subsampleid, subgroup, subcluster, subnew_celltype,predicate, script_path, shell_run=True):
        self.input = input
        self.outdir = outdir
        self.step = step
        if species not in VALID_SPECIES:
            raise ValueError(f"Invalid species '{species}'. Supported: {VALID_SPECIES}")
        self.species = species
        self.subsampleid = subsampleid
        self.subgroup = subgroup
        self.subcluster =subcluster
        self.subnew_celltype = subnew_celltype
        self.shell_run = shell_run
        self.script_path = script_path
        self.predicate = predicate

        if self.step == 'analysis':
            self.environment = 'source /gpfs/oe-scrna/liuchenglong/anaconda3/bin/activate /gpfs/oe-scrna/liuchenglong/envs/pyscenic'
            self.num_workers = 20
        elif self.step == 'datapreparation':
            self.environment = 'module load OESingleCell/3.0.d'
            self.num_workers = 10
        else:
            self.environment = 'module load OESingleCell/3.0.d'
            self.num_workers = 10

        if self.subsampleid == "no":
            self.subsampleid = 'all'
        if self.subgroup == "no":
            self.subgroup = 'all'
        if self.subcluster == "no":
            self.subcluster = 'all'
        if self.subnew_celltype == "no":
            self.subnew_celltype = 'all'
        if self.predicate == "no":
            self.predicate = 'all'
        
        # Ensure the output directory exists
        os.makedirs(f"{os.path.dirname(self.outdir)}/script", exist_ok=True)
        os.makedirs(f"{self.outdir}", exist_ok=True)

    @abc.abstractmethod
    def run(self):
        sys.exit("Please implement run() method.")

    @property
    def extension(self):
        return os.path.splitext(self.input)[1]

    def save_script(self, shell_script_content):
        shell_path = f"{os.path.dirname(self.outdir)}/script/{self.step}.sh"
        with open(
            shell_path,
            "w",
        ) as file:
            file.write(shell_script_content)
        
        if self.shell_run:
            subprocess.run(shell_script_content, shell=True)
    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        print(f"Init shell script for {self.step} finished.")

class Datapreparation(Step):
    def __init__(
        self,
        input: str,
        outdir: str,
        species: str,
        shell_run,
        subsampleid: str = "no",
        subgroup: str = "no",
        subcluster: str = "no",
        subnew_celltype: str = "no",
        predicate: str = 'no',
        step: str = 'datapreparation',
        script_path: str = None,
    ):
        # Initialize attributes
        super().__init__(input, outdir, step, species, subsampleid, subgroup, subcluster, subnew_celltype, predicate, script_path, shell_run)

    @property
    def shell_script(self):
        shell_script_content = f"""
#!/bin/bash
{self.environment}
Rscript {self.script_path} \\
    --input {self.input} \\
    --outdir {self.outdir} \\
    --subnew_celltype {self.subnew_celltype} \\
    --subsampleid {self.subsampleid}  \\
    --subgroup {self.subgroup}  \\
    --subcluster {self.subcluster} \\
    --predicate  {self.predicate}
"""
        return shell_script_content

    def create_loom_file(self):
        """Function to create loom file from processed counts data."""
        counts = sc.read_csv(f"{self.outdir}/counts.csv")
        row_attrs = {"Gene": np.array(counts.var_names)}
        col_attrs = {"CellID": np.array(counts.obs_names)}
        loom_file_path = f"{self.outdir}/for_scenic.loom"
        lp.create(loom_file_path, counts.X.transpose(), row_attrs, col_attrs)

    def run(self):
        if np.isin(self.extension, ['.rds', '.h5seurat']):
            self.save_script(self.shell_script)
        else:
            print(f"Unsupported file extension: {self.extension}")
            exit(1)
        if self.shell_run:
            self.create_loom_file()


class Analysis(Step):
    def __init__(
        self,
        input: str,
        outdir: str,
        species: str,
        shell_run,
        subsampleid: str = "no",
        subgroup: str = "no",
        subcluster: str = "no",
        subnew_celltype: str = "no",
        predicate: str = 'no',
        method: str = 'grnboost2',
        database: str = '/gpfs/oe-scrna/liuchenglong/RaD/pyscenic/data',
        tfs: str = 'allTFs.txt',
        genes_vs_motifs: str = 'auto',
        motifs_tbl: str = 'motifs-v10nr_clust-nr-m0.001-o0.0.tbl',
        rank_threshold: int = 5000,
        auc_threshold: float = 0.05,
        nes_threshold: float = 3.0,
        min_orthologous_identity: float = 0.0,
        max_similarity_fdr: float = 0.001,
        chunk_size: int = 100,
        thresholds: str = "0.75 0.90",
        top_n_targets: str = "5 10 50",
        min_genes: int = 20,
        all_modules: bool = True,
        step = 'analysis',
        script_path: str = None,
    ):
        super().__init__(input, outdir, step, species, subsampleid, subgroup, subcluster, subnew_celltype, predicate, script_path, shell_run)
        self.method = method
        self.database = database
        self.tfs = tfs
        self.genes_vs_motifs = genes_vs_motifs
        self.motifs_tbl = motifs_tbl

        self.rank_threshold = rank_threshold
        self.auc_threshold = auc_threshold
        self.nes_threshold = nes_threshold
        self.min_orthologous_identity = min_orthologous_identity
        self.max_similarity_fdr = max_similarity_fdr
        self.chunk_size = chunk_size
        self.thresholds = thresholds
        self.top_n_targets = top_n_targets
        self.min_genes = min_genes
        self.all_modules = all_modules
        
        if self.all_modules:
            self.min_genes = f"{self.min_genes}  --all_modules"

        if self.genes_vs_motifs == 'auto':
            if self.species == 'human':
                self.genes_vs_motifs = '500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
            else self.species == 'mouse':
                self.genes_vs_motifs = '500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather'

    @property
    def shell_script(self):
        shell_script_content = f"""
#!/bin/bash
{self.environment}

database={self.database}
input={self.input}
outdir={self.outdir}
species={self.species}
num_workers={self.num_workers}
method={self.method}

# Step 1: GRN inference
pyscenic grn \\
    --num_workers $num_workers \\
    --output $outdir/expr_mat.adjacencies.tsv \\
    --method $method \\
    $input/for_scenic.loom \\
    $database/$species/{self.tfs}

# Step 2: Context enrichment
pyscenic ctx \\
    $outdir/expr_mat.adjacencies.tsv \\
    $database/$species/{self.genes_vs_motifs} \\
    --annotations_fname $database/$species/{self.motifs_tbl} \\
    --expression_mtx_fname $input/for_scenic.loom \\
    --mode "custom_multiprocessing" \\
    --output $outdir/regulons.csv \\
    --num_workers $num_workers \\
    --rank_threshold {self.rank_threshold}  \\
    --auc_threshold {self.auc_threshold} \\
    --nes_threshold {self.nes_threshold} \\
    --min_orthologous_identity {self.min_orthologous_identity} \\
    --max_similarity_fdr {self.max_similarity_fdr} \\
    --chunk_size {self.chunk_size} \\
    --thresholds {self.thresholds} \\
    --top_n_targets {self.top_n_targets} \\
    --min_genes {self.min_genes}

# Step 3: AUCell scoring
pyscenic aucell \\
    --num_workers $num_workers \\
    --output $outdir/sce_SCENIC.loom \\
    $input/for_scenic.loom \\
    $outdir/regulons.csv
"""
        return shell_script_content

    def run(self):
        self.save_script(self.shell_script)

class Visualize(Step):
    def __init__(
        self,
        input: str,
        outdir: str,
        species: str,
        shell_run,
        rds_filepath: str,
        binmethod: str,
        groupby: str,
        threshold: int,
        topGenes: int,
        extended: str,
        nclust: int,
        utils_path: str,
        subsampleid: str = "no",
        subgroup: str = "no",
        subcluster: str = "no",
        subnew_celltype: str = "no",
        predicate: str = 'no',
        step = 'visualize',
        script_path: str = None,
    ):
        # Initialize attributes
        super().__init__(input, outdir, step, species, subsampleid, subgroup, subcluster, subnew_celltype, predicate, script_path, shell_run)
        self.rds_filepath = rds_filepath
        self.binmethod = binmethod
        self.groupby = groupby
        self.threshold = threshold
        self.topGenes = topGenes
        self.extended = extended
        self.nclust = nclust
        self.utils_path =utils_path
        
    @property
    def shell_script(self):
        shell_script_content = f"""
#!/bin/bash
{self.environment}

Rscript {self.script_path} \\
  --input {self.input} \\
  --rds_filepath {self.rds_filepath} \\
  --outdir {self.outdir} \\
  --subnew_celltype {self.subnew_celltype} \\
  --subsampleid {self.subsampleid}  \\
  --subgroup {self.subgroup}  \\
  --subcluster {self.subcluster} \\
  --predicate  {self.predicate} \\
  --binmethod {self.binmethod} \\
  --groupby {self.groupby} \\
  --threshold {self.threshold} \\
  --topGenes {self.topGenes} \\
  --extended {self.extended} \\
  --nclust  {self.nclust} \\
  --utils_path {self.utils_path}

"""
        return shell_script_content

    def run(self):
        self.save_script(self.shell_script)

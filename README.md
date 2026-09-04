# Genome-wide profiling of chromatin accessibility during differentiation and following SETDB1 depletion

**Maintainer**: Marek Gierlinski
**Collaborators**: Jane Wright, Constance Alabert

## Usage

The software in this repository is designed to be run in two stages. The first stage, including FASTQ quality control, cleaning, mapping, and peak detection is run on a computing cluster. The second stage, including all downstream analysis, is designed to be run in R on a local laptop or desktop.

### On a Linux Cluster

Copy the contents of the `cluster` folder to the cluster filesystem. Change to the `cluster` directory and create and activate a `conda` environment:

```bash
cd rna_seq
conda create --name h3me --file spec-file.txt
conda activate h3me
```

Ensure that the FASTQ files are in the `./fastq` subdirectory. They can be downloaded from [BioStudies](https://www.ebi.ac.uk/biostudies/studies/S-BSST3319).

Run Snakemake:

```bash
./run.sh
```

This will trim adapters, perform quality control, download genome files, map reads to the reference, and detect peaks.
**Note**: This shell script contains a call to *snakemake* with arguments configured for Slurm. If you're using a different cluster setup, you will need to modify the arguments accordingly.

### On a Local Laptop/Desktop

Once Snakemake has finished, we recommend using Positron or RStudio for downstream analysis. If you’re working on a different machine (e.g., I run Positron on a laptop), you will need to copy some data using:

```bash
./scripts/get_data.sh
```

**Note**: Before using this script, edit it to change the remote location, as it currently points to my directory.

Once the mapping results are available on your local machine, open an R console in Positron or RStudio and create the environment using `renv`:

```r
install.packages("renv")
renv::restore()
```

This will install all the required packages. If anything is missing, install additional packages locally using `renv::install()`.

Once the environment is set up, run the `targets` pipeline:

```r
targets::tar_make()
```

This will perform all calculations, generate data objects, figures, and tables (as [targets](https://books.ropensci.org/targets/) objects). An HTML report will be generated in the `./doc` directory.

**Note**: Creating the *DiffBind* object is a step requiring a lot of memory. You might need to run it on cluster node. Please read the instructions in the R/tar_main.R file.

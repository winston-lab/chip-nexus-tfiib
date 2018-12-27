
# TFIIB ChIP-nexus analysis pipeline

## description

A pipeline for the following TFIIB-specific ChIP-nexus analyses:

- classification of peaks into genomic categories (genic, intragenic, intergenic)
- gene ontology analysis of genes with peaks
- motif enrichment analysis near peaks

## requirements

### required software

- Unix-like operating system (tested on CentOS 7.2.1511)
- Git
- [conda](https://conda.io/docs/user-guide/install/index.html)
- [ChIP-nexus pipeline](https://github.com/winston-lab/chip-nexus)

### required files

- required only if you want to run gene ontology analyses:
    - a gene ontology mapping file in a three column, tab delimited format where the columns are common name, systematic name, and GO category:

    |      |             |            |
    | ---  | ---         | ---        |
    | ypf1 | SPAC25B8.17 | GO:1990578 |
    | nhe1 | SPAC977.10  | GO:1990578 |

- required only if you want to run motif enrichment analyses:
    - motif databases in [MEME](http://meme-suite.org/doc/meme-format.html) format

## instructions
**0**. If you haven't already done so, clone the separate ['build-annotations' pipeline](https://github.com/winston-lab/build-annotations), make a copy of the `config_template.yaml` file called `config.yaml`, and edit `config.yaml` as needed so that it points to the experimental genome FASTA file, ORF annotation BED file, transcript annotation BED file, and motif databases to be used for the ChIP-nexus pipeline. The 'build-annotations' pipeline will be used to create annotation files needed for classifying ChIP-nexus peaks into different genomic categories.

```bash

# clone the repository
git clone https://github.com/winston-lab/build-annotations.git

# move into the build-annotations pipeline directory
cd build-annotations

# make a copy of the configuration template file
cp config_template.yaml config.yaml

# edit the configuration file
vim config.yaml         # or use your favorite editor
```

**1**. Clone this repository.

```bash
git clone https://github.com/winston-lab/chip-nexus.git
```

**2**. Create and activate the `snakemake_default` virtual environment for the ChIP-nexus pipeline using conda. The virtual environment creation can take a while. If you've already created the `snakemake_default` environment from another one of my pipelines, this is the same environment, so you can skip creating the environment and just activate it.

```bash
# navigate into the pipeline directory
cd chip-nexus

# create the snakemake_default environment
conda env create -v -f envs/default.yaml

# activate the environment
source activate snakemake_default

# to deactivate the environment
# source deactivate
```

**3**. Make a copy of the configuration file template `config_template.yaml` called `config.yaml`, and edit `config.yaml` to suit your needs.

```bash
# make a copy of the configuration template file
cp config_template.yaml config.yaml

# edit the configuration file
vim config.yaml    # or use your favorite editor
```

**4**. With the `snakemake_default` environment activated, do a dry run of the pipeline to see what files will be created.

```bash
snakemake -p --use-conda --dry-run
```

**5**. If running the pipeline on a local machine, you can run the pipeline using the above command, omitting the `--dry-run` flag. You can also use N cores by specifying the `--cores N` flag. The first time the pipeline is run, conda will create separate virtual environments for some of the jobs to operate in. Running the pipeline on a local machine can take a long time, especially for many samples, so it's recommended to use an HPC cluster if possible. On the HMS O2 cluster, which uses the SLURM job scheduler, entering `sbatch slurm.sh` will submit the pipeline as a single job which spawns individual subjobs as necessary. This can be adapted to other job schedulers and clusters by modifying `slurm.sh` and `cluster.yaml`, which specifies the resource requests for each type of job.


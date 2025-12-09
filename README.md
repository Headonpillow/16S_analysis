# 16s analysis

Denoising pipeline using DADA2 algorithm to process raw .gz sequencing files from paired-end MiSeq Illumina sequencing.  

The pipeline currently takes care of trimming standard Illumina adapters from reads, filter and merge the reads (denoising), and determine the taxonomy of the different identified ASVs, as well as aligning them.  

Moreover, the pipeline results in the creation of a Phyloseq object, containing the processed samples and their metadata, for further downstream analysis in R.

The intended contents of each directory is explained in separate README.md files.

## How to run

### Option 1: Per-run isolated workflow (Recommended)

This approach keeps each analysis run isolated in its own directory.

1) Use the helper script to set up a new run:

``` {bash}
./setup_run.sh my_run
```

Or manually create the structure:

``` {bash}
mkdir -p runs/my_run/data/{raw_external,db,meta}
cp config_templates/basic.yaml runs/my_run/config.yaml
```

2) Put your raw data, database and metadata in:

* `runs/my_run/data/raw_external/` - Your .fastq.gz files
* `runs/my_run/data/db/` - SILVA database file
* `runs/my_run/data/meta/` - metadata.tsv file

> Metadata file should be in **.tsv** format, the names of the raw files should follow the convention "**{your_sample}**.R1.fastq.gz" to work.

3) Edit `runs/my_run/config.yaml` to customize settings if needed.

4) Run the pipeline from the repository root:

``` {bash}
snakemake --configfile runs/my_run/config.yaml --directory runs/my_run --use-conda --cores all all
```

This will create all outputs (`results/`, `intermediate/`, etc.) inside `runs/my_run/`, keeping your runs isolated.

### Option 2: Legacy repository-root workflow

For backwards compatibility, you can still run from the repository root:

1) Put your raw data, database and metadata in:

* data/raw_external
* data/db
* data/meta

2) Run the pipeline:

``` {bash}
snakemake --use-conda --cores all all
```

### Configuration Options

The config file supports the following options:

* `preprocess: "yes"` or `"no"` - Run preprocessing steps (FastQC, MultiQC, Trim_galore)
* `phylogeny: "yes"` or `"no"` - Build phylogenetic alignment and tree

You can also override these via command line:

``` {bash}
--config preprocess="no" phylogeny="yes"
```

> Just note that trimming doesn't happen in the pipeline when using DADA2, so if no pre-processing
takes place, the sequences will not be trimmed at all.

> the **--cores** flag specifies the amount of cores to use, you can select what you think works best.


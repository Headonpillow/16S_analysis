# 16s analysis

Denoising pipeline using DADA2 algorithm to process raw .gz sequencing files from paired-end MiSeq Illumina sequencing.  

The pipeline currently takes care of trimming standard Illumina adapters from reads, filter and merge the reads (denoising), and determine the taxonomy of the different identified ASVs, as well as aligning them.  

Moreover, the pipeline results in the creation of a Phyloseq object, containing the processed samples and their metadata, for further downstream analysis in R.

The intended contents of each directory is explained in separate README.md files.

## How to run

1) After cloning the repository, put your raw data, database and metadata in:

* data/raw_internal
* data/raw_internal/db
* data/meta

> Metadata file should be in **.tsv** format, the names of the raw files should follow the convention "**{your_sample}**.R1.fastq.gz" to work.

2) Have snakemake installed and working on your machine.

3) Run the pipeline with:

``` {bash}
snakemake --use-conda --cores all all
```

There are some parameters that you might want to use in the command line,
one of them is:

``` {bash}
--config preprocess="value"
```

Where value can be either "yes" or "no", the choice indicates if you want snakemake
to use preprocessing steps that include:

1. Generating fastqc files
2. Summarizing fastqc files with MultiQC
3. Trimming sequences with Trim_galore

Moreover an option for performing optionally the phylogenetic tree was added, as

``` {bash}
--config phylogeny="value"
```

> Just note that trimming doesn't happen in the pipeline when using DADA2, so if no pre-processing
takes place, the sequences will not be trimmed at all.

> the **--cores** flag just specify the amount of cores to use, you can select what you think works best.


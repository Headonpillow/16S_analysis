# 16s analysis

Denoising pipeline using DADA2 algorithm to process raw .gz sequencing files from paired-end MiSeq Illumina sequencing.  

The pipeline currently takes care of trimming standard Illumina adapters from reads, filter and merge the reads (denoising), and determine the taxonomy of the different identified ASVs, as well as aligning them.  

Moreover, the pipeline results in the creation of a Phyloseq object, containing the processed samples and their metadata, for further downstream analysis in R.  

The intended contents of each directory is explained in separate README.md files.

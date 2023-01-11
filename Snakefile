#IDS, = glob_wildcards("{id,[^/]+}.fastq.gz")

IDS, = glob_wildcards("data/raw_internal/{id}.fastq.gz")

# there is a problem because glob_wildcards read the whole
# directory + subdirectory, maybe to solve later implementing 
# other directories, for now excluding "/" works

rule all: 
  input:
    "results/multiQC/report_R1.html", 
    "results/multiQC/report_R2.html",
    "results/multiQC_trimmed/report_R1.html",
    "results/multiQC_trimmed/report_R2.html",
    "results/denoising/read_count_tracking.tsv",
    "results/denoising/qc.pdf",
    "results/asv/ASVs.fa", 
    "results/asv/ASVs_counts.tsv",
    "results/asv/ASVs_taxonomy.tsv",
    "results/asv/ASV_tree.nwk",
    "results/phyloseq/starting_phyla_table.tsv",
    "results/phyloseq/prevalence_graph.png",
    "results/phyloseq/Phyloseq.RData"


#################### RULES FOR QUALITY CONTROL AND TRIMMING

rule FastQC:
  conda: "16s_analysis.yml"
  input:
    expand("data/raw_internal/{id}.fastq.gz", id=IDS)
  output:
    expand("results/fastqc/{id}_fastqc.zip", id=IDS)
  shell:
    """
    [ ! -d results/fastqc ] && mkdir results/fastqc
    fastqc -o results/fastqc data/raw_internal/*.gz
    """

rule MultiQC:
  conda: "16s_analysis.yml"
  input:
    expand("results/fastqc/{id}_fastqc.zip", id=IDS)
  output:
    "results/multiQC/report_R1.html",
    "results/multiQC/report_R2.html"
  shell:
    """
    [ ! -d results/multiQC ] && mkdir results/multiQC
    multiqc -n report_R1 results/fastqc/*R1_fastqc.zip -o results/multiQC
    multiqc -n report_R2 results/fastqc/*R2_fastqc.zip -o results/multiQC
    """

rule Trim_galore:
  conda: "16s_analysis.yml"
  input:
    expand("data/raw_internal/{id}.fastq.gz", id=IDS)
  output:
    expand("intermediate/trimmed/{id}.fastq.gz", id=IDS)
  shell:
    """
    trim_galore --length 200 --paired data/raw_internal/*fastq.gz -o intermediate/trimmed
    rm intermediate/trimmed/*report.txt
    for f in intermediate/trimmed/*_val_1.fq.gz; do mv -- "$f" "${{f%_val_1.fq.gz}}.fastq.gz"; done
    for f in intermediate/trimmed/*_val_2.fq.gz; do mv -- "$f" "${{f%_val_2.fq.gz}}.fastq.gz"; done
    """

rule FastQC_trimmed:
  conda: "16s_analysis.yml"
  input: 
    expand("intermediate/trimmed/{id}.fastq.gz", id=IDS)
  output:
    expand("results/fastqc_trimmed/{id}_fastqc.zip", id=IDS)
  shell:
    """
    [ ! -d results/fastqc_trimmed ] && mkdir results/fastqc_trimmed
    fastqc -o results/fastqc_trimmed intermediate/trimmed/*.gz
    """

rule MultiQC_trimmed:
  conda: "16s_analysis.yml"
  input:
    expand("results/fastqc_trimmed/{id}_fastqc.zip", id=IDS)
  output:
    "results/multiQC_trimmed/report_R1.html",
    "results/multiQC_trimmed/report_R2.html"
  shell:
    """
    [ ! -d results/multiQC_trimmed ] && mkdir results/multiQC_trimmed
    multiqc -n report_R1 results/fastqc_trimmed/*R1_fastqc.zip -o results/multiQC_trimmed
    multiqc -n report_R2 results/fastqc_trimmed/*R2_fastqc.zip -o results/multiQC_trimmed
    """

#################### RULES FOR DENOISING AND TAX ASSIGNMENT 

# ------------ Some comments -------------
# if you wanna run the analysis parallel on two sets of data,
# that requires that you make some changes in the code
# so you don't consider only the files in the current pwd
# but maybe leave the path selection to snakemake

# OR you can implement the path inside the script as a parameter
# that you can pass to snakemake two times so to run the script
# on both paths
# ----------------------------------------

# In this rule some parameters are passed to the rule for helping with the directories
# from inside the script. These parameters can actually be easily accessed inside
# the external Python or R script.
# Remember "pwd" is always the one from where the Snakefile is running.

rule retrieve_samplenames:
  conda: "16s_analysis.yml"
  input: 
    expand("intermediate/trimmed/{id}.fastq.gz", id=IDS)
  output:
    "intermediate/trimmed/samples.txt"
  params:
    trimmed_files_loc = "intermediate/trimmed"
  script:
    "code/retrieve_names.py"

# As you can see, Snakemake also creates the necessary folders when they are 
# not present during the execution of the script, no need of a rule for them.

# In params there are some parameters that you might be finding useful to
# change in the DADA2 pipeline, the two truncLen parameters for trimming FWD
# and REV reads. They might be changed inside the script, but I prefer 
# keeping "interactive" params changeable here.
# More parameters could be added in the future.

rule denoise_reads:
  conda: "16s_analysis.yml"
  input:
    "intermediate/trimmed/samples.txt" 
  output:
    "results/denoising/read_count_tracking.tsv",
    "results/denoising/qc.pdf",
    "results/denoising/seqtab.RData"
  params:
    sample_file_loc = "intermediate/trimmed",
    results_dir = "results/denoising"
  threads: 16
  script: 
    "code/DADA2_2.0.R"

rule assign_taxonomy:
  conda: "16s_analysis.yml"
  input: 
    "intermediate/trimmed/samples.txt",
    "results/denoising/seqtab.RData"
  output:
    "results/asv/ASVs.fa", 
    "results/asv/ASVs_counts.tsv",
    "results/asv/ASVs_taxonomy.tsv"
  params:
    database = "data/raw_internal/db/SILVA_SSU_r138_2019.RData",
    dada_files_dir = "results/denoising",
    results_dir = "results/asv"
  script:
    "code/Taxonomic_assignment.R"

rule align_seqs:
  conda: "16s_analysis.yml"
  input:
    "results/asv/ASVs.fa"
  output:
    "results/asv/ASV_tree.nwk"
  shell:
    """
    mafft --auto results/asv/ASVs.fa > results/asv/ASV_tree.nwk
    """

#################### RULES FOR PHYLOSEQ ANALYSIS 

# Here the only thing saved is the phyloseq object to be reused
rule run_phyloseq:
  conda: "16s_analysis.yml"
  input:
    "intermediate/trimmed/samples.txt",
    "results/asv/ASVs_counts.tsv",
    "results/asv/ASVs_taxonomy.tsv",
    "results/asv/ASV_tree.nwk",
    "data/meta/metadata.tsv"
  output:
    "results/phyloseq/starting_phyla_table.tsv",
    "results/phyloseq/prevalence_graph.png",
    "results/phyloseq/Phyloseq.RData"
  params:
    sample_file_loc = "intermediate/trimmed",
    asv_dir = "results/asv",
    metadata_dir = "data/meta",
    results_dir = "results/phyloseq"
  script:
    "code/Phyloseq.R"
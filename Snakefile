# Default config file for backward compatibility (overridden by --configfile)
configfile: "config.yml"

# Configure directory paths from config
DATA_DIR = config.get("data_dir", "data")
RAW_DATA_SUBDIR = config.get("raw_data_subdir", "raw_external")
DB_SUBDIR = config.get("db_subdir", "db")
META_SUBDIR = config.get("meta_subdir", "meta")
INTERMEDIATE_DIR = config.get("intermediate_dir", "intermediate")
RESULTS_DIR = config.get("results_dir", "results")
LOGS_DIR = config.get("logs_dir", "logs")

# Construct full paths
RAW_DATA_PATH = f"{DATA_DIR}/{RAW_DATA_SUBDIR}"
DB_PATH = f"{DATA_DIR}/{DB_SUBDIR}"
META_PATH = f"{DATA_DIR}/{META_SUBDIR}"
TRIMMED_PATH = f"{INTERMEDIATE_DIR}/trimmed"
FILTERED_PATH = f"{INTERMEDIATE_DIR}/dada_filtered"

# Get sample IDs from raw data
IDS, = glob_wildcards(f"{RAW_DATA_PATH}/{{id}}.fastq.gz")

# there is a problem because glob_wildcards read the whole
# directory + subdirectory, maybe to solve later implementing 
# other directories, for now excluding "/" works

# Here's a conditional statement to check whether the user
# running the pipeline wants the QC and primer removal to happen
# or not.
# Based on that the input list for rule all will be different.
myoutput = [f"{RESULTS_DIR}/denoising/read_count_tracking.tsv",
f"{RESULTS_DIR}/denoising/qc.pdf",
f"{RESULTS_DIR}/asv/ASVs.fa", 
f"{RESULTS_DIR}/asv/ASVs_counts.tsv",
f"{RESULTS_DIR}/asv/ASVs_taxonomy.tsv",
f"{RESULTS_DIR}/phyloseq/starting_phyla_table.tsv",
f"{RESULTS_DIR}/phyloseq/prevalence_graph.png",
f"{RESULTS_DIR}/phyloseq/Phyloseq.RData",
f"{RESULTS_DIR}/phyloseq/ASVs_good.fasta",
f"{RESULTS_DIR}/phyloseq/plots/plot_1.tiff"]

if config['preprocess'] in ["yes"]:
  extended = [f"{RESULTS_DIR}/multiQC/report_R1.html", 
  f"{RESULTS_DIR}/multiQC/report_R2.html",
  f"{RESULTS_DIR}/multiQC_trimmed/report_R1.html",
  f"{RESULTS_DIR}/multiQC_trimmed/report_R2.html"]
  myoutput = myoutput + extended

print(myoutput)

if config['phylogeny'] in ["yes"]:
  extended = [f"{RESULTS_DIR}/phyloseq/ASV_alignment.mafft",
  f"{RESULTS_DIR}/phyloseq/ASV_alignment.mafft.treefile"]
  myoutput = myoutput + extended

rule all: 
  input:
    myoutput

#################### RULES FOR QUALITY CONTROL AND TRIMMING

# Here's a second conditional statement to check according
# to what the user chose to run or not the preprocessing

if config['preprocess'] in ["yes"]:

  rule FastQC:
    conda: "16s_analysis.yml"
    input:
      expand(f"{RAW_DATA_PATH}/{{id}}.fastq.gz", id=IDS)
    output:
      expand(f"{RESULTS_DIR}/fastqc/{{id}}_fastqc.zip", id=IDS)
    params:
      outdir = f"{RESULTS_DIR}/fastqc",
      indir = RAW_DATA_PATH
    log:
      f"{LOGS_DIR}/fastqc.log"
    shell:
      """
      [ ! -d {params.outdir} ] && mkdir -p {params.outdir}
      mkdir -p {LOGS_DIR}
      fastqc -o {params.outdir} {params.indir}/*.gz > {log} 2>&1
      """

  rule MultiQC:
    conda: "16s_analysis.yml"
    input:
      expand(f"{RESULTS_DIR}/fastqc/{{id}}_fastqc.zip", id=IDS)
    output:
      f"{RESULTS_DIR}/multiQC/report_R1.html",
      f"{RESULTS_DIR}/multiQC/report_R2.html"
    params:
      outdir = f"{RESULTS_DIR}/multiQC",
      fastqc_dir = f"{RESULTS_DIR}/fastqc"
    log:
      f"{LOGS_DIR}/multiqc.log"
    shell:
      """
      [ ! -d {params.outdir} ] && mkdir -p {params.outdir}
      mkdir -p {LOGS_DIR}
      multiqc -n report_R1 {params.fastqc_dir}/*R1_fastqc.zip -o {params.outdir} > {log} 2>&1
      multiqc -n report_R2 {params.fastqc_dir}/*R2_fastqc.zip -o {params.outdir} >> {log} 2>&1
      """

  rule Trim_galore:
    conda: "16s_analysis.yml"
    input:
      expand(f"{RAW_DATA_PATH}/{{id}}.fastq.gz", id=IDS)
    output:
      expand(f"{TRIMMED_PATH}/{{id}}.fastq.gz", id=IDS)
    params:
      indir = RAW_DATA_PATH,
      outdir = TRIMMED_PATH
    log:
      f"{LOGS_DIR}/trim_galore.log"
    shell:
      """
      mkdir -p {params.outdir}
      mkdir -p {LOGS_DIR}
      trim_galore --illumina --clip_R1 5 --clip_R2 5 --length 200 --paired {params.indir}/*fastq.gz -o {params.outdir} > {log} 2>&1
      rm {params.outdir}/*report.txt
      for f in {params.outdir}/*_val_1.fq.gz; do mv -- "$f" "${{f%_val_1.fq.gz}}.fastq.gz"; done
      for f in {params.outdir}/*_val_2.fq.gz; do mv -- "$f" "${{f%_val_2.fq.gz}}.fastq.gz"; done
      """

  rule FastQC_trimmed:
    conda: "16s_analysis.yml"
    input: 
      expand(f"{TRIMMED_PATH}/{{id}}.fastq.gz", id=IDS)
    output:
      expand(f"{RESULTS_DIR}/fastqc_trimmed/{{id}}_fastqc.zip", id=IDS)
    params:
      outdir = f"{RESULTS_DIR}/fastqc_trimmed",
      indir = TRIMMED_PATH
    log:
      f"{LOGS_DIR}/fastqc_trimmed.log"
    shell:
      """
      [ ! -d {params.outdir} ] && mkdir -p {params.outdir}
      mkdir -p {LOGS_DIR}
      fastqc -o {params.outdir} {params.indir}/*.gz > {log} 2>&1
      """

  rule MultiQC_trimmed:
    conda: "16s_analysis.yml"
    input:
      expand(f"{RESULTS_DIR}/fastqc_trimmed/{{id}}_fastqc.zip", id=IDS)
    output:
      f"{RESULTS_DIR}/multiQC_trimmed/report_R1.html",
      f"{RESULTS_DIR}/multiQC_trimmed/report_R2.html"
    params:
      outdir = f"{RESULTS_DIR}/multiQC_trimmed",
      fastqc_dir = f"{RESULTS_DIR}/fastqc_trimmed"
    log:
      f"{LOGS_DIR}/multiqc_trimmed.log"
    shell:
      """
      [ ! -d {params.outdir} ] && mkdir -p {params.outdir}
      mkdir -p {LOGS_DIR}
      multiqc -n report_R1 {params.fastqc_dir}/*R1_fastqc.zip -o {params.outdir} > {log} 2>&1
      multiqc -n report_R2 {params.fastqc_dir}/*R2_fastqc.zip -o {params.outdir} >> {log} 2>&1
      """

else:
  rule move_raw_files_to_trimmed:
    input:
      expand(f"{RAW_DATA_PATH}/{{id}}.fastq.gz", id=IDS)
    output:
      expand(f"{TRIMMED_PATH}/{{id}}.fastq.gz", id=IDS)
    params:
      indir = RAW_DATA_PATH,
      outdir = TRIMMED_PATH
    shell:
      """
      mkdir -p {params.outdir}
      cp {params.indir}/*.gz {params.outdir}
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
    expand(f"{TRIMMED_PATH}/{{id}}.fastq.gz", id=IDS)
  output:
    f"{TRIMMED_PATH}/samples.txt"
  params:
    trimmed_files_loc = TRIMMED_PATH
  log:
    f"{LOGS_DIR}/retrieve_samplenames.log"
  script:
    "code/retrieve_names.py"

# As you can see, Snakemake also creates the necessary folders when they are 
# not present during the execution of the script, no need of a rule for them.

rule denoise_reads:
  conda: "r.yml"
  input:
    samples = f"{TRIMMED_PATH}/samples.txt"
  output:
    track = f"{RESULTS_DIR}/denoising/read_count_tracking.tsv",
    qc = f"{RESULTS_DIR}/denoising/qc.pdf",
    seqtab = f"{RESULTS_DIR}/denoising/seqtab.RData"
  params:
    # output directories (script will create/use this dir)
    results_dir = f"{RESULTS_DIR}/denoising",
    intermediate_filtered_dir = FILTERED_PATH
  log:
    f"{LOGS_DIR}/denoise_reads.log"
  script: 
    "code/DADA2_2.0.R"

rule assign_taxonomy:
  conda: "r.yml"
  input: 
    samples = f"{TRIMMED_PATH}/samples.txt",
    seqtab = f"{RESULTS_DIR}/denoising/seqtab.RData"
  output:
    fa = f"{RESULTS_DIR}/asv/ASVs.fa", 
    counts = f"{RESULTS_DIR}/asv/ASVs_counts.tsv",
    tax = f"{RESULTS_DIR}/asv/ASVs_taxonomy.tsv"
  params:
    database = f"{DB_PATH}/{config.get('silva_db', 'SILVA_SSU_r138_2_2024.RData')}",
  log:
    f"{LOGS_DIR}/assign_taxonomy.log"
  script:
    "code/Taxonomic_assignment.R"

# This rule is now doing the filtering of the NA at the phyla level
# the other thing is that the script code stops at saving the ASVs
# identified as NA or not and actually these are the ones that need to
# get aligned for phylogenetic distances and go back in R in the
# downstream analysis.
# The rule outputs the fasta files, plus an .RData file to reupload in R
# with two Physeq objects with both normalized and un-normalized counts.

# TODO: correct the location where the good and bad ASVs are out, it doesn't make sense they're 
# in the Phyloseq directory.
rule filter_taxa_and_normalization:
  conda: "r.yml"
  input:
    samples = f"{TRIMMED_PATH}/samples.txt",
    counts = f"{RESULTS_DIR}/asv/ASVs_counts.tsv",
    tax = f"{RESULTS_DIR}/asv/ASVs_taxonomy.tsv",
    fa = f"{RESULTS_DIR}/asv/ASVs.fa",
    metadata = f"{META_PATH}/{config.get('metadata_file', 'metadata.tsv')}"
  output:
    starting_phyla = f"{RESULTS_DIR}/phyloseq/starting_phyla_table.tsv",
    prevalence = f"{RESULTS_DIR}/phyloseq/prevalence_graph.png",
    phyloseq = f"{RESULTS_DIR}/phyloseq/Phyloseq.RData",
    asv_good = f"{RESULTS_DIR}/phyloseq/ASVs_good.fasta"
  log:
    f"{LOGS_DIR}/filter_taxa_and_normalization.log"
  script:
    "code/Taxa_filtering.R"

if config['phylogeny'] in ["yes"]:

  rule align_seqs:
    conda: "16s_analysis.yml"
    input:
      f"{RESULTS_DIR}/phyloseq/ASVs_good.fasta"
    output:
      f"{RESULTS_DIR}/phyloseq/ASV_alignment.mafft"
    params:
      infile = f"{RESULTS_DIR}/phyloseq/ASVs_good.fasta",
      outfile = f"{RESULTS_DIR}/phyloseq/ASV_alignment.mafft"
    log:
      f"{LOGS_DIR}/align_seqs.log"
    shell:
      """
      mkdir -p {LOGS_DIR}
      mafft --auto {params.infile} > {params.outfile} 2> {log}
      """

  rule build_tree:
    conda: "16s_analysis.yml"
    input:
      f"{RESULTS_DIR}/phyloseq/ASV_alignment.mafft"
    output:
      f"{RESULTS_DIR}/phyloseq/ASV_alignment.mafft.treefile"
    params:
      infile = f"{RESULTS_DIR}/phyloseq/ASV_alignment.mafft"
    log:
      f"{LOGS_DIR}/build_tree.log"
    shell:
      """
      mkdir -p {LOGS_DIR}
      iqtree -s {params.infile} -m GTR -B 1000 -alrt 1000 -T AUTO --redo-tree > {log} 2>&1
      """

#################### RULES FOR DOWNSTREAM PHYLOGENETIC ANALYSIS

rule run_phyloseq_analysis:
  conda: "r.yml"
  input:
    phyloseq = f"{RESULTS_DIR}/phyloseq/Phyloseq.RData"
  params:
    # full path to output plots directory
    out_dir = f"{RESULTS_DIR}/phyloseq/plots"
  output:
    f"{RESULTS_DIR}/phyloseq/plots/plot_1.tiff"
  log:
    f"{LOGS_DIR}/phyloseq_analysis.log"
  script:
    "code/Phyloseq.R"

# TODO: add a feature for users to use their own ASV and TAX table?
# TODO: add error messages for stuff like sequencing not merging etc.
# TODO: add logs for some tools

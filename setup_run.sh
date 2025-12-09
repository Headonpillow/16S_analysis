#!/bin/bash
# setup_run.sh - Helper script to set up a new run directory

# repository root (script directory)
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# default shared conda-prefix (can be overridden by exporting CONDA_PREFIX before running)
CONDA_PREFIX="${CONDA_PREFIX:-${REPO_ROOT}/.snakemake/conda}"

if [ -z "$1" ]; then
    echo "Usage: $0 <run_name>"
    echo "Example: $0 my_experiment"
    exit 1
fi

RUN_NAME=$1
RUN_DIR="runs/${RUN_NAME}"

echo "Creating run directory: ${RUN_DIR}"

# Create directory structure
mkdir -p "${RUN_DIR}/data/"{raw_external,db,meta}

# Copy config template
cp config_templates/basic.yaml "${RUN_DIR}/config.yaml"

echo "Run directory created successfully!"
echo ""
echo "Next steps:"
echo "1. Place your data files in:"
echo "   - ${RUN_DIR}/data/raw_external/  (*.fastq.gz files)"
echo "   - ${RUN_DIR}/data/db/             (SILVA database)"
echo "   - ${RUN_DIR}/data/meta/           (metadata.tsv)"
echo ""
echo "2. Edit ${RUN_DIR}/config.yaml if needed"
echo ""
echo "3. Run the pipeline:"
echo "   snakemake --configfile ${RUN_DIR}/config.yaml --directory ${RUN_DIR} --use-conda --conda-prefix ${CONDA_PREFIX} --cores all all"
#!/bin/bash

# Description: Launch a single WDL workflow as a Slurm head job using Toil.
#              Unlike toil_sbatch_slurm.sh, this is for single (non-array) jobs
#              where all input is in one JSON file (e.g. ntsm_eval_workflow).

#SBATCH --job-name=toil-run
#SBATCH --cpus-per-task=4
#SBATCH --threads-per-core=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --time=7-00:00
#SBATCH --partition=long
#SBATCH --output=slurm_logs/submission_%x_%j.log


set -e

while [ $# -gt 0 ]; do
  case "$1" in
    -h|--help)
     cat << 'EOF'

      Launch a single WDL workflow as a Slurm head job using Toil.

      Usage:

      sbatch toil_sbatch_single_job.sh \
        --wdl /path/to/workflow.wdl \
        --input_json /path/to/inputs.json \
        --output_dir /path/to/output/dir \
        --output_file /path/to/outputs.json

      Options:

      -h, --help          Show brief help
      -w, --wdl           Path to the WDL file
      -i, --input_json    Path to the input JSON file
      -o, --output_dir    Directory for workflow outputs
      -f, --output_file   Path for the outputs JSON file

EOF
      exit 0
      ;;
    -w|--wdl)
      shift; WDL_PATH=$1; shift ;;
    -i|--input_json)
      shift; INPUT_JSON=$1; shift ;;
    -o|--output_dir)
      shift; OUTPUT_DIR=$1; shift ;;
    -f|--output_file)
      shift; OUTPUT_FILE=$1; shift ;;
    *)
      break ;;
  esac
done

set -x

###############################################################################
##                             Prepare For Run                               ##
###############################################################################

WDL_NAME=$(basename ${WDL_PATH%%.wdl})

export SHARED_FILESYSTEM_RUNFOLDER=$(pwd)

mkdir -p "${SHARED_FILESYSTEM_RUNFOLDER}/toil_logs"
mkdir -p "${OUTPUT_DIR}"

SINGULARITY_CACHEDIR=${SINGULARITY_CACHEDIR:-/data/tmp/$(whoami)/cache/singularity}
MINIWDL__SINGULARITY__IMAGE_CACHE=${MINIWDL__SINGULARITY__IMAGE_CACHE:-/data/tmp/$(whoami)/cache/miniwdl}
TOIL_COORDINATION_DIR=${TOIL_COORDINATION_DIR:-/data/tmp}
TOIL_SLURM_ARGS="--time=3-0:00 --partition=${SLURM_JOB_PARTITION}"

export SINGULARITY_CACHEDIR
export MINIWDL__SINGULARITY__IMAGE_CACHE
export TOIL_COORDINATION_DIR
export TOIL_SLURM_ARGS

echo "This job is running in the ${SLURM_JOB_PARTITION} partition."

set -o pipefail

###############################################################################
##                             Launch Workflow                               ##
###############################################################################

toil-wdl-runner \
    --jobStore "${WDL_NAME}_jobstore" \
    --stats \
    --clean=never \
    --batchSystem slurm \
    --batchLogsDir "${SHARED_FILESYSTEM_RUNFOLDER}/toil_logs" \
    "${WDL_PATH}" \
    "${INPUT_JSON}" \
    --outputDirectory "${OUTPUT_DIR}" \
    --outputFile "${OUTPUT_FILE}" \
    --runLocalJobsOnWorkers \
    --retryCount 1 \
    --disableProgress \
    --caching=false \
    2>&1 | tee "${WDL_NAME}_log.txt"

toil stats --outputFile "${WDL_NAME}_stats.txt" "${WDL_NAME}_jobstore"
toil clean "${WDL_NAME}_jobstore"

###############################################################################
##                                    DONE                                   ##
###############################################################################

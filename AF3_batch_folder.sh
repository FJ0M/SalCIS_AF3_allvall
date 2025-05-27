#!/bin/bash

#### Hardware resources
#SBATCH --cpus-per-task=12  # 12 CPUs for each job
#SBATCH --mem=64G           # 64G CPU memory
#SBATCH --gres=gpu:1        # 1 GPU 
#SBATCH --partition=a100    

#### Array job configuration
#SBATCH --array=0-19%8

#### Optional
#SBATCH --time=3-00:00:00      # estimated running time; adjust as needed
#SBATCH --job-name=AF3-array
#SBATCH --output=slurm-%x-%A_%a.log   # %x=job name, %A=array job ID, %a=array task ID

#SBATCH --mail-type=end           # send email when job ends
#SBATCH --mail-user=phr361@ku.dk

set -o errexit # job exit on errors

# Load required modules
module load miniconda3/23.5.2
conda activate /projects/cpr_software/apps/condaenvs/23.5.2/alphafold3_venv
module load hmmer/3.4

# GPU and XLA flags
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95

# Change to AlphaFold3 script directory
cd /projects/cpr_software/apps/software-src/alphafold3/ 

set -o errexit # job exit on errors
# Workdir and output configuration
WORKDIR="/projects/cpr_sbmm/people/wqf929/Alpha3Pulldown/SalCIS1"

mkdir -p "${WORKDIR}/output"
# JSON files directory
JSON_DIR="${WORKDIR}/JSONs"

# Create an array of JSON files
mapfile -t JSON_FILES < <(find "${JSON_DIR}" -maxdepth 1 -name "*.json" | sort)

# Calculate the total number of JSON files
TOTAL_FILES=${#JSON_FILES[@]}

# Check if the current array task is within the number of files
if [[ ${SLURM_ARRAY_TASK_ID} -le ${TOTAL_FILES} ]]; then

    # Adjust index (subtract 1 since bash arrays are 0-indexed)
    FILE_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
    
    # Selected JSON file for this array task
    SELECTED_JSON="${JSON_FILES[${FILE_INDEX}]}"

    # Extract filename without extension for output directory
    BASENAME=$(basename "${SELECTED_JSON}" .json)

    OUTPUT_DIR="${WORKDIR}/output/${BASENAME}"
    # Create output directory
    mkdir -p "${OUTPUT_DIR}"

    # Run AlphaFold3 for the selected JSON file
    python ./run_alphafold_cprome.py \
        --json_path="${SELECTED_JSON}" \
        --output_dir="${OUTPUT_DIR}" \
        --jackhmmer_n_cpu=${SLURM_CPUS_PER_TASK} \
        --nhmmer_n_cpu=${SLURM_CPUS_PER_TASK}
else
    echo "Array task ID ${SLURM_ARRAY_TASK_ID} is out of bounds for the number of JSON files (${TOTAL_FILES})"
    exit 1
fi
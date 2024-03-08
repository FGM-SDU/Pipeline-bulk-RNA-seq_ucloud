#!/bin/bash -l

#SBATCH -A mandrup    # Names the account for tracking usage, will stay as your lab name
#SBATCH --output slurm-%x.%A.%a.log   # Names an output file for what would go to STDOUT, %x %A %a represent jobname jobid jobarrayindex  
#SBATCH --mail-user victorg@bmb.sdu.dk   # Names an address to be emailed when the job finishes
#SBATCH --mail-type END,FAIL,ARRAY_TASKS  # Specifies when the job finishes (either correctly or failed)
#SBATCH --job-name this_job   # Gives the job a name, so you can find its log file & see it in the queue status, etc
#SBATCH --nodes 1         # How many nodes to be used, ie one compute nodes
#SBATCH --mem 90G        # The job can use up to __GB of ram. It is mutually exclusive with --mem-per-cpu.
#SBATCH -p CLOUD       # Names to use the serial partition
#SBATCH --cpus-per-task 16    # How many cores on that node
##SBATCH --mem-per-cpu 2500M   # he job can use up to 2.5 GB (non-integers are not allowed) of ram per cpu, i.e. 80 GB ram. NOT IN USE
#SBATCH -t 1:23:30       # Means to use up to hour:minutes:seconds of run time before it will get killed


# run in u1-standard-16 node

# Start runtime
START=$(date +%s)
echo -e "\nStarting processing"

#we set OMP_NUM_THREADS to the number of available cores
echo "Running on $SLURM_CPUS_ON_NODE CPU cores"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# If using #SBATCH --cpus-per-task then
#we set MEM_PER_THREADS to the max memory per CPU
#echo "Running on $SLURM_MEM_PER_CPU M max RAM per CPU"
#export MEM=$((${SLURM_CPUS_PER_TASK}*${SLURM_MEM_PER_CPU}))
#echo "Running on ${MEM} M max RAM"


# else using #SBATCH --mem then
echo "Running on $SLURM_MEM_PER_NODE M max RAM"

MEM=$(((${SLURM_MEM_PER_NODE}/1024)-5)) # Leave 5 GB ram aside of the available RAM
MEM_PER_THREAD=$((${MEM}/${SLURM_CPUS_PER_TASK}))

echo "${MEM} G max total mem"
echo "${MEM_PER_THREAD} G max mem per thread"

# Input Files
 
RAW_FASTQ_FILE=${1?Missing input fastq file}
#RAW_FASTQ_FILE_2=${2?Missing input fastq file}
RAW_FASTQ_FILE_2=$(basename "${RAW_FASTQ_FILE}" _R1_001.fastq.gz)_R2_001.fastq.gz

RAW_FASTQ_FILE=$(readlink -f "${RAW_FASTQ_FILE}")
RAW_FASTQ_FILE_2=$(readlink -f "${RAW_FASTQ_FILE_2}")

# ===============================================
# Create Main Output Directory
# ================================================

NAME=$(basename "${RAW_FASTQ_FILE}" _R1_001.fastq.gz)

if [ -f "${NAME}" ]; then
  >&2 echo "Error: Output location (${NAME}) already exists as a file"
  exit 1
fi

if [ -d "${NAME}" ]; then
  echo "Warning: Output location (${NAME}) is already a directory, reusing, could overwrite"
  # If you don't want to reuse, you could make it exit 1 here to kill the script if
  # the folder already exists
else
  mkdir "${NAME}"
fi

cd "${NAME}"

QC_FASTQC=QC_FASTQC && mkdir -p "${QC_FASTQC}"
QC_AdapterRemoval=QC_AdapterRemoval && mkdir -p "${QC_AdapterRemoval}"
QC_fastq_Screen=QC_fastq_Screen && mkdir -p "${QC_fastq_Screen}"
 
# ===============================================
# Fastqc
# Adapter contamination screening
# Fastq_screen
# ================================================
 
# Output Files
 
 
# Commands
module load FastQC AdapterRemoval
 
echo "FASTQC..."
START_SUBPROCESS=$(date +%s)
 
fastqc -o ${QC_FASTQC} -t "${OMP_NUM_THREADS}" "${RAW_FASTQ_FILE}"
fastqc -o ${QC_FASTQC} -t "${OMP_NUM_THREADS}" "${RAW_FASTQ_FILE_2}"

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

cd "${QC_AdapterRemoval}"

echo "AdapterRemoval..."
START_SUBPROCESS=$(date +%s)

AdapterRemoval --threads "${OMP_NUM_THREADS}" --file1 "${RAW_FASTQ_FILE}" --file2 "${RAW_FASTQ_FILE_2}" --basename "${NAME}" --gzip

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

cd ..
cd "${QC_fastq_Screen}"


module load FastQ_Screen
echo "FastQ_Screen..."
START_SUBPROCESS=$(date +%s)

fastq_screen --aligner bwa --threads "${OMP_NUM_THREADS}" --conf /work/References/fastq_screen.conf "${RAW_FASTQ_FILE}"
fastq_screen --aligner bwa --threads "${OMP_NUM_THREADS}" --conf /work/References/fastq_screen.conf "${RAW_FASTQ_FILE_2}"

module purge

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

chmod 770 ./*
chmod 770 ./*/*
chmod 770 ./*/*/*

# Finalize
END=$(date +%s)
RUNTIME=$((END-START))
H=$((RUNTIME / 3600 ))  # Calculate hours
M=$(( (RUNTIME / 60 ) % 60 ))  # Calculate minutes
S=$(( RUNTIME % 60 ))  # Calculate seconds
echo -e "\tProcessing completed. Total run time: ${H} hours, ${M} minutes, and ${S} seconds."

echo "Done!"
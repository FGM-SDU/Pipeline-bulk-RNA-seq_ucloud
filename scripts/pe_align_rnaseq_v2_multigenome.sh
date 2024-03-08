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
#SBATCH -t 20:30:00       # Means to use up to hours:minutes:seconds of run time before it will get killed


# run in u1-standard-16 node

# Start runtime
START=$(date +%s)
echo -e "\nStarting processing"

## CPU and MEMORY settings
#we set OMP_NUM_THREADS to the number of available cores
echo "Running on $SLURM_CPUS_ON_NODE CPU cores"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
SAM_THREADS=8 # Use only 8 threads for samtools sort 

# If using #SBATCH --cpus-per-task then
#we set MEM_PER_THREADS to the max memory per CPU
#echo "Running on $SLURM_MEM_PER_CPU M max RAM per CPU"
#export MEM=$((${SLURM_CPUS_PER_TASK}*${SLURM_MEM_PER_CPU}))
#echo "Running on ${MEM} M max RAM"

# else using #SBATCH --mem then
echo "Running on $SLURM_MEM_PER_NODE M max RAM"

# One should only specify 80-90% of available memory 
MEM=$((${SLURM_MEM_PER_NODE}/1024)) # Memory in GB
MEM=$((9*${MEM}/10)) # Leave 19 GB ram aside of the available RAM
MEM_PER_THREAD=$((${MEM}/${SLURM_CPUS_PER_TASK}))
MEM_SAMTOOLS=$((${MEM}/${SAM_THREADS})) # Memory per thread for samtools sort

echo "${OMP_NUM_THREADS} CPUs per task"
echo "${MEM} G total mem"
echo "${MEM_PER_THREAD} G mem per thread"
echo "Samtools sorting conditions: ${SAM_THREADS} threads and ${MEM_SAMTOOLS} G Mem Per Thread"

# Read optional arguments
INCLUDE_SALMON=false  # Default value for the -t flag

## Read optional arguments    
function usage {
  echo -e "\n Usage:$(basename $0) -g <genome> <input_files>"
  echo "Options:"
  echo " -g <genome>   - Specify the genome (mm10, mm39, hg38)"
  echo "Optional argument:"
  echo " -t            - Additionally run Salmon for isoform level DGE analysis  "
  echo " -h            - Display this help message"
  exit 1
}

while getopts g:th opt; do
    case "${opt}" in
      g) GENOME="${OPTARG}"
      ;;
      t)
        INCLUDE_SALMON=true
      ;;
      h)          
        usage     
        exit 0              
      ;;
    	\?) ## Invalid option
      echo "Invalid option: -${OPTARG}" >&2
      usage
      exit 1
      ;;
      :) ## Argument required
      echo "Option -${OPTARG} requires an argument." >&2
      usage
      exit 1
      ;;
    esac
done
shift "$((OPTIND-1))"  #This tells getopts to move on to the next argument.

if [ -z "${GENOME}" ]; then
  echo "Genome option (-g) is required."
  usage
  exit 1
fi

# PATH to the reference genomes
if [ "${GENOME}" == "mm10" ]; then
  REFERENCE="/work/References/Mouse/mm10/rna_star2.7.9a/"
  TRANSCRIPTOME="/work/References/GENCODE/Gencode_mouse/release_M25/salmon1.9.0_k21"

elif [ "${GENOME}" == "mm39" ]; then
  REFERENCE="/work/References/Mouse/mm39/star2.7.9a_rnaseq/"
  TRANSCRIPTOME="/work/References/GENCODE/Gencode_mouse/release_M33/salmon1.9.0_k21"

elif [ "${GENOME}" == "hg38" ]; then
  REFERENCE="/work/References/Human/hg38_analysisSet/star2.7.9a_rnaseq/"
  TRANSCRIPTOME="/work/References/GENCODE/Gencode_human/release_44/salmon1.9.0_k21"

else
  echo "Invalid genome option: ${GENOME}"
  usage
  exit 1
fi

echo "Using reference: ${REFERENCE}"

## Input files
if [ $# -eq 0 ]; then
  echo "No input files provided."
  usage
fi

INPUT1=${1?Missing input R1_001.fastq.gz file}
INPUT2=$(basename "${INPUT1}" _R1_001.fastq.gz)_R2_001.fastq.gz
#INPUT2=${2?Missing input R2_001.fastq.gz file}

echo "Input file: ${INPUT1} ${INPUT2}"

PREFIX=$(basename "${INPUT1}" _R1_001.fastq.gz)

## Commands
module load STAR picard SAMtools

echo "Aligning..."
START_SUBPROCESS=$(date +%s)
STAR --runMode alignReads --runThreadN "${OMP_NUM_THREADS}" --genomeDir ${REFERENCE} --genomeLoad LoadAndKeep --readFilesIn "${INPUT1}" "${INPUT2}" \
--readFilesCommand zcat --outFileNamePrefix ${PREFIX} --outSAMunmapped Within --outFilterType BySJout \
--outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 10 --outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 --alignSJDBoverhangMin 1
END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "SAM to BAM and sorting..."
START_SUBPROCESS=$(date +%s)
samtools view -@ "${OMP_NUM_THREADS}" -Shu "${PREFIX}Aligned.out.sam" | samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" -O bam -o "${PREFIX}.tmp.bam"

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Picard AddOrReplaceReadGroups ..."
START_SUBPROCESS=$(date +%s)

ID=$(zcat "${INPUT1}" | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4) ## read group identifier 
LB=$(echo "${INPUT1}" | cut -d "_" -f1,2)                                ## library ID
PL="illumina"                                                           ## platform (e.g. illumina, solid)
PU=$(echo "${ID}"."${LB}")                                                            ##Platform Unit
SM=$(echo "${INPUT1}" | cut -d"_" -f1)                                          ##sample ID

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups -I "${PREFIX}.tmp.bam" -O "${PREFIX}.bam" --RGID "${ID}" --RGLB "${LB}" --RGPL "${PL}" --RGPU "${PU}" --RGSM "${SM}" --VALIDATION_STRINGENCY LENIENT --QUIET true --VERBOSITY ERROR

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."


echo "Indexing..."
samtools index -b "${PREFIX}.bam" "${PREFIX}.bai"


rm "${PREFIX}Aligned.out.sam" "${PREFIX}Log.progress.out" "${PREFIX}.tmp.bam"

LOGS=Align_logs && mkdir -p "${LOGS}"

mv *Log* ./Align_logs
mv *SJ.out* ./Align_logs

module unload STAR SAMtools picard

# Salmon alignment
if [ "${INCLUDE_SALMON}" = true ]; then
  # ===============================================
  # Create Main Output Directory
  # ================================================
  if [ -f "${PREFIX}" ]; then
    >&2 echo "Error: Output location (${PREFIX}) already exists as a file"
    exit 1
  fi

  if [ -d "${PREFIX}" ]; then
    echo "Warning: Output location (${PREFIX}) is already a directory, reusing, could overwrite"
    # If you don't want to reuse, you could make it exit 1 here to kill the script if
    # the folder already exists
  else
    mkdir "${PREFIX}"
  fi

  echo "Salmon..."
  START_SUBPROCESS=$(date +%s)
  # --seqBias: learn and correct for sequence-specific biases in the input data
  # --useVBOpt: use variational Bayesian EM algorithm.
  # --numBootstraps: compute bootstrapped abundance estimates. Required for isoform level DGE analysis for estimation of technical variance.

  module load Salmon

  salmon quant -i ${TRANSCRIPTOME} -p "${OMP_NUM_THREADS}" -l A \
  -1 "${INPUT1}" -2 "${INPUT2}" -o "${PREFIX}/${PREFIX}_salmon" \
  --seqBias --useVBOpt --numBootstraps 30

  END_SUBPROCESS=$(date +%s)
  RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
  H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
  M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
  S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
  echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."
fi

# Finalize
END=$(date +%s)
RUNTIME=$((END-START))
H=$((RUNTIME / 3600 ))  # Calculate hours
M=$(( (RUNTIME / 60 ) % 60 ))  # Calculate minutes
S=$(( RUNTIME % 60 ))  # Calculate seconds
echo -e "\tProcessing completed. Total run time: ${H} hours, ${M} minutes, and ${S} seconds."

echo "Done!"
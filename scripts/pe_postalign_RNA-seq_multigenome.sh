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
#SBATCH -t 1-20:6:30       # Means to use up to days-hours:minutes:seconds of run time before it will get killed


# run in u1-standard-16 node

# Start runtime
START=$(date +%s)
echo -e "\nStarting processing"

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

## Read optional arguments    
function usage {
  echo -e "\n Usage:$(basename $0) -g <genome>  <input_files>"
  echo "Options:"
  echo " -g <genome>   - Specify the genome (mm10, mm39, hg38)"
  echo " -h            - Display this help message"
  exit 1
}

while getopts g:h opt; do
    case "${opt}" in
      g) GENOME="${OPTARG}"
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
  REFERENCE=/work/References/Mouse/mm10/mm10.fa

  BLACKLIST=/work/References/Blacklist/Blacklist_Mouse/mm10.blacklist.v2_merge1000.bed
  GENE=/work/References/GENCODE/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf
  GENE_RNASEQC=/work/References/GENCODE/Gencode_mouse/release_M25/gencode.vM25.annotation.genes.gtf
  REF_FLAT=/work/References/GENCODE/Gencode_mouse/release_M25/gencode.vM25.annotation.refFlat.txt
  RIBOSOMAL_INTERVALS=/work/References/GENCODE/Gencode_mouse/release_M25/gencode.vM25.annotation.rRNA.bed.interval_list
  REF_GENE_MODEL=/work/References/RSeQC/mm10_GENCODE_vm25.bed
  REF_HOUSEKEEPING_GENES=/work/References/RSeQC/mm10.HouseKeepingGenes.bed

elif [ "${GENOME}" == "mm39" ]; then
  BLACKLIST=/work/References/Blacklist/Blacklist_Mouse/mm39.blacklist.v2_merge1000.bed
  REFERENCE="/work/References/Mouse/mm39/mm39.fa"
  GENE=/work/References/GENCODE/Gencode_mouse/release_M33/gencode.vM33.annotation.gtf
  GENE_RNASEQC=/work/References/GENCODE/Gencode_mouse/release_M33/gencode.vM33.annotation.genes.gtf
  REF_FLAT=/work/References/GENCODE/Gencode_mouse/release_M33/gencode.vM33.annotation.refFlat.txt
  RIBOSOMAL_INTERVALS=/work/References/GENCODE/Gencode_mouse/release_M33/gencode.vM33.annotation.rRNA.bed.interval_list
  REF_GENE_MODEL=/work/References/RSeQC/mm39_GENCODE_VM27.bed
  REF_HOUSEKEEPING_GENES=/work/References/RSeQC/mm39.HouseKeepingGenes.bed

elif [ "${GENOME}" == "hg38" ]; then
  REFERENCE=/work/References/Human/hg38_analysisSet/hg38.analysisSet.fa
  BLACKLIST=/work/References/Blacklist/Blacklist_Human/hg38.blacklist.v2_merge1000.bed
  GENE=/work/References/GENCODE/Gencode_human/release_44/gencode.v44.annotation.gtf
  GENE_RNASEQC=/work/References/GENCODE/Gencode_human/release_44/gencode.v44.annotation.genes.gtf
  REF_FLAT=/work/References/GENCODE/Gencode_human/release_44/gencode.v44.annotation.refFlat.txt
  RIBOSOMAL_INTERVALS=/work/References/GENCODE/Gencode_human/release_44/gencode.v44.annotation.rRNA.bed.interval_list
  REF_GENE_MODEL=/work/References/RSeQC/hg38_GENCODE_V42_Basic.bed
  REF_HOUSEKEEPING_GENES=/work/References/RSeQC/hg38.HouseKeepingGenes.bed

else
  echo "Invalid genome option: ${GENOME}"
  usage
  exit 1
fi

echo "Using reference: ${REFERENCE}"

# Input Files
if [ $# -eq 0 ]; then
  echo "No input files provided."
  usage
fi

RAW_BAM_FILE=${1?Missing input bam file}
RAW_BAM_FILE=$(readlink -f "${RAW_BAM_FILE}")

echo "Input file: ${RAW_BAM_FILE}"

# ===============================================
# Create Main Output Directory
# ================================================

NAME=$(basename "${RAW_BAM_FILE}" .bam)

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

# Aditional input files
HEADER_PBC=/work/References/Multiqc/pbc_libcomp_header.txt

# Output Files

INITIAL_BAM_FILE_SAMSTATS="$(basename "${RAW_BAM_FILE}" bam)samstat.qc.txt" # QC file for multiQC
INITIAL_BAM_FILE_IDXSTATS="$(basename "${RAW_BAM_FILE}" bam)idxstat.qc.txt" # QC file for multiQC
INITIAL_BAM_FILE_MAPSTATS="$(basename "${RAW_BAM_FILE}" bam)flagstat.qc.txt" # QC file

ALIGMENTSUMMARY_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)alignmentSummary.picardMetrics.qc.txt" # QC file
QSCOREDIST_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)QScoreDist.picardMetrics.qc.txt" # QC file
CHART_QSCOREDIST_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)QScoreDist.picardMetrics.qc.pdf" # QC graph
QBYCYCLE_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)QByCycle.picardMetrics.qc.txt" # QC file
CHART_QBYCYCLE_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)QByCycle.picardMetrics.qc.pdf" # QC graph
RNASEQ_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)CollectRnaSeqMetrics.picardMetrics.qc.txt"	# QC file
CHART_RNASEQ_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)CollectRnaSeqMetrics.picardMetrics.qc.pdf"	# QC graph


# Only for Pair-End
INSERTSIZE_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)InsertSize.picardMetrics.qc.txt" # QC file
CHART_INSERTSIZE_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)InsertSize.picardMetrics.qc.pdf" # QC graph
#LIBRARYCOMPLEXITY_PICARDMETRICS="$(basename "${RAW_BAM_FILE}" bam)libraryComplexity.picardMetrics.qc.txt" # QC file

# Library complexity
PBC_FILE_QC="$(basename "${RAW_BAM_FILE}" bam)qc.pbc_mqc.tsv" # QC file
# Library complexity description
# TotalReadPairs(MT)  DistinctReadPairs(M1) OneReadPair TwoReadPairs(M2)  NRF=Distinct/Total  PBC1=OnePair/Distinct PBC2=OnePair/TwoPair

LIBRARYCOMPLEXITY_PRESEQ="$(basename "${RAW_BAM_FILE}" bam)libraryComplexity.preseq.qc.tab" # QC file

# Output Files of Deeptools
FINGERPRINTS="$(basename "${RAW_BAM_FILE}" bam)qc.fingerprints.txt" # QC file
CHART_FINGERPRINTS="$(basename "${RAW_BAM_FILE}" bam)qc.fingerprints.png" # QC graph
CHART_COVERAGE="$(basename "${RAW_BAM_FILE}" bam)qc.coveragePlot.png" # QC graph
ESTIMATE_READ_FILTERING="$(basename "${RAW_BAM_FILE}" bam)qc.estimatereadfiltering.txt" # QC file
BIGWIG_COVERAGE="$(basename "${RAW_BAM_FILE}" bam)dupmark.bs10.bw" # Use this file with IGV browser to inspect the signal
MATRIX_SCALED_GENES="$(basename "${RAW_BAM_FILE}" bam)qc.matrix.scaled.gz"
PLOTPROFILE_DATA="$(basename "${RAW_BAM_FILE}" bam)qc.plotprofile.txt" # QC file
PLOTPROFILE_PLOT="$(basename "${RAW_BAM_FILE}" bam)qc.plotprofile.png" # QC file

# Output Files for Markduplicates
MARKDUP_BAM_FILE="$(basename "${RAW_BAM_FILE}" bam)raw.dupmark.bam"	# Intermediate file later removed
NMSRT_MARKDUP_QFILT_BAM_FILE="$(basename "${RAW_BAM_FILE}" bam)qfilt.dupmark.nmsrt.bam" # Intermediate file later removed	
FINAL_MARKDUP_BAM_FILE="$(basename "${RAW_BAM_FILE}" bam)dupmark.bam"	# Bam file with duplicates Flagged for allow a downstream program to remove/ignore duplicates or keep them
FINAL_MARKDUP_INDEX_FILE="$(basename "${RAW_BAM_FILE}" bam)dupmark.bai"	# Index file
DUP_FILE_QC="$(basename "${RAW_BAM_FILE}" bam)filt.dup.qc.txt"	# QC file
FINAL_NMSRT_MARKDUP_BAM_FILE="$(basename "${RAW_BAM_FILE}" bam)dupmark.nmsrt.bam" 	# Name sorted bam file used by FeatureCounts

# Temporary (Tmp) files: transitory files in the proccesses of fixing mates in pair-end data and filtering contigs (remove of chrMT, chrUn, etc.)
TMP_QFILT_BAM_FILE="$(basename "${RAW_BAM_FILE}" bam)qfilt.bam"
TMP_QFILT_BED_FILE="$(basename "${RAW_BAM_FILE}" bam)qfilt.bed"
TMP_NMSRT_BAM_FILE="$(basename "${RAW_BAM_FILE}" bam)qfilt.nmsrt.bam"
SAM_FILE_TMP="$(basename "${RAW_BAM_FILE}" bam)filt_contigs.sam"
BAM_FILE_CHR_FILETERED_TMP="$(basename "${RAW_BAM_FILE}" bam)filt_contigs.bam"

# Output Files of Removing Duplicates
NODUP_BAM_FILE="$(basename "${RAW_BAM_FILE}" bam)nodup.bam"	# Coordinate sorted bam file
FINAL_NMSRT_BAM_FILE="$(basename "${RAW_BAM_FILE}" bam)nodup.nmsrt.bam"	# Name sorted bam file
FINAL_BAM_INDEX_FILE="$(basename "${RAW_BAM_FILE}" bam)nodup.bai"	# Index File

# Output Files of featureCounts
COUNTS_FULL_TABLE="$(basename "${RAW_BAM_FILE}" bam)featurecounts.txt" # Table of multiples columns containing geneIDs, several annotated features and raw counts
FINAL_COUNTS="$(basename "${RAW_BAM_FILE}" bam)RawCounts.txt"	# Two columns table: GeneID/ RawCounts


# ===============================================
# Flagstat
# Remove low MAPQ reads
# Compute Library Complexity
# ================================================
#Sort by name
#convert to bedPE and obtain fragment coordinates
#sort by position and strand
#Obtain unique count statistics

# Commands

module load SAMtools R picard BEDTools

echo "Samtools QC..."
START_SUBPROCESS=$(date +%s)

samtools stats ${RAW_BAM_FILE} > ${INITIAL_BAM_FILE_SAMSTATS}
samtools idxstats ${RAW_BAM_FILE} > ${INITIAL_BAM_FILE_IDXSTATS}
samtools flagstat ${RAW_BAM_FILE} > ${INITIAL_BAM_FILE_MAPSTATS}

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Picard QC..."
START_SUBPROCESS=$(date +%s)

java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics -R ${REFERENCE} -I ${RAW_BAM_FILE} -O ${ALIGMENTSUMMARY_PICARDMETRICS} --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --QUIET true --VERBOSITY ERROR 2> "${NAME}.picard.alignmentSummary.error.log"
java -jar $EBROOTPICARD/picard.jar QualityScoreDistribution -CHART ${CHART_QSCOREDIST_PICARDMETRICS} -I ${RAW_BAM_FILE} -O ${QSCOREDIST_PICARDMETRICS} --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --QUIET true --VERBOSITY ERROR 2> "${NAME}.picard.QScoreDist.error.log"
java -jar $EBROOTPICARD/picard.jar MeanQualityByCycle -I ${RAW_BAM_FILE} -O ${QBYCYCLE_PICARDMETRICS} -CHART ${CHART_QBYCYCLE_PICARDMETRICS} --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --QUIET true --VERBOSITY ERROR 2> "${NAME}.picard.QByCycle.error.log"

java -jar $EBROOTPICARD/picard.jar CollectRnaSeqMetrics -I ${RAW_BAM_FILE} -O ${RNASEQ_PICARDMETRICS} -STRAND NONE --REF_FLAT ${REF_FLAT} --RIBOSOMAL_INTERVALS ${RIBOSOMAL_INTERVALS} -CHART ${CHART_RNASEQ_PICARDMETRICS} --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --QUIET true --VERBOSITY ERROR 2> "${NAME}.picard.RnaSeqMetrics.error.log"

#Only Pair-End
java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics -H ${CHART_INSERTSIZE_PICARDMETRICS} -M 0.5 -I ${RAW_BAM_FILE} -O ${INSERTSIZE_PICARDMETRICS} --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --QUIET true --VERBOSITY ERROR 2> "${NAME}.picard.InsertSizeMetrics.error.log"

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

# Filtering by MAPQ
START_SUBPROCESS=$(date +%s)
MAPQ_THRESH=30
echo "Quality filtering q "${MAPQ_THRESH}" and indexing..."

samtools view -q ${MAPQ_THRESH} -Sbh ${RAW_BAM_FILE} -@ "${OMP_NUM_THREADS}" | samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" - -o ${TMP_QFILT_BAM_FILE}

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Library complexity.."
START_SUBPROCESS=$(date +%s)
# Get PBC bottlenecking metrics
  # mt: number of reads (TotalReadPairs)
  # m0: number of all genomic locations where reads mapped (DistinctReadPairs)
  # m1: number of genomic locations where only one read maps uniquely (OneReadPair)
  # m2: number of genomic locations where 2 reads map uniquely (TwoReadPairs)
  # NRF: Non-Redundant Fraction
  # PBC1: PCR Bottlenecking Coefficient 1
  # PBC2: PCR Bottlenecking Coefficient 2  

#echo -e "TotalReadPairs \tDistinctReadPairs \tOneReadPair \tTwoReadPairs \tNRF=Distinct/Total \tPBC1=OnePair/Distinct \tPBC2=OnePair/TwoPair" > ${PBC_FILE_QC}
samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" -n ${TMP_QFILT_BAM_FILE} -o ${TMP_NMSRT_BAM_FILE} # Will produce name sorted BAM

# When using bamtobed -bedpe option, it is required that the BAM file is sorted/grouped by the read name
# Convert bam to bedpe, exclude chrM, sort and count number of occurrences of each unique entry. Finally awk calculte the metrics.
bedtools bamtobed -i ${TMP_NMSRT_BAM_FILE} -bedpe | awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' | grep -v 'chrM' | \
#sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}'| cat ${HEADER_PBC} - > "${PBC_FILE_QC}" 
sort | uniq -c | awk -v NAME="${NAME}" 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%s\t%d\t%d\t%d\t%d\t%f\t%f\t%f\n",NAME,mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}'| cat ${HEADER_PBC} - > "${PBC_FILE_QC}" 

#Only Pair-End
#java -jar $EBROOTPICARD/picard.jar EstimateLibraryComplexity -I ${TMP_QFILT_BAM_FILE} -O ${LIBRARYCOMPLEXITY_PICARDMETRICS} --VALIDATION_STRINGENCY LENIENT --QUIET true --VERBOSITY ERROR 2> "${NAME}.picard.EstimateLibraryComplexity.error.log"

# Create a bed file for Preseq run more stable than from bam
(bedtools bamtobed -bedpe -i ${TMP_NMSRT_BAM_FILE} 2>/dev/null) | sort -k 1,1 -k 2,2n -k 3,3n -k 6,6 - > "${TMP_QFILT_BED_FILE}"

# Preseq
module load preseq 

preseq lc_extrap -pe -o ${LIBRARYCOMPLEXITY_PRESEQ} ${TMP_QFILT_BED_FILE} # computation in Bam files is less stable

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

module unload SAMtools R picard BEDTools preseq

# ================================
# RNA-seq QC 
# ================================
module load Qualimap

echo "qualimap..."
START_SUBPROCESS=$(date +%s)

qualimap rnaseq -bam ${RAW_BAM_FILE} -gtf ${GENE} -outdir ./${NAME}.qc.qualimap -pe --java-mem-size="${MEM}G"
qualimap bamqc -bam ${RAW_BAM_FILE} -c -outdir ./${NAME}.qc.qualimap.bamqc -nt "${OMP_NUM_THREADS}" --java-mem-size="${MEM}G"

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

module unload Qualimap
module load RNA-SeQC

echo "RNA-SeQC..."
START_SUBPROCESS=$(date +%s)

rnaseqc ${GENE_RNASEQC} ${RAW_BAM_FILE} ./${NAME}.qc.rnaseqc

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

module unload RNA-SeQC
module load SAMtools RSeQC

echo "RSeQC..."
START_SUBPROCESS=$(date +%s)

read_duplication.py -i ${RAW_BAM_FILE} -o "${NAME}.qc.rseqc"   # output.dup.pos.DupRate.xls, output.dup.seq.DupRate.xls, output.DupRate_plot.r, output.DupRate_plot.pdf
read_GC.py -i ${RAW_BAM_FILE} -o "${NAME}.qc.rseqc"  # output.GC.xls, output.GC_plot.r, output.GC_plot.pdf
inner_distance.py -i ${RAW_BAM_FILE} -o "${NAME}.qc.rseqc" -r ${REF_GENE_MODEL}
junction_annotation.py -i ${RAW_BAM_FILE} -o "${NAME}.qc.rseqc" -r ${REF_GENE_MODEL}
junction_saturation.py -i ${RAW_BAM_FILE} -o "${NAME}.qc.rseqc" -r ${REF_GENE_MODEL} # *.junctionSaturation_plot.r"

bam_stat.py -i ${RAW_BAM_FILE} > "${NAME}.qc.rseqc.bam_stat.txt"
infer_experiment.py -i ${RAW_BAM_FILE} -r ${REF_GENE_MODEL} > "${NAME}.qc.rseqc.infer_experiment.txt"
read_distribution.py -i ${RAW_BAM_FILE} -r ${REF_GENE_MODEL} > "${NAME}.qc.rseqc.read_distribution.txt"

samtools index ${RAW_BAM_FILE} "${RAW_BAM_FILE}.bai"
geneBody_coverage.py -i ${RAW_BAM_FILE} -o "${NAME}.qc.rseqc" -r ${REF_HOUSEKEEPING_GENES} # *qc.rseqc.geneBodyCoverage.txt"
tin.py -i ${RAW_BAM_FILE} -r ${REF_GENE_MODEL} 2> "${NAME}.qc.rseqc.tin.summary.txt" # tin rna-quality score

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

module unload RSeQC

# ===================================
# Mark duplicates
# Insert Size
# ====================================

module load SAMtools R picard

START_SUBPROCESS=$(date +%s)

echo "Fixing Coodinates and sorting"
samtools view -@ "${OMP_NUM_THREADS}" -Sh ${TMP_QFILT_BAM_FILE} | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' > ${SAM_FILE_TMP}
samtools view -@ "${OMP_NUM_THREADS}" -Shb ${SAM_FILE_TMP} | samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" - -o ${BAM_FILE_CHR_FILETERED_TMP}

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Marking duplicates and remove unmapped, mate unmapped not primary aligment and reads failing platform/vendor quiality checks..."
START_SUBPROCESS=$(date +%s)

java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I ${BAM_FILE_CHR_FILETERED_TMP} -O ${MARKDUP_BAM_FILE} -M ${DUP_FILE_QC} --VALIDATION_STRINGENCY LENIENT --ASSUME_SORTED true --REMOVE_DUPLICATES false --QUIET true --VERBOSITY ERROR 2> "${NAME}.picard.MarkDuplicates.error.log"

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Filtering 780, Fixming mate, Sorting and indexing..."
# ===========================================================================================
# Remove unmapped, mate unmapped
# not primary alignment, reads failing platform
# Obtain name sorted BAM file
# Remove orphan reads (pair was removed) and read pairs mapping to different chromosomes
# Obtain position sorted BAM with flagged duplicates
# Index final position sorted BAM
# Create final name sorted BAM with flagged duplicates for featureCounts
# ===================================================================================================

START_SUBPROCESS=$(date +%s)

samtools view -F 780 -Shb ${MARKDUP_BAM_FILE} -@ "${OMP_NUM_THREADS}" | samtools sort -n -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" - -o ${NMSRT_MARKDUP_QFILT_BAM_FILE}
samtools fixmate -r -@ "${OMP_NUM_THREADS}" ${NMSRT_MARKDUP_QFILT_BAM_FILE} - | samtools view -@ "${OMP_NUM_THREADS}" -F 780 -Shb - | samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" - -o ${FINAL_MARKDUP_BAM_FILE}
samtools index ${FINAL_MARKDUP_BAM_FILE} ${FINAL_MARKDUP_INDEX_FILE}

rm ${TMP_NMSRT_BAM_FILE} ${TMP_QFILT_BAM_FILE} ${TMP_QFILT_BED_FILE} ${SAM_FILE_TMP} ${BAM_FILE_CHR_FILETERED_TMP} ${NMSRT_MARKDUP_QFILT_BAM_FILE} ${MARKDUP_BAM_FILE} "${RAW_BAM_FILE}.bai"

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."


echo "Filtering 1804: remove duplicates, remove unmapped, mate unmapped not primary aligment and reads failing platform/vendor quiality checks..."
# ===========================================================================================
# Remove unmapped, mate unmapped
# not primary alignment, reads failing platform
# Only keep properly paired reads
# Obtain name sorted BAM file
# Remove orphan reads (pair was removed) and read pairs mapping to different chromosomes
# Obtain position sorted BAM
# Remove duplicates
# Index final position sorted BAM
# Create final name sorted BAM
# ===================================================================================================

START_SUBPROCESS=$(date +%s)

samtools view -@ "${OMP_NUM_THREADS}" -F 1804 -f 2 -Sbh ${FINAL_MARKDUP_BAM_FILE} | samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" - -o ${NODUP_BAM_FILE} # Will produce coordinate sorted BAM / Final NODUP BAM file

echo "Indexing..."
samtools index ${NODUP_BAM_FILE} ${FINAL_BAM_INDEX_FILE}		# Index Final NODUP BAM file

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Name Sorting..."
START_SUBPROCESS=$(date +%s)

samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" -n ${NODUP_BAM_FILE} -o ${FINAL_NMSRT_BAM_FILE}	# Create final name sorted NODUP_BAM
samtools sort -@ "${SAM_THREADS}" -m "${MEM_SAMTOOLS}G" -n ${FINAL_MARKDUP_BAM_FILE} -o ${FINAL_NMSRT_MARKDUP_BAM_FILE}	# Create final name sorted MARKDUP_BAM

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."


module unload SAMtools R picard
module load deepTools

echo "Deeptools QC..."
START_SUBPROCESS=$(date +%s)

plotFingerprint -b ${FINAL_MARKDUP_BAM_FILE} -plot ${CHART_FINGERPRINTS} --outRawCounts ${FINGERPRINTS} --skipZeros -T "Fingerprints"
plotCoverage -b ${FINAL_MARKDUP_BAM_FILE} -o ${CHART_COVERAGE} --ignoreDuplicates -T "CoveragePlot"
estimateReadFiltering -b ${RAW_BAM_FILE} -o ${ESTIMATE_READ_FILTERING} -p "${OMP_NUM_THREADS}" # -bl ${BLACKLIST}  do not use blacklist in RNA-seq

echo "Deetools Coverages Profiles: RPKM (per bin) = number of reads per bin / (number of mapped reads (in millions) * bin length (kb))..."
# bin size 10 bp
# The smooth length defines a window, larger than the binSize, to average the number of reads. For example, if the –binSize is set to 20 and the –smoothLength is set to 60, then, for each bin, the average of the bin and its left and right neighbors is considered.

bamCoverage -b ${FINAL_MARKDUP_BAM_FILE} -o ${BIGWIG_COVERAGE} -bs 10 -e --normalizeUsing RPKM  -p "${OMP_NUM_THREADS}" --samFlagExclude 780 # -bl ${BLACKLIST} do not use blacklist in RNA-seq

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

echo "Deeptools plotProfile..."
START_SUBPROCESS=$(date +%s)

computeMatrix scale-regions -R ${GENE} -S ${BIGWIG_COVERAGE} -o ${MATRIX_SCALED_GENES} --regionBodyLength 1000 -b 3000 -a 3000 --missingDataAsZero  --skipZeros -p "${OMP_NUM_THREADS}" -q # -bl ${BLACKLIST} do not use blacklist in RNA-seq
plotProfile -m ${MATRIX_SCALED_GENES} -out ${PLOTPROFILE_PLOT} --outFileNameData ${PLOTPROFILE_DATA} --plotTitle "Read Distribution Profile"

RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

module unload deepTools
module load SAMtools Homer

echo "HOMER: Making TagDirectory..."  
START_SUBPROCESS=$(date +%s)

makeTagDirectory ./${NAME}_TagDirectory ${FINAL_MARKDUP_BAM_FILE} -single -tbp 1 -genome ${REFERENCE} -checkGC

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

module unload SAMtools Homer
module load Subread

# Counting reads with featureCounts
echo "****Running featureCounts****"

START_SUBPROCESS=$(date +%s)
featureCounts -p --countReadPairs -t exon -T "${OMP_NUM_THREADS}" -a "${GENE}" -o "${COUNTS_FULL_TABLE}" ${FINAL_NMSRT_MARKDUP_BAM_FILE}
# awk '{print $1"\t"$7}' ${COUNTS_FULL_TABLE} > ${FINAL_COUNTS}	# Two Columns GeneID/Raw Counts by FeatureCounts

cut --complement -f2-6 ${COUNTS_FULL_TABLE} | tail -n +2 > "${FINAL_COUNTS}" # Two Columns GeneID/Raw Counts by FeatureCounts and remove the first line (Feature Count header)

END_SUBPROCESS=$(date +%s)
RUNTIME_SUBPROCESS=$((END_SUBPROCESS-START_SUBPROCESS))
H=$((RUNTIME_SUBPROCESS / 3600 ))  # Calculate hours
M=$(((RUNTIME_SUBPROCESS / 60 ) % 60 ))  # Calculate minutes
S=$((RUNTIME_SUBPROCESS % 60 ))  # Calculate seconds
echo -e "Status: Done! Used ${H} hours, ${M} minutes, and ${S} seconds."

module unload Subread

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

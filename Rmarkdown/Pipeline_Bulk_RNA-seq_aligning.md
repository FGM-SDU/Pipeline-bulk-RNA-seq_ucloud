Bulk RNA-seq Alignment Pipeline
================
Victor Enrique Goitea
2024-03-04

- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Usage](#usage)
- [Output](#output)

# Overview

This bash script is designed for running a RNA-seq alignment on a SLURM
cluster. The pipeline includes read alignment using STAR,
post-processing with SAMtools and Picard. Optionally it can also do
psuedoalignment and quantification using Salmon with bootstrap for
isoform level DGE analysis.

# Prerequisites

- SLURM scheduler
- Modules: STAR, SAMtools, picard, Salmon
- Reference genomes (provided in /work/References)

# Usage

1.  **Modify SLURM Parameters (Optional):** Open the script
    (**pe_align_rnaseq_v2_multigenome.sh**) and modify SLURM parameters
    at the beginning of the file, such as account, output file, email
    notifications, nodes, memory, CPU cores, and runtime. Alternatively,
    you can modify these parameters on-the-fly when executing the
    script.

2.  On UCloud, start a **Terminal Ubuntu** run:

    - Enable **Slurm cluster**
    - - To process several samples consider requesting nodes \> 1
    - Set the modules path to **FGM \> Utilities \> App \> easybuild**

    ![](../Img/terminal_slurm.png)

    - Include the References folder **FGM \> References \> References**

    ![](../Img/terminal_folders.png)

    - Include your Scripts folder and the folder with the fastq.gz
      files.

    - **Notes:**

      - Match the job CPUs to the amounts requested in the script.
      - Make sure the scripts have executing permission. If not run:
        `chmod 700 script.sh`
      - If you modify the memory parameter in the script, specify 5-10%
        less than the memory available in the terminal run.
      - Although it is not necessary to enable **tmux**, it is a good
        practise to always do it.

3.  **Run the Script:** Submit the script to the SLURM cluster:

        sbatch -J <job_name> path_to/Scripts_folder/pe_align_rnaseq_v2_multigenome.sh -g <mm10|mm39|hg38> <input-R1_001.fastq.gz-file> -t

    **Required Arguments**

    - **-g:** specify the genome to use (mm10, mm39, or hg38).
    - Replace **input-R1_001.fastq.gz-file** with the full path to your
      input FASTQ file (R1).

    **Optional Arguments**

    - **-t:** include Salmon alignment for isoform-level differential
      gene expression. Include the flag “-t” after the input file.

    For several samples you can use a for loop:

        for i in *R1_001.fastq.gz; do sbatch -J <job_name> path_to/Scripts_folder/pe_align_rnaseq_v2_multigenome.sh -g <mm10|mm39|hg38> $i; sleep 1; done

4.  **Monitor Job:** You can monitor the job using the SLURM commands,
    such as squeue, scontrol show job <job-id>, and check the log files
    generated.

# Output

The script generates various output files, including aligned BAM files,
log files, and optional Salmon output for isoform-level analysis.

**Notes:**  
- Ensure that the necessary modules are available on your cluster.  
- The script includes Slurm directives to specify resource requirements.
Review and customize the script based on your specific requirements.  
- For additional information on individual tools and parameters, refer
to the documentation for STAR, SAMtools, Picard, and Salmon.

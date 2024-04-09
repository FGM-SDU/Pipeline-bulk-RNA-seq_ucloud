# Bulk RNA-seq pipeline UCloud/Slurm 1.1

* The pipeline evaluates first strandedness using RSeQC (infer_experiment.py) and based on the results it automatically proceeds accordingly, i.e. reverse-stranded, forward-stranded or non-stranded setups. Tools strandness aware include Qualimap RNA-SeQC and featureCounts.


# Bulk RNA-seq pipeline UCloud/Slurm 1.0

* Initial release of bulk-RNA-seq pipeline

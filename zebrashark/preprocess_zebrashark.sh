#!/bin/bash
#SBATCH --partition=batch
#SBATCH --mail-user=samantha.bock@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=fastqc_zebrashark_5.0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/work/dir/scripts/fastqc_zebrashark_5.0.o
#SBATCH --error=/work/dir/scripts/fastqc_zebrashark_5.0.e


## go to directory with raw reads
cd /work/dir/work/ZEBRA_SHARK/RAW_READS/METHYLSEQ_121923

## load fastqc
module load FastQC/0.11.9-Java-11
module load MultiQC/1.14-foss-2022a

for i in GN22895 GN22896 GN22897 GN22898 GN22899 GN22900 GN22901 GN22902 GN22903 GN22904 GN22905 GN22906 GN22907 GN22908 GN22909 GN22910 GN22911 GN22912 GN22913 GN22914 GN22915 GN22916 GN22917 GN22918 GN22919 GN22920 GN22921 GN22922 GN22923 GN22924 GN22925 GN22926 GN22927 GN22928
do
  fastqc -t 8 ${i}_*_L00*_R1_001.fastq.gz -o /work/dir/work/ZEBRA_SHARK/Pretrim_fastqc/METHYLSEQ_121923
  fastqc -t 8 ${i}_*_L00*_R2_001.fastq.gz -o /work/dir/work/ZEBRA_SHARK/Pretrim_fastqc/METHYLSEQ_121923
done


cd /work/dir/work/ZEBRA_SHARK/Pretrim_fastqc/METHYLSEQ_121923
multiqc -d *_fastqc.zip -o /work/dir/work/ZEBRA_SHARK/Pretrim_fastqc/METHYLSEQ_121923



# guide to FastQC interpretation: https://training.galaxyproject.org/training-material/topics/sequence-analysis/tutorials/quality-control/tutorial.html

#!/bin/bash
#SBATCH --partition=batch
#SBATCH --mail-user=samantha.bock@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=GN22898_trim_zebrashark_1.0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/sb61937/scripts/GN22898_trim_zebrashark_1.0.o
#SBATCH --error=/scratch/sb61937/scripts/GN22898_trim_zebrashark_1.0.e


module load Trim_Galore/0.6.7-GCCcore-11.2.0
cd /scratch/sb61937/work/ZEBRA_SHARK/RAW_READS/METHYLSEQ_121923

# trim adapter sequences from emseq reads, setting strigency to 1 since we've used this in the past for methylation data
# adding an additional trim of 4 bp based on past m-bias plots

for i in GN22898_S35_L005 GN22898_S35_L007 GN22898_S4_L004 GN22898_S4_L006
do
  trim_galore --fastqc --fastqc_args "--outdir /scratch/sb61937/work/ZEBRA_SHARK/TrimGalore_Fastqc/METHYLSEQ_121923 -f fastq" -stringency 1 --cores 8 --clip_r1 4 --clip_r2 4 --three_prime_clip_r1 4 --three_prime_clip_r2 4 -o /scratch/sb61937/work/ZEBRA_SHARK/TrimGalore_Fastqc/METHYLSEQ_121923 --keep --paired --retain_unpaired /scratch/sb61937/work/ZEBRA_SHARK/RAW_READS/METHYLSEQ_121923/${i}_R1_001.fastq.gz /scratch/sb61937/work/ZEBRA_SHARK/RAW_READS/METHYLSEQ_121923/${i}_R2_001.fastq.gz
done


## load fastqc
#module load FastQC/0.11.9-Java-11
#module load MultiQC/1.14-foss-2022a
## go to output directory
#cd /scratch/sb61937/work/ZEBRA_SHARK/TrimGalore_Fastqc/METHYLSEQ_120523
## run multiqc on the output files from fastqc
#multiqc -d *_fastqc.zip -o /scratch/sb61937/work/ZEBRA_SHARK/TrimGalore_Fastqc/METHYLSEQ_120523

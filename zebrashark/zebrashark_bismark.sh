#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --mail-user=samantha.bock@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=GN22898_align_dedup_st4_1.0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=120:00:00
#SBATCH --output=/scratch/sb61937/scripts/GN22898_align_dedup_st4_1.0.o
#SBATCH --error=/scratch/sb61937/scripts/GN22898_align_dedup_st4_1.0.e

cd /scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome

# load samtools
module load SAMtools/1.16.1-GCC-11.3.0

# load Bismark
module load Bismark/0.24.1-GCC-11.3.0

# load Bowtie2
module load Bowtie2/2.4.5-GCC-11.3.0

#running Bismark
cd /scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/METHYLSEQ_121923

for i in GN22898_S35_L005 GN22898_S35_L007 GN22898_S4_L004 GN22898_S4_L006
do
  bismark --parallel 4 --bowtie2 --maxins 1000 --genome /scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome -1 /scratch/sb61937/work/ZEBRA_SHARK/TrimGalore_Fastqc/METHYLSEQ_121923/${i}_R1_001_val_1.fq.gz -2 /scratch/sb61937/work/ZEBRA_SHARK/TrimGalore_Fastqc/METHYLSEQ_121923/${i}_R2_001_val_2.fq.gz
done

cd /scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome

for i in GN22898
do
  deduplicate_bismark -o ${i} --bam --paired --multiple /scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/METHYLSEQ_120523/${i}_S*_L004_R1_001_val_1_bismark_bt2_pe.bam /scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/METHYLSEQ_121923/${i}_S*_L004_R1_001_val_1_bismark_bt2_pe.bam /scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/METHYLSEQ_121923/${i}_S*_L005_R1_001_val_1_bismark_bt2_pe.bam /scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/METHYLSEQ_121923/${i}_S*_L006_R1_001_val_1_bismark_bt2_pe.bam /scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/METHYLSEQ_121923/${i}_S*_L007_R1_001_val_1_bismark_bt2_pe.bam
done

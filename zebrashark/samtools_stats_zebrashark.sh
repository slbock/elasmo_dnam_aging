#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --mail-user=samantha.bock@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=GN22895-GN22928_samtools_zebrashark_1.0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=168:00:00
#SBATCH --output=/scratch/sb61937/scripts/GN22895-GN22928_samtools_zebrashark_1.0.o
#SBATCH --error=/scratch/sb61937/scripts/GN22895-GN22928_samtools_zebrashark_1.0.e


cd /scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/

# load samtools
module load SAMtools/1.16.1-GCC-11.3.0 # note you can find the most recent version of samtools on the cluster by using the command 'module spider samtools'

# run AFTER deduplication
# sort bam files
for i in GN22895 GN22896 GN22897 GN22898 GN22899 GN22900 GN22901 GN22902 GN22903 GN22904 GN22905 GN22906 GN22907 GN22908 GN22909 GN22910 GN22911 GN22912 GN22913 GN22914 GN22915 GN22916 GN22917 GN22918 GN22920 GN22921 GN22922 GN22923 GN22926 GN22927 GN22928
do
  samtools sort /scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/${i}.multiple.deduplicated.bam -o /scratch/sb61937/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome/Sorted_Bam/${i}_sort.bam
done

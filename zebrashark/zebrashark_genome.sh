#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --mail-user=samantha.bock@uga.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=bismark_genomeprep_zebrashark_1.0
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=72G
#SBATCH --time=120:00:00
#SBATCH --output=/dir/scripts/bismark_genomeprep_zebrashark_1.0.o
#SBATCH --error=/dir/scripts/bismark_genomeprep_zebrashark_1.0.e

cd /dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome

# download  Stegostoma tigrinum genome and annotation
# curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_030684315.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCF_030684315.1.zip" -H "Accept: application/zip"
# unzip GCF_030684315.1.zip

# load Bismark
module load Bismark/0.24.1-GCC-11.3.0

# load Bowtie2
module load Bowtie2/2.4.5-GCC-11.3.0

# Bismark performs alignments of bisulfite-treated reads to a reference genome
# peforms cytosine methylation calls at the same time

#first, the genome needs to be bisulfite converted and indexed to allow Bowtie alignments
#prepare genome, note bowtie1 and bowtie2 require distinct indexing steps since their indexes are not compatible
#genome preparation step only needs to be preformed once
bismark_genome_preparation --bowtie2 --verbose /dir/work/ZEBRA_SHARK/Bismark/sSteTig4_Genome

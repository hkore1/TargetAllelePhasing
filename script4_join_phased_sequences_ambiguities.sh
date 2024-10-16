# List of genes and samples
genelist=/data/gpfs/projects/punim1533/corymbia_phylogeny/6_analyses/nDNA/4_allele_phasing/genelist.txt
namelist=/data/gpfs/projects/punim1533/corymbia_phylogeny/6_analyses/nDNA/4_allele_phasing/_sample_list.txt


# Check to make sure script is being run from correct directory
if [ ! -d 2_Phased_Sequences/ ]; then
  echo "Directory '2_Phased_Sequences' does not exist! Exiting..."
  exit 1
fi

#set -eo pipefail
module load java/17.0.6
module load parallel/20220722

# Create new directory.
if [ ! -d 3_alignments/ ]; then
	mkdir -p 3_alignments
	rm -f 3_alignments/*
fi


###### Exon sequences generated from HybPiper output:

mkdir -p 3_alignments/exon
rm -f 3_alignments/exon/* 
parallel "cat 2_Phased_Sequences/iupac/{1}_IUPAC/{1}_exon/{2}.iupac.FNA >> 3_alignments/exon/{2}.exon.iupac.fasta" :::: $namelist :::: $genelist
wait


###### Intron sequences generated from HybPiper output:

mkdir -p 3_alignments/intron
rm -f 3_alignments/intron/* 
parallel "cat 2_Phased_Sequences/iupac/{1}_IUPAC/{1}_intron/{2}.intron.iupac.fasta >> 3_alignments/intron/{2}.intron.iupac.fasta" :::: $namelist :::: $genelist
wait
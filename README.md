# TargetAllelePhasing
Pipeline for generating phased allele sequences from target capture data after processing with HybPiper.
Based on the workflow detailed here (https://github.com/mossmatters/phyloscripts/tree/master/alleles_workflow), but updated to use newer software versions which break the original workflow, and with re-formatted scripts to make it easier to run on any given batch of data.

## Dependencies:
WhatsHap: https://whatshap.readthedocs.io/en/latest/

GATK: https://gatk.broadinstitute.org/hc/en-us

Picard: https://github.com/broadinstitute/picard

BCFtools: https://github.com/samtools/bcftools

SAMtools: https://github.com/samtools/samtools

MACSE: https://bioweb.supagro.inra.fr/macse/

trimAl: http://trimal.cgenomics.org

----------------------------------------------------
All of these programs (except for WhatsHap) are packaged and compiled to be downloaded from here.

WhatsHap should be installed using Anaconda, and the lines to activate conda environment containing WhatsHap should be altered in the wrapper 'script3_phase_extract_align.sh' or 'subscript3.1_phase_alleles_with_whatshap.sh' (the latter if you're running WhatsHap independently of other scripts).

## Before running the pipeline
HybPiper must have already been run, using the --run_intronerate function to generate supercontigs with introns (supercontig.fasta files) and an intronerate.gff file (that stores exon/intron boundary information) for each gene.

We use a wrapper script for HybPiper2 (based on https://github.com/Royal-Botanic-Gardens-Victoria/GAP_hybpiper_helpers/blob/main/hp2.py) that compresses these two files (plus many others) into a 'tar.gz' file.

If you have run HybPiper using a method that does not tar these files, then comment out line 37 of the script 'script1_concat_hp2_supercontigs_and_introns_rename.sh' ("`tar -xf "${dir%/}/$prefix.tar.gz`").



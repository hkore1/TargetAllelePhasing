# TargetAllelePhasing
Pipeline for generating phased allele or ambiguity sequences from target capture data after processing with HybPiper.
Based on the workflow detailed here (https://github.com/mossmatters/phyloscripts/tree/master/alleles_workflow), but updated to use newer software versions which break the original workflow, and with re-formatted scripts to make it easier to run on any given batch of data.

## Dependencies:
WhatsHap (v1.7): https://whatshap.readthedocs.io/en/latest/

GATK (v4.3.0.0): https://gatk.broadinstitute.org/hc/en-us

Picard (v2.27.5): https://github.com/broadinstitute/picard

BCFtools (v1.15): https://github.com/samtools/bcftools

SAMtools (v1.16.1): https://github.com/samtools/samtools

trimAl (v1.2rev59): http://trimal.cgenomics.org

----------------------------------------------------
WhatsHap should be installed using Anaconda, and the lines to activate conda environment containing WhatsHap should be altered in the wrapper 'script3_phase_extract_align.sh' or 'subscript3.1_phase_alleles_with_whatshap.sh' (the latter is only necessary if you're running WhatsHap independently of other scripts).

Please note that some scripts contain lines for initialising modules specific to the HPC system that we used. These will likely need to be adjusted for other systems.


## Before running the pipeline
HybPiper must have already been run, using the `--run_intronerate` function to generate supercontigs with introns (supercontig.fasta files) and an intronerate.gff file (that stores exon/intron boundary information) for each gene.

We use a wrapper script for HybPiper2 (based on https://github.com/Royal-Botanic-Gardens-Victoria/GAP_hybpiper_helpers/blob/main/hp2.py) that compresses these two files (plus many others) into a 'tar.gz' file.

If you have run HybPiper using a method that does not tar these files, then comment out line 37 of the script 'script1_concat_hp2_supercontigs_and_introns_rename.sh' ("`tar -xf "${dir%/}/$prefix.tar.gz`").

## Step 1: Concatenate HybPiper2 supercontig and intronerate files
Run `bash script1_concat_hp2_supercontigs_and_introns_rename.sh <hp2_results_directory>`. 

Note: 'hp2_results_directory' is a directory containing HP2 output folders for each sample in your dataset.

## Step 2: Map reads to supercontig sequences to generate new IUPAC-coded ambiguity sequences
Run `bash script2_map_to_supercontigs.sh <supercontig_directory> <read_directory> <number_of_CPUs>`.

Note: 'supercontig_directory' is the 'supercontigs' folder generated during Step 1. 'read_directory' is the directory containing reads used for the initial HybPiper assembly.

This step requires multiple CPUs and is rather slow (~ 1 hr per sample average).

## Step 3: Generate phased sequences and extract exons, introns and supercontigs
Run `bash script3_phase_extract.sh <full/path/to/gene_list_file>`

This is a wrapper script to perform functions in 4 subscripts that:
1. Phase alleles from IUPAC reference supercontigs (with WhatsHap)
2. Generate separate fasta files for phased sequences (with bcftools)
3. Replace variant sites outside the longest phase block with ambiguity codes (with haplonerate.py)
4. Generate separate files for the intron, exon, and supercontig sequences for samples for default HybPiper, IUPAC-coded, and phased-allele data (with intron_exon_extractor.py)

Notes: "gene_list_file" is either a .txt file (pre made by the user) containing the list of gene names, or the user can specify the string "Angiosperms353" to generate this file automatically if they are using the Angiosperms353 genes. Make sure to set the variable `Step4_subscript_path` as a path to 'subscript3.4_intron_exon_extractor.py'

It is probably preferable to run these subscripts individually whilst setting up and/or troubleshooting the pipeline.

## Step 4: Combine, align and trim
Run 'script4_join_phased_sequences_ambiguities.sh' from the directory containing the '2_Phased_Sequences' folder. This joins the IUPAC-coded sequences of all samples into multifasta files for each locus.

These files can then be aligned and trimmed/cleaned however is preferred. The steps we used for this are detailed in the bloodwood eucalypt paper.

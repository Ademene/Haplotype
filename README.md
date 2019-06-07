# Haplotype

This text describe the method for the « Identification of the different haplotypes among the French clonal lineages » mentionned in the article "Whole-genome sequencing reveals recent and frequent genetic recombination between clonal lineages of Cryphonectria parasitica in western Europe" Demené et al., 2019 (https://doi.org/10.1016/j.fgb.2019.06.002). This method is not automated. 

## 1) Remove singletons from the initial vcf using bcftools{1}.

$ bcftools view -e 'HOM=1' input.vcf.gz -o input_WithoutSingletons.vcf

## 2) Extract the SNPs number between each couple of individuals considered along sliding windows of the desired size using vcftools{2}. 
### Example for three individuals ind1, ind2 and ind3 ; for a window size of 10kb. Note : For n individuals, number of files produced is (n-1)! . 
### I did not automated this step.

$ vcftools --vcf input_WithoutSingletons.vcf --indv ind1 --indv ind2 --window-pi 10000 --out ind1_ind2

$ vcftools --vcf input_WithoutSingletons.vcf --indv ind1 --indv ind3 --window-pi 10000 --out ind1_ind3

$ vcftools --vcf input_WithoutSingletons.vcf --indv ind2 --indv ind3 --window-pi 10000 --out ind2_ind3

## 3) Create a genome file as follow: 

#### Contig Length
#### MS1-1 3855066
#### MS1-2 554633

### Here is an example of python script wich works for a single line multifasta file:
### You first need to unwrap multi-lines fasta: 
$ awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < YourGenome.fasta
### Python script (extension of your fasta file must be .fasta):
$ python Fasta_to_GenomeFile.py YourGenome.fasta

### The file lengthCOntigPacbio.genome located in the folder Haplotype_Test is the genome file I used.

## 4) Run R{3} script to make the plot. Here, I provide the script designed to plot the haplotype figure S3 (10 individuals, 45 files). This script should be adapted manualy for use with other data sets.

### Download the folder Haplotype_Test and the script Haplotypes_script.R
### Change the work directory in the R script to link to Haplotype_Test folder.
### Run the whole script in R to obtain the Figure S3.

## References

{1} Copyright (c) Genome Research Ltd. http://samtools.github.io/bcftools/ 
Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. [PMID: 19505943] 

{2} https://github.com/vcftools/vcftools
The Variant Call Format and VCFtools, Petr Danecek, Adam Auton, Goncalo Abecasis, Cornelis A. Albers, Eric Banks, Mark A. DePristo, Robert Handsaker, Gerton Lunter, Gabor Marth, Stephen T. Sherry, Gilean McVean, Richard Durbin and 1000 Genomes Project Analysis Group, Bioinformatics, 2011

{3} R Core Team (2014). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.   URL http://www.R-project.org/.

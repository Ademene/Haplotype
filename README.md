# Haplotype

This text describe the method for the « Identification of the different haplotypes among the French clonal lineages » mentionned in the article "Whole-genome sequencing reveals recent and frequent genetic recombination between clonal lineages of Cryphonectria parasitica in western Europe" Demené et al., 2019. This method is not automated. 

### 1) Remoove singletons from the initial vcf using bcftools.

$ bcftools view -e 'HOM=1' input.vcf.gz -o input_WithoutSingletons.vcf

### 2) Extract the SNPs number between each couple of individuals considered along sliding windows of the desired size using vcftools. Example for three individuals ind1, ind2 and ind3 ; for a window size of 10kb. Note : For n individuals, number of files produced is (n-1)! . I did not automated this step.

$ vcftools --vcf input_WithoutSingletons.vcf --indv ind1 --indv ind2 --window-pi 10000 --out ind1_ind2
$ vcftools --vcf input_WithoutSingletons.vcf --indv ind1 --indv ind3 --window-pi 10000 --out ind1_ind3
$ vcftools --vcf input_WithoutSingletons.vcf --indv ind2 --indv ind3 --window-pi 10000 --out ind2_ind3

### 3) Create a genome file as follow : 

## Contig Length
## MS1-1 3855066
## MS1-2 554633
## …
## Here is an example of python script wich works for a single line multifasta file :
# If you need to unwrap your multiple lines fasta file you can use : 
$ awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < file.fa
## Python script

### 4) Open all the files with the R script to make the plot. Here is the script designed to plot the haplotype figure S3 (10 individuals, 45 files)

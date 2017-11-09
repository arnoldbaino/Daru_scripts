# Daru_scripts
Scripts for manuplating raw illumina data 

The following tools should be installed and added to the system path
  
Trimmomatic
  
BWA
  
Samtools
  
bcftools

Picard
  
GATK
  
VCFTools
  
SnpEff
  
vcf_filter_module.py
  
  
The following commands can be for each sample data or adjusted to run all the sample data

1. Trim the adapters and low quality data at the 3' 

java -jar /home/trimmomatic/current/bin/trimmomatic-0.32.jar PE -phred33 -trimlog file_log.txt file_R1.fastq.gz file_R2.fastq.gz file_paired_R1.fastq.gz file_unpaired_R1.fastq.gz file_paired_R2.fastq.gz file_unpaired_R2.fastq.gz ILLUMINACLIP:/home/trimmomatic/current/bin/adapters/NexteraPE-PE.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36

2. Select a reference genome (i.e H37Rv, NC_000962.3) and map trimmed readsto the reference. Known variable sites (PE/PPE etc) in the genome can be masked prior to mapping raw reads

bwa mem -t 4 -M -R "@RG\tID:OOOO\tSM:file.name\tPL:Illumina\tLB:001\tPU:001" reference.fasta file_R1.fastq.gz file_R2.fastq.gz > file.sam

3. Change Sam file to Bam file 

samtools view -bS -F 4 file.sam > file.bam

4.Sort the bam file

java -Xms1000M -Djava.io.tmpdir=./tmp. -jar picard-tools-1.119/SortSam.jar I=file.bam O=file_sorted.bam SO=coordinate

5. Mark duplicate reads

java -Xms1000M -Djava.io.tmpdir=./tmp. -jar picard-tools-1.119/MarkDuplicates.jar I=file_sorted.bam O=file_dup.bam M=file.txt

6\. Generate index for each bam file
  
  java -Xms1000M -Djava.io.tmpdir=./tmp. -jar picard-tools-1.119/BuildBamIndex.jar I=file_dup.bam

7. Define targets for local realignment

java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R reference.fasta -I file_dup.bam -o file.intervals

8. Perform realignment around Indels

java -jar GenomeAnalysisTK.jar -T IndelRealigner -R reference.fasta -I file_dup.bam -targetIntervals file.intervals -o file_dup_alig.bam

9. Call variants on all bam files (files.list) 

java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R reference.fasta -I files.list -A AlleleBalance -pnrm EXACT_GENERAL_PLOIDY -ploidy 1 -glm SNP -o files.vcf

10. Filters SNPs with (i) less than 10 unambigous read alignments (ii) 80% or more of the read depth had ambigous mappings (iii) fewer than 80% of reads supported the alternative allele

java -jar GenomeAnalysisTK.jar -T VariantFiltration -R reference.fasta -V files.vcf --filterExpression "((DP-MQ0)<10) || ((MQ0/(1.0*DP))>=0.8) || (ABHom <0.8) || (QUAL > 90)" --filterName LowConfidence -o files_filtered.vcf

11. Remove SNPs with 'lowconfidence'

vcftools --vcf files_filtered.vcf --recode --keep-INFO-all

12. To filter SNPs within 10 base distance of each other, the 'master vcf' is split up into individual sample vcf and then use vcf_filter_module.py on each file

for file in files_filtered.vcf;do for sample in `/home/bcftools/bcftools-git/bcftools view -h $file | grep "^#CHROM" | cut  -f10-`; do /home/bcftools/bcftools-git/bcftools view -c1 -s $sample -o ${file/.vcf*/.$sample.vcf} $file;done;done

then use: vcf_filter_module.py 9 in.vcf out.vcf

13. Combine the vcfs into one

java -jar GenomeAnalysisTK.jar -R reference.fasta -T CombineVariants –V vcf1 vcf2 vcf3.. -genotypeMergeOptions UNIQUIFY –o master.vcf

14. Annotate the SNPs in the master vcf file and predict there effect on the genes

java -Xmx4g -jar snpEff.jar -c snpEff.config -v m_tuberculosis_H37Rv file.vcf > annotated_file.vcf

15. Change annotated_file.vcf into a fasta file by concatenating





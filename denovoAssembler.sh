# This is a complete chloroplast genome assembler that takes two whole genome shotgun .fastq files as input.  It generates a short-read
# alignment map (.sam) of these short-reads aligned to the reference, with SOAP kmer value as input.  (NOTE!  This script can be used for
# other assemblies, but the soapdenovo2-63mer step must be modified to change the estimated genome size.  It is currently set for chloroplast
# size of 150000)

# example usage (NOTE!  Input fastqs must be lacking their .fastq or .fq extention!):
# bash denovoAssembler.sh I_cyaneum170_dan1_ATTCAGAA-CAGGACGT_L001_R1_001 I_cyaneum170_dan1_ATTCAGAA-CAGGACGT_L001_R2_001 Solanum_tuberosum.fasta 31

#left- and right- fastq reads
library_input1=$1
library_input2=$2

#reference
ref=$3

#soapdenovo2-63mer kmer value
kmer_value=$4

#Index the reference
bwa index $ref

#Get reference ID for to pull out reads from the .sam file that map to the reference
refID=`head -n 1 Solanum_tuberosum_cp.fasta | sed 's/>//g' | cut -f1 -d " "`

#generate first-pass alignment of whole-genome shotgun data onto potato chloroplast genome.  Most reads in our 
#library are nuclear genome reads, and will not align in any way to the potato chloroplast.  We want to determine 
#the reads that do.
bwa mem -k 13 -B 2 -O 2 -L 3 $ref ${library_input1}.fastq ${library_input2}.fastq > ${library_input1}.sam

echo "bwa mem step 1: Complete"

#now create left and right .fastq files for the mapped chloroplast reads
cat ${library_input1}.sam | grep $refID | sort -k10,1 | awk '{print $1}' | sort | uniq -c | sort -nr | awk '$1 == 2 {print $2}' > mappedReadsList.txt
grep -A3 --no-filename --no-group-separator -Ff mappedReadsList.txt ${library_input1}.fastq > mappedReads_1.fastq
grep -A3 --no-filename --no-group-separator -Ff mappedReadsList.txt ${library_input2}.fastq > mappedReads_2.fastq

echo "Extracting chloroplast reads:  Complete"


#####Legacy code to generate SNP table#####
#re-do the bwa mem alignment to generate a binary call file (bcf) with mpileup.  This then generates a variant call file (.vcf)
#bwa mem $ref mappedReads_1.fastq mappedReads_2.fastq > ${library_input1}_mapped.sam
#samtools view -b -o ${library_input1}_mapped.bam -S ${library_input1}_mapped.sam
#samtools sort -m 5000000000 ${library_input1}_mapped.bam ${library_input1}_mapped.sorted
#samtools index ${library_input1}_mapped.sorted.bam
#samtools faidx Solanum_tuberosum_cp.fasta
#samtools mpileup -uf Solanum_tuberosum_cp.fasta ${library_input1}_mapped.sorted.bam | bcftools view -bvcg -> ${library_input1}_mapped.bcf
#bcftools view ${library_input1}_mapped.bcf > ${library_input1}_snps_indels.vcf

#echo "Successfully generated .vcf"
#####End legacy code#####

path=$PWD

#Reconfigure SOAP config file with path to mapped reads
sed "s;QQQQQ1;$path\/mappedReads_1.fastq;g" soap.config | sed "s;QQQQQ2;$path\/mappedReads_2.fastq;g" > soap.config.temp

#Run SOAP!
soapdenovo2-63mer all -s soap.config.temp -K $kmer_value -R -o K${kmer_value} -z 150000 1>assembly.log 2>assembly.err

#Move newly created files into own directory
mkdir K${kmer_value}
mv K${kmer_value}* K${kmer_value}/

#####Legacy code that probably no longer works#####
#soapdenovo2-63mer sparse_pregraph -s soap.config4 -K $kmer_value -R -o K$kmer_value -z 150000 > soaplog_sparsegraphK$kmer_value
#soapdenovo2-63mer contig -R -g K$kmer_value -M 2 -E > soaplog_contigK$kmer_value
#soapdenovo2-63mer map -s soap.config4 -g K$kmer_value -f -k $kmer_value > soaplog_mapK$kmer_value
#soapdenovo2-63mer scaff -g K$kmer_value -F -k $kmer_value -G 500 > soaplog_scaffK$kmer_value

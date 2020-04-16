#!/bin/bash

#R1="MM30-Barcode-ad1-ligated-CN18_R1.fastq.gz"
#R2="MM30-Barcode-ad1-ligated-CN18_R2.fastq.gz"
R1="test_R1.fastq.gz"
R2="test_R2.fastq.gz"

cutadapt \
    -e 0.12 \
    -a GCGGCCGCT \
    --rc \
    --action=none \
    --untrimmed-output untrimmed_R1.fastq.gz \
    -o output_R1.fastq.gz \
    $R1 2>&1>report_R1.txt
    #MM30-Barcode-ad1-ligated-CN18_R1.fastq.gz 2>&1>report_R1.txt

cutadapt \
    -e 0.12 \
    -a GCGGCCGCT \
    --rc \
    --action=none \
    --untrimmed-output untrimmed_R2.fastq.gz \
    -o output_R2.fastq.gz \
    $R2 2>&1>report_R2.txt
    #MM30-Barcode-ad1-ligated-CN18_R2.fastq.gz 2>&1>report_R2.txt


# get the union set of reads that have tags in them and then split out the original FASTQ:
python3 extract_readnames.py -r1 output_R1.fastq.gz -r2 output_R2.fastq.gz -o readnames.list
/opt/software/bbmap/filterbyname.sh in=$R1 names=readnames.list include=t out=tagged_original_seqs_R1.fastq
/opt/software/bbmap/filterbyname.sh in=$R2 names=readnames.list include=t out=tagged_original_seqs_R2.fastq
/opt/software/bbmap/filterbyname.sh in=$R1 names=readnames.list include=f out=untagged_original_seqs_R1.fastq.gz
/opt/software/bbmap/filterbyname.sh in=$R2 names=readnames.list include=f out=untagged_original_seqs_R2.fastq.gz


# don't comparess yet since we will jus thave to uncompress n the next step
python3 split_seqs.py -i tagged_original_seqs_R1.fastq -o split_seqs_R1.fastq
python3 split_seqs.py -i tagged_original_seqs_R2.fastq -o split_seqs_R2.fastq

gzip -f split_seqs_R1.fastq
gzip -f split_seqs_R2.fastq

/opt/software/bwa-0.7.17/bwa mem -t 4 /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa split_seqs_R1.fastq.gz \
        | /opt/software/samtools/bin/samtools view -bht /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa - \
        | /opt/software/samtools/bin/samtools sort -o final_out_R1.bam -;
        /opt/software/samtools/bin/samtools index final_out_R1.bam

/opt/software/bwa-0.7.17/bwa mem -t 4 /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa split_seqs_R2.fastq.gz \
        | /opt/software/samtools/bin/samtools view -bht /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa - \
        | /opt/software/samtools/bin/samtools sort -o final_out_R2.bam -;
        /opt/software/samtools/bin/samtools index final_out_R2.bam

/opt/software/bwa-0.7.17/bwa mem -t 4 /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa untagged_original_seqs_R1.fastq.gz \
        | /opt/software/samtools/bin/samtools view -bht /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa - \
        | /opt/software/samtools/bin/samtools sort -o tagless_out_R1.bam -;
        /opt/software/samtools/bin/samtools index tagless_out_R1.bam

/opt/software/bwa-0.7.17/bwa mem -t 4 /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa untagged_original_seqs_R2.fastq.gz \
        | /opt/software/samtools/bin/samtools view -bht /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa - \
        | /opt/software/samtools/bin/samtools sort -o tagless_out_R2.bam -;
        /opt/software/samtools/bin/samtools index tagless_out_R2.bam

# add read groups:
java -jar /opt/software/picard/picard.jar AddOrReplaceReadGroups \
    I=final_out_R1.bam \
    O=final_out_R1.grp.bam \
    RGID=R1_tagged \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=NULL \
    RGSM=MM30

java -jar /opt/software/picard/picard.jar AddOrReplaceReadGroups \
    I=final_out_R2.bam \
    O=final_out_R2.grp.bam \
    RGID=R2_tagged \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=NULL \
    RGSM=MM30

java -jar /opt/software/picard/picard.jar AddOrReplaceReadGroups \
    I=tagless_out_R1.bam \
    O=tagless_out_R1.grp.bam \
    RGID=R1_untagged \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=NULL \
    RGSM=MM30

java -jar /opt/software/picard/picard.jar AddOrReplaceReadGroups \
    I=tagless_out_R2.bam \
    O=tagless_out_R2.grp.bam \
    RGID=R2_untagged \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=NULL \
    RGSM=MM30

# merge the BAMs and index:
/opt/software/samtools/bin/samtools merge xyz.bam final_out_R1.grp.bam final_out_R2.grp.bam tagless_out_R1.grp.bam tagless_out_R2.grp.bam
/opt/software/samtools/bin/samtools index xyz.bam

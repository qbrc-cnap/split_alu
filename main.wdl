
workflow TaggedReadSplitWorkflow{

    Array[File] r1_files
    Array[File] r2_files
    String barcode
    File bwa_fa
    File bwa_amb
    File bwa_ann
    File bwa_bwt
    File bwa_fai
    File bwa_pac
    File bwa_sa
    File bwa_dict

    String git_repo_url
    String git_commit_hash

    Array[Pair[File, File]] fastq_pairs = zip(r1_files, r2_files)

    scatter(item in fastq_pairs){

        call trim_reads {
            input:
                r1 = item.left,
                r2 = item.right
        } 

        call split_and_align_reads {
            input:
                r1 = trim_reads.trimmed_r1,
                r2 = trim_reads.trimmed_r2,
                barcode = barcode,
                bwa_fa = bwa_fa,
                bwa_amb = bwa_amb,
                bwa_ann = bwa_ann,
                bwa_bwt = bwa_bwt,
                bwa_fai = bwa_fai,
                bwa_pac = bwa_pac,
                bwa_sa = bwa_sa,
                bwa_dict = bwa_dict
        }   
    }

    output {
        Array[File] trimmed_r1 = trim_reads.trimmed_r1
        Array[File] trimmed_r2 = trim_reads.trimmed_r2
        Array[File] bams = split_and_align_reads.bam
        Array[File] bais = split_and_align_reads.bai
        Array[File] beds = split_and_align_reads.bed
        Array[File] flagstats = split_and_align_reads.flagstat
    }

    meta {
        workflow_title : "Split tagged reads and align with BWA"
        workflow_short_description : "Split tagged reads and align with BWA"
        workflow_long_description : "Split tagged reads and align with BWA"
    }

}

task trim_reads {

    File r1
    File r2

    # Extract the samplename from the fastq filename
    String sample_name = basename(r1, "_R1.fastq.gz")

    Int disk_size = 200


    command {
        java -jar /opt/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
            -trimlog ${sample_name}.trim.log \
            -summary ${sample_name}.trim_summary.log \
            ${r1} ${r2} \
            -baseout ${sample_name}.trimmed.fastq.gz \
            ILLUMINACLIP:/opt/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:8:true
    }

    output {
        File trimmed_r1 = "${sample_name}.trimmed_1P.fastq.gz"
        File trimmed_r2 = "${sample_name}.trimmed_2P.fastq.gz"
    }

    runtime {
        docker: "docker.io/blawney/alu_split:v1"
        cpu: 4
        memory: "12 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

task split_and_align_reads {

    File r1
    File r2
    String barcode
    File bwa_fa
    File bwa_amb
    File bwa_ann
    File bwa_bwt
    File bwa_fai
    File bwa_pac
    File bwa_sa
    File bwa_dict

    Int disk_size = 200

    String sample_name = basename(r1, ".trimmed_1P.fastq.gz")

    command {

        cutadapt \
            -e 0.12 \
            -a ${barcode} \
            --rc \
            --action=none \
            --untrimmed-output /dev/null \
            -o output_R1.fastq.gz \
            ${r1} 2>&1>${sample_name}.R1_barcode_report.txt

        cutadapt \
            -e 0.12 \
            -a ${barcode} \
            --rc \
            --action=none \
            --untrimmed-output /dev/null \
            -o output_R2.fastq.gz \
            ${r2} 2>&1>${sample_name}.R2_barcode_report.txt


        # get the union set of reads that have tags in them and then split out the original FASTQ:
        python3 /opt/software/extract_readnames.py \
            -r1 output_R1.fastq.gz \
            -r2 output_R2.fastq.gz \
            -o ${sample_name}.readnames.list

        /opt/software/bbmap/filterbyname.sh in=${r1} names=${sample_name}.readnames.list include=t out=${sample_name}.tagged_original_seqs_R1.fastq
        /opt/software/bbmap/filterbyname.sh in=${r2} names=${sample_name}.readnames.list include=t out=${sample_name}.tagged_original_seqs_R2.fastq
        /opt/software/bbmap/filterbyname.sh in=${r1} names=${sample_name}.readnames.list include=f out=${sample_name}.untagged_original_seqs_R1.fastq.gz
        /opt/software/bbmap/filterbyname.sh in=${r2} names=${sample_name}.readnames.list include=f out=${sample_name}.untagged_original_seqs_R2.fastq.gz

        # don't compress yet since we will just have to uncompress in the next step
        python3 /opt/software/split_seqs.py -i ${sample_name}.tagged_original_seqs_R1.fastq -o ${sample_name}.split_seqs_R1.fastq -b ${barcode}
        python3 /opt/software/split_seqs.py -i ${sample_name}.tagged_original_seqs_R2.fastq -o ${sample_name}.split_seqs_R2.fastq -b ${barcode}

        gzip -f ${sample_name}.split_seqs_R1.fastq
        gzip -f ${sample_name}.split_seqs_R2.fastq

        /opt/software/bwa-0.7.17/bwa mem -t 4 /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${sample_name}.split_seqs_R1.fastq.gz \
                | /opt/software/samtools/bin/samtools view -bht /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa - \
                | /opt/software/samtools/bin/samtools sort -o ${sample_name}.tagged.R1.bam -;
                /opt/software/samtools/bin/samtools index ${sample_name}.tagged.R1.bam

        /opt/software/bwa-0.7.17/bwa mem -t 4 /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${sample_name}.split_seqs_R2.fastq.gz \
                | /opt/software/samtools/bin/samtools view -bht /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa - \
                | /opt/software/samtools/bin/samtools sort -o ${sample_name}.tagged.R2.bam -;
                /opt/software/samtools/bin/samtools index ${sample_name}.tagged.R2.bam

        /opt/software/bwa-0.7.17/bwa mem -t 4 /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${sample_name}.untagged_original_seqs_R1.fastq.gz \
                | /opt/software/samtools/bin/samtools view -bht /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa - \
                | /opt/software/samtools/bin/samtools sort -o ${sample_name}.untagged.R1.bam -;
                /opt/software/samtools/bin/samtools index ${sample_name}.untagged.R1.bam

        /opt/software/bwa-0.7.17/bwa mem -t 4 /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${sample_name}.untagged_original_seqs_R2.fastq.gz \
                | /opt/software/samtools/bin/samtools view -bht /workspace/full_bwa_index/Homo_sapiens.GRCh38.dna.primary_assembly.fa - \
                | /opt/software/samtools/bin/samtools sort -o ${sample_name}.untagged.R2.bam -;
                /opt/software/samtools/bin/samtools index ${sample_name}.untagged.R2.bam

        # add read groups:
        java -jar /opt/software/picard/picard.jar AddOrReplaceReadGroups \
            I=${sample_name}.tagged.R1.bam \
            O=${sample_name}.tagged.R1.w_grps.bam \
            RGID=R1_tagged \
            RGLB=lib1 \
            RGPL=ILLUMINA \
            RGPU=NULL \
            RGSM=${sample_name}

        java -jar /opt/software/picard/picard.jar AddOrReplaceReadGroups \
            I=${sample_name}.tagged.R2.bam \
            O=${sample_name}.tagged.R2.w_grps.bam \
            RGID=R2_tagged \
            RGLB=lib1 \
            RGPL=ILLUMINA \
            RGPU=NULL \
            RGSM=${sample_name}

        java -jar /opt/software/picard/picard.jar AddOrReplaceReadGroups \
            I=${sample_name}.untagged.R1.bam \
            O=${sample_name}.untagged.R1.w_grps.bam \
            RGID=R1_untagged \
            RGLB=lib1 \
            RGPL=ILLUMINA \
            RGPU=NULL \
            RGSM=${sample_name}

        java -jar /opt/software/picard/picard.jar AddOrReplaceReadGroups \
            I=${sample_name}.untagged.R2.bam \
            O=${sample_name}.untagged.R2.w_grps.bam \
            RGID=R2_untagged \
            RGLB=lib1 \
            RGPL=ILLUMINA \
            RGPU=NULL \
            RGSM=${sample_name}


        # merge the BAMs and index:
        /opt/software/samtools/bin/samtools merge ${sample_name}.merged.bam \
            ${sample_name}.tagged.R1.w_grps.bam \
            ${sample_name}.tagged.R2.w_grps.bam \
            ${sample_name}.untagged.R1.w_grps.bam \
            ${sample_name}.untagged.R2.w_grps.bam \
        /opt/software/samtools/bin/samtools index ${sample_name}.merged.bam

        # Run coverage and flagstat
        /opt/software/bedtools2/bin/bedtools genomecov -ibam ${sample_name}.merged.bam -bga > ${sample_name}.bed
        /opt/software/samtools-1.10/samtools flagstat ${sample_name}.merged.bam > ${sample_name}.flagstat
    }

    output {
        File bam = "${sample_name}.merged.bam"
        File bai = "${sample_name}.merged.bam.bai"
        File bed = "${sample_name}.bed"
        File flagstat = "${sample_name}.flagstat"
    }

    runtime {
        docker: "docker.io/blawney/alu_split:v1"
        cpu: 4
        memory: "24 G"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 0
    }
}

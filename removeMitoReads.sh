#!/bin/bash

if [ $# -ne 5 ]; then
        echo -e "\nusage: checkInsertSize.sh <mitochondria.fasta> <forward_reads.fastq> <reverse_reads.fastq> <forward_noMito.fastq> <reverse_noMito.fastq>\n"
        exit
else
        if [ -e "$1.amb" -a -e "$1.ann" -a -e "$1.bwt" -a -e "$1.pac" -a -e "$1.sa" ]; then
                echo -e "\nMitochondria already indexed by bwa\nMapping now!\n"
                date
                bwa mem -t 64 $1 $2 $3 > montagem.sam
                echo -e "\nConverting SAM to BAM\n"
                date
                samtools view -b -S montagem.sam > montagem.bam
                echo -e "\nGenerating BAM file with only unmapped reads...\n"
                date
                samtools view -u -f 12 -F 256 montagem.bam > both_unmapped.bam
        else
                echo -e "\nIndexing genome...\n"
                date
                bwa index $1
                echo -e "\nMitochondria indexed by bwa\nMapping now!\n"
                date
                bwa mem -t 64 $1 $2 $3 > montagem.sam
                echo -e "\nConverting SAM to BAM\n"
                date
                samtools view -b -S montagem.sam > montagem.bam
                echo -e "\nGenerating BAM file with only unmapped reads...\n"
                date
                samtools view -u -f 12 -F 256 montagem.bam > both_unmapped.bam
        fi
        echo -e "\nSorting BAM file...\n"
        date
        samtools sort -n both_unmapped.bam > montagem_sorted.bam
        if [ -e montagem_sorted.bam ]; then
                echo -e "\nRemoving original SAM/BAM file to save space\n"
                date
                rm montagem.bam
                rm montagem.sam
                rm both_unmapped.bam
                echo -e "\nConvertng BAM to fastq files...\n"
                date
                bamToFastq -i montagem_sorted.bam -fq $4 -fq2 $5
        else
                exit
        fi
        echo -e "\nRemoving sorted BAM file\n"
        date
        rm montagem_sorted.bam
        if [[ $2 =~ \.gz$ ]]; then
		count1=$(cat $4 | wc -l)
		count2=$(zcat $2 | wc -l)
		percentage=$(echo "$count1" "$count2" | awk '{div = $1 / $2 * 100; final = 100 - div; print final}')
		echo -e "\n$percentage% of your reads are mitochondrial"
        else
		count1=$(cat $4 | wc -l)
                count2=$(cat $2 | wc -l)
                percentage=$(echo "$count1" "$count2" | awk '{div = $1 / $2 * 100; final = 100 - div; print final}')
                echo -e "\n$percentage% of your reads are mitochondrial"
	fi
	echo -e "\nDONE!\n"
        date
fi

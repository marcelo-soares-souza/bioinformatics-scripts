#!/bin/bash

if [ $# -ne 4 ]; then
        echo -e "\nusage: checkInsertSize.sh <genome.fasta> <forward_reads.fastq> <reverse_reads.fastq> <metrics_output>\n"
        exit
else
	if [ -e "$1.amb" -a -e "$1.ann" -a -e "$1.bwt" -a -e "$1.pac" -a -e "$1.sa" ]; then
        	echo -e "\nGenome already indexed by bwa\nMapping now!\n"
		bwa mem -t 64 $1 $2 $3 > montagem.sam
	else
		echo -e "\nIndexing genome...\n"
		bwa index $1
		echo -e "\nGenome indexed by bwa\nMapping now!\n"
		bwa mem -t 64 $1 $2 $3 > montagem.sam
	fi
	echo -e "\nSorting SAM file...\n"
        samtools sort montagem.sam > montagem_sorted.sam
        if [ -e montagem_sorted.sam ]; then
		echo -e "\nRemoving original SAM file to save space\n"
                rm montagem.sam
		echo -e "\nCollecting metrics using picard-tools...\n"
                picard-tools CollectInsertSizeMetrics I=montagem_sorted.sam R=$1 O=$4 H=histogram.out
        else
                exit
        fi
	echo -e "\nRemoving sorted SAM file\n"
	rm montagem_sorted.sam
	echo -e "\nDONE!\nOpen $4 to check alignment metrics\n"
fi

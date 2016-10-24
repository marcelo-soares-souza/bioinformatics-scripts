#!/bin/bash

for FILE in `ls *.vcf`; do
  FILENAME="${FILE%.*}"
  echo -e "\n${FILENAME}\n"
  grep -v "^#" $FILE | cut -d$'\t' -f1,2,4,5 | sort -u -t$'\t' -k1,1 -k 2,2 | uniq > ${FILENAME}.csv
done

#!/bin/bash

for FILE in `ls *.vcf`; do
  FILENAME="${FILE%.*}"
  echo -e "\n${FILENAME}\n"
  grep -v "^#" $FILE | cut -d$'\t' -f1,2,4,5 | sort | uniq > ${FILENAME}.csv
done

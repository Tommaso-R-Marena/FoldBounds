#!/bin/bash

# Assumes that hpfold has been build and exists in the dist directory

benchmarkInputFile="benchmark.txt"
resultFilePath="last_results.txt"
backupResultFilePath="old_last_result.txt"
hpfold="./../dist/hpfold"

rm -f $backupResultFilePath
if [[ -f $resultFilePath ]]; then
   mv $resultFilePath $backupResultFilePath
fi

if [[ $# -gt 0 && -f $1 ]]; then
	echo $1
   benchmarkInputFile=$1
fi


counter=0

while read -r chain || [[ -n $chain ]]
do
  (time echo "$chain" | eval $hpfold) >> $resultFilePath 2>&1
  printf '\n\n' >> $resultFilePath
  printf '=%.0s' {1..100} >> $resultFilePath
  printf '\n\n\n' >> $resultFilePath
  ((counter++))
  echo "${counter} chain(s) processed."
done < "$benchmarkInputFile"


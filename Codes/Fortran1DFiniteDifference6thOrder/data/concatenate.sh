#!/usr/bin/bash

touch results.txt

for i in 128 512 2048 8192 32768
do
  for tt in 128 512 2048 8192 32768
  do
    mv results.txt temp.txt
    cat  errortime${i}and${tt}.dat temp.txt > results.txt
  done
done


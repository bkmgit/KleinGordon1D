#!/usr/bin/bash

for i in 128 512 2048 8192 32768
do
  for tt in 128 512 2048 8192 32768
  do
    rm KgSemiImp1d${i}and${tt}.f90 
    rm makefileKabuki${i}and${tt}
    rm subscriptKabuki${i}and${tt}
    rm kg${i}and${tt}
    rm FFT${i}and${tt}.*
  done
done

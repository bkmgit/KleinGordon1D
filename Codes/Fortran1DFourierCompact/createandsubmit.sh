#!/usr/bin/bash

for i in 128 512 2048 8192 32768
do
  for tt in 128 512 2048 8192 32768
  do
    cp KgSemiImp1d.f90 KgSemiImp1d${i}and${tt}.f90 
    sed -i "88s/Nx=4096/Nx=${i}/"  KgSemiImp1d${i}and${tt}.f90
    sed -i "89s/Nt=512/Nt=${tt}/"  KgSemiImp1d${i}and${tt}.f90
    cp makefileKabuki makefileKabuki${i}and${tt}
    sed -i "s/kg/kg${i}and${tt}/" makefileKabuki${i}and${tt}
    sed -i "s/kg.o/kg${i}and${tt}.o/" makefileKabuki${i}and${tt}
    sed -i "s/KgSemiImp1d.f90/KgSemiImp1d${i}and${tt}.f90/" makefileKabuki${i}and${tt}
    sed -i "s/-o kg.o/-o kg${i}and${tt}.o/" makefileKabuki${i}and${tt}
    sed -i "s/-o kg \$(LIBS)/-o kg${i}and${tt} \$(LIBS)/" makefileKabuki${i}and${tt}
    cp subscriptKabuki subscriptKabuki${i}and${tt}
    sed -i "2s/#PBS -N FFT/#PBS -N FFT${i}and${tt}/" subscriptKabuki${i}and${tt}
    sed -i "12s/kg/kg${i}and${tt}/" subscriptKabuki${i}and${tt}
    make -f makefileKabuki${i}and${tt}
    qsub subscriptKabuki${i}and${tt}
  done
done

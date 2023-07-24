#!/bin/bash

read -p "Enter last cycle: " Cycle
NextCycle=Cycle$(($Cycle+1))

mkdir -pv $NextCycle/{1.lammps,2.QM,3.fit}

cp -vi Cycle$Cycle/3.fit/lammps.txt $NextCycle/1.lammps/in.minimization

cp -vi Cycle$Cycle/1.lammps/{submit.sh,bodies.txt} $NextCycle/1.lammps/

cp -vi Cycle$Cycle/1.lammps/mix.restart.30000 $NextCycle/1.lammps/mix.restart.30000.anterior

cd $NextCycle/1.lammps

sbatch submit.sh



#/bin/bash

read -p "Enter current cycle: " Cycle

cd ~/NANCY/augusto/Cycle$Cycle/1.lammps

DME_mapeo_augusto_7_23

cp -v ~/NANCY/augusto/Cycle$Cycle/1.lammps/*.com ~/NANCY/augusto/Cycle$Cycle/2.QM/

cp -vi ~/NANCY/augusto/submit_dilithium.sh ~/NANCY/augusto/Cycle$Cycle/2.QM/

cd ~/NANCY/augusto/Cycle$Cycle/2.QM

sbatch submit_dilithium.sh



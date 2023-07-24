#/bin/bash

read -p "Enter current cycle: " Cycle

PrevCycle=$(($Cycle-1))

cd ~/NANCY/augusto/Cycle$Cycle/2.QM

Check=$(grep -cl "Normal termination" *.log | wc -l)

Files=(Archivo*.com)
NumberFiles=${#Files[@]}

if [ $Check != $NumberFiles ]; then
  echo "Some cycles ended wrong"
  exit 1
else
  go.picky.recovernrg b3lyp
fi

cd ~/NANCY/augusto/Cycle$PrevCycle/3.fit

echo $PWD

/home/amusetti/NANCY/augusto/generar_pickyfit_inp_sig.x < $PrevCycle

cd -

cp -v ./geo.list.dat ../3.fit

cp -v ~/NANCY/augusto/Cycle$PrevCycle/3.fit/{cycle$PrevCycle.QM.dat,pickyfit.DME.cycle$Cycle.inp} ~/NANCY/augusto/Cycle$Cycle/3.fit

sed ~/NANCY/augusto/Cycle$Cycle/3.fit/pickyfit.DME.cycle$Cycle.inp -i -re '4,19 s/$/  1  0/g; 20,21 s/$/  2  0/g; 85 s/[0-9]+/'"$Cycle"'/'

# cd ~/NANCY/augusto/Cycle$Cycle/3.fit
# go.pickyfit

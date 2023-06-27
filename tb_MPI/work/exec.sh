#!/bin/bash

#cd src

#Para compilação e execução otimizada com intel
#for i in $(seq 1 4)
#do
#mpiifort -complex-limited-range -O3 -w -xHost -mp1 -ipo -fast=2 *.f90 -o tbhpc.x
#mv tbhpc.x ../.
#cd ../work
#(/usr/bin/time -f "%e" mpirun -np $i ./../tbhpc.x) &> ../out_cores$i.dat
#done

# Ciclo para a grade 24 x 15
for i in $(seq 1 3)
do
(/usr/bin/time -f "%e" mpirun -np $i ./../tbhpc.x)  &>> ../out_g1_np$i.dat
done
mv autovalores.dat autovalores_g1.dat

# Ciclo para a grade 20 x 48
cp input_malha2.dat input.dat
for i in $(seq 1 4)
do
(/usr/bin/time -f "%e" mpirun -np $i ./../tbhpc.x)  &>> ../out_g2_np$i.dat
done
mv autovalores.dat autovalores_g2.dat

# Ciclo para a grade 15 x 192

cp input_malha3.dat input.dat
for i in $(seq 1 4)
do
(/usr/bin/time -f "%e" mpirun -np $i ./../tbhpc.x)  &>> ../out_g3_np$i.dat
done
mv autovalores.dat autovalores_g3.dat

echo "Trabalho Concluido!"

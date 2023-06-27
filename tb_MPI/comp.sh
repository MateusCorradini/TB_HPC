#!/bin/bash

cd src

#Para compilação e execução otimizada com intel
mpiifort -complex-limited-range -O3 -w -xHost -mp1 -ipo -fast=2 *.f90 -o tbhpc.x
mv tbhpc.x ../.

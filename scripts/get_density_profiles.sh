#!/bin/bash

T=20000
dT=5000
step=1000

path=dens
mkdir ${path}
for t_start in $(seq 0 ${step} $(( T-dT ))); do
    t_end=$(( t_start+dT ))
    folder=${path}/dens_${t_end}

    mkdir ${folder}

    echo 3 3 | gmx density -f cal_dec_tip4p.xtc -s cal_dec_tip4p.tpr -o ${folder}/dens.xvg -dens number -sl 200 -d Z -b ${t_start} -e ${t_end} -center
    echo 0 | gmx trjconv -f cal_dec_tip4p.xtc -s cal_dec_tip4p.tpr -o ${folder}/dump_${t_end}.gro -dump ${t_end}
done

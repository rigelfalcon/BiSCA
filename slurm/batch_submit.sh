#!/bin/bash

# para for taskmng Ntask
Ntask=111;#110;#35;#1;#105;
npertask=16;#1772#3;#105;#1;
slurmfile="./slurm/step3_tri_3model_multiple_node.sl"

# loop for submit task 
for ((itask=1; itask<=$Ntask; itask++)); do
# for ((itask=1; itask<=2; itask++)); do
    sed -i "s/Ntask=[0-9]\+/Ntask=$Ntask/" "$slurmfile"
    sed -i "s/itask=[0-9]\+/itask=$itask/" "$slurmfile"
    sed -i "s/npertask=[0-9]\+/npertask=$npertask/" "$slurmfile"
    sbatch "$slurmfile"
done

#!/bin/bash

source /etc/profile.d/modules.sh

module load papi

#text
srun_ps ./gccg text input_files/tjunc.dat tjunk1;       
srun_ps ./gccg text input_files/tjunc.dat tjunk2;
srun_ps ./gccg text input_files/tjunc.dat tjunk3;

srun_ps ./gccg text input_files/drall.dat drall1;
srun_ps ./gccg text input_files/drall.dat drall2;
srun_ps ./gccg text input_files/drall.dat drall3;

srun_ps ./gccg text input_files/pent.dat pent1;
srun_ps ./gccg text input_files/pent.dat pent2;
srun_ps ./gccg text input_files/pent.dat pent3;

srun_ps ./gccg text input_files/cojack.dat cojack1;
srun_ps ./gccg text input_files/cojack.dat cojack2;
srun_ps ./gccg text input_files/cojack.dat cojack3;

exit;



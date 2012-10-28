#!/bin/bash
#SBATCH -o /home/hpc/h039u/h039uag/workspace/pos2012ws/A1/code/output.%j.%N.out
#SBATCH -D /home/hpc/h039u/h039uag/workspace/pos2012ws/A1/code
#SBATCH -J mrcsrv
#SBATCH --get-user-env
#SBATCH --clusters=mpp1
#SBATCH --mail-type=end
#SBATCH --mail-user=marco.seravalli@tum.de
#SBATCH --export=NONE
#SBATCH --time=00:09:00

source /etc/profile.d/modules.sh

module load papi

#text
./gccg text input_files/tjunc.dat tjunk1
./gccg text input_files/tjunc.dat tjunk2
./gccg text input_files/tjunc.dat tjunk3

./gccg text input_files/drall.dat drall1
./gccg text input_files/drall.dat drall2
./gccg text input_files/drall.dat drall3

./gccg text input_files/pent.dat pent1
./gccg text input_files/pent.dat pent2
./gccg text input_files/pent.dat pent3

./gccg text input_files/cojack.dat cojack1
./gccg text input_files/cojack.dat cojack2
./gccg text input_files/cojack.dat cojack3

#binary
./gccg bin input_files/tjunc.bin tjunk1
./gccg bin input_files/tjunc.bin tjunk2
./gccg bin input_files/tjunc.bin tjunk3

./gccg bin input_files/drall.bin drall1
./gccg bin input_files/drall.bin drall2
./gccg bin input_files/drall.bin drall3

./gccg bin input_files/pent.bin pent1
./gccg bin input_files/pent.bin pent2
./gccg bin input_files/pent.bin pent3

./gccg bin input_files/cojack.bin cojack1
./gccg bin input_files/cojack.bin cojack2
./gccg bin input_files/cojack.bin cojack3


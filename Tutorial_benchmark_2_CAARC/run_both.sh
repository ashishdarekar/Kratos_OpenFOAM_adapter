#!/bin/sh
set -e -u

# To Run OpenFOAM
blockMesh
decomposePar
# change name of 0 folder, make 0 dummy folder
mpirun -np 2 snappyHexMesh -overwrite -parallelrwrite -parallel
reconstructParMesh -constant
# delete 0 folder and reane the o folder to old

. "${WM_PROJECT_DIR}/bin/tools/RunFunctions"
solver=$(getApplication)
${solver} > logopenfoam &

# To Run KRATOS
export KRATOS_SOURCE=~/Documents/MS/Kratos/Kratos
alias startkratos="setupkratosenv /home/ashish/Documents/MS/Kratos/Kratos"
#startkratos = /home/ashish/Documents/MS/Kratos/Kratos
startkratos

runkratosmpi MainKratosCoSim.py 1 > logkratos &



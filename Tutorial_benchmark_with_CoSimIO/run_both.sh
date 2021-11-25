#!/bin/sh
set -e -u

# To Run OpenFOAM
blockMesh
. "${WM_PROJECT_DIR}/bin/tools/RunFunctions"
solver=$(getApplication)
${solver} > logopenfoam &

# To Run KRATOS
export KRATOS_SOURCE=~/Documents/MS/Kratos/Kratos
alias startkratos="setupkratosenv /home/ashish/Documents/MS/Kratos/Kratos"
#startkratos = /home/ashish/Documents/MS/Kratos/Kratos
startkratos

runkratosmpi MainKratosCoSim.py 1 > logkratos &



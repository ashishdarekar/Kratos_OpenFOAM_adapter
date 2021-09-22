#!/bin/sh
set -e -u

# To Run OpenFOAM
blockMesh
. "${WM_PROJECT_DIR}/bin/tools/RunFunctions"
solver=$(getApplication)
${solver} > logopenfoam &

# To Run KRATOS
#startkratos = /home/ashish/Documents/MS/Kratos/Kratos
#$startkratos

python3 MainKratosCoSim.py > logkratos &



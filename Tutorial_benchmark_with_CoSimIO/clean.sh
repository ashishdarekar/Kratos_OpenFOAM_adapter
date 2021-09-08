#!/bin/sh
set -e -u

echo "--- Cleaning CoSimulation result files in $(pwd)"
rm -rfv .CoSimIOFileComm_Openfoam_Adapter_Openfoam_Kratos_Wrapper/
rm -rfv vtk_output/
rm -rfv vtk_output_coupling/

echo "--- Cleaning OpenFOAM result files in $(pwd)"
if [ -n "${WM_PROJECT:-}" ] || error "No OpenFOAM environment is active."; then
    # shellcheck disable=SC1090 # This is an OpenFOAM file which we don't need to check
    . "${WM_PROJECT_DIR}/bin/tools/CleanFunctions"
    cleanCase
    rm -rfv 0/uniform/functionObjects/functionObjectProperties
fi

echo "--- Cleaning Post-processing files in $(pwd)"
rm disp.txt
rm forces.txt

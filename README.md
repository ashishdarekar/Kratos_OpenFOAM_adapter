# OpenFoam_Kratos_adapter
This is for my Master Thesis

# How to run a first tutorial: Cavity
# Requirements:
1. OpenFOAM 7 (https://openfoam.org/download/7-ubuntu/)
2. CoSimIO (https://github.com/KratosMultiphysics/CoSimIO)

# Steps to be followed:
1. Building KratosOpenfomaAdapterFunctionObject:
    1. Go to the folder .../KratosOpenfoamAdapterfunctionObject (https://github.com/ashishdarekar/OpenFoam_Kratos_adapter/tree/main/KratosOpenfoamAdapterFunctionObject)
    2. Go to SuperUser: *sudo -E bash*
    3. *wclean*
    4. *wmake*

2. To compile export and import data files
    1. Go to the folder .../Tutorial_cavity (https://github.com/ashishdarekar/OpenFoam_Kratos_adapter/tree/main/Tutorial_cavity)
    2. Go to SuperUser: *sudo -E bash*
    3. "g++ trial_exporting_data_using_CoSimIO.cpp -o export"
    4. "g++ trial_importing_data_using_CoSimIO.cpp -o import"

2. To run the case:
    1. Go to the folder .../Tutorial_cavity (https://github.com/ashishdarekar/OpenFoam_Kratos_adapter/tree/main/Tutorial_cavity)
    2. Go to SuperUser: *sudo -E bash*
    3. *blockMesh*
    4. *icoFoam*
    5. Simultaneously open another terminal to export and import data in OpenFoam using CoSimIO:
        1. *./export*
        2. *./import*

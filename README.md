# OpenFoam_Kratos_adapter
This is for my Master Thesis

# How to run a tutorial on Cavity
# requirement:
1. OpenFOAM 7
2. CoSimIO

# Steps to be followed
1. Building KratosOpenfomaAdapterFunctionObject:
    Go to the folder .../KratosOpenfoamAdapterfunctionObject
    Go to SuperUser by typing sudo -E bash
    wclean
    wmake
2. To run a case:
    Go to the folder .../Tutorial_cavity
    Go to SuperUser by typing sudo -E bash
    blockMesh
    icoFoam

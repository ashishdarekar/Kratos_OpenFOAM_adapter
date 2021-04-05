# OpenFoam_Kratos_adapter
For more information about this Master Thesis: [Abstract of this Master thesis](https://github.com/ashishdarekar/OpenFoam_Kratos_adapter/blob/main/Abstract_of_Master_Thesis_ashish_darekar.pdf)

# Requirements:
1. [OpenFOAM-7](https://openfoam.org/download/7-ubuntu/)
2. [CoSimIO](https://github.com/KratosMultiphysics/CoSimIO)

# How to run the first tutorial: Cavity
1. **To Build the KratosOpenfomaAdapterFunctionObject:**
    1. Add Path of the **CoSimulationAdapter.H** ($HOME/OpenFoam_Kratos_adapter/Adapter) into your include path.
    2. Go to the folder .../KratosOpenfoamAdapterfunctionObject (https://github.com/ashishdarekar/OpenFoam_Kratos_adapter/tree/main/KratosOpenfoamAdapterFunctionObject)
    3. Go to SuperUser: *sudo -E bash*
    4. *wclean*
    5. *wmake*

3. **To run the first tutorial case:**
    1. Go to the folder .../Tutorial_cavity (https://github.com/ashishdarekar/OpenFoam_Kratos_adapter/tree/main/Tutorial_cavity)
    2. In one terminal run the OpneFOAM commands:
        1. Go to SuperUser: *sudo -E bash*
        2. *blockMesh*
        3. *pimpleFoam*
    3. Simultaneously open another terminal to export and import data in OpenFoam using CoSimIO:
        1. Go to SuperUser: *sudo -E bash*
        2. *g++ trial_export_import_data_using_CoSimIO.cpp -o export_import*

4. **To run the FSI-Benchmarking case:**
    1. Go to the folder .../Tutorial_benchmark (https://github.com/ashishdarekar/OpenFoam_Kratos_adapter/tree/main/Tutorial_benchmark)
    2. In one terminal run the OpneFOAM commands:
        1. Go to SuperUser: *sudo -E bash*
        2. *blockMesh*
        3. *pimpleFoam*
    3. Simultaneously open another terminal to export and import data in OpenFoam using CoSimIO:
        1. Go to SuperUser: *sudo -E bash*
        2. *g++ trial_export_import_data_using_CoSimIO.cpp -o export_import*

**Note:** To compile the export and import data files, Dont forget to add Path of the **co_sim_io.hpp** ($HOME/CoSimIO) into your include path.

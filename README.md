# OpenFoam_Kratos_adapter
For more information about this Master Thesis: [Abstract of this Master thesis](https://github.com/ashishdarekar/OpenFoam_Kratos_adapter/blob/main/Abstract_of_Master_Thesis_ashish_darekar.pdf)

# Requirements:
1. [OpenFOAM-7](https://openfoam.org/download/7-ubuntu/)
2. [KRATOS](https://github.com/KratosMultiphysics/Kratos)
3. [Scripts](https://github.com/philbucher/bash_scripts) - Bash scripts for consistent & convenient use of Kratos
4. [CoSimIO](https://github.com/KratosMultiphysics/CoSimIO)

# How to run the first tutorial: Cavity
1. **To Build the KratosOpenfomaAdapterFunctionObject:**
    1. Go to the folder > KratosOpenfoamAdapterfunctionObject
    2. Make linker related modifications in the file KratosOpenfoamAdapterfunctionObject/Make/options
        1. ```
            -I$<Path_to_directory_CoSimIO>/co_sim_io \
            -I$<Path_to_directory_CoSimIO> \
            -L$<Path_to_directory_CoSimIO>/bin
           ```

        2. To MPI related error *MPI_ROOT* should direct to *openmpi/include* in your system:
        e.g. using following commad (For future use its better to add this in *bashrc*)
           ```
           export MPI_ROOT=/usr/lib/x86_64-linux-gnu/openmpi/include
           ```

    3.  Give following Commands to compile:
        ```
        wclean
        wmake
        ```

2. **To run the FSI-Benchmarking case without MPI:**
    1. Go to the folder > Tutorial_benchmark_with_CoSimIO
    2. In one terminal run the OpneFOAM commands:
        ```
        blockMesh
        pimpleFoam
        ```
    3. Simultaneously, open another terminal to run KRATOS:
        ```
        startkratos
        runkratos MainKratosCoSim.py
        ```
    4. To see the results in ParaView:
        1.  ```
            paraFoam -case .&
            ```
        2. Load the KRATOS results from the folder > vtk_structure_output
    5. Displacement of the tip can be seen using:
        ```
        python plot_displacement.py
        ```
    6. To clean all generated files during simulation:
        ```
        ./clean.sh
        ```
2. **To run the FSI-Benchmarking case with MPI:**
    1. Go to the folder > Tutorial_benchmark_with_CoSimIO
    2. In one terminal run the OpneFOAM commands:
        ```
        blockMesh
        decomposePar
        mpirun -np <number_of_processes> pimpleFoam -parallel
        ```
    3. Simultaneously, open another terminal to run KRATOS:
        ```
        startkratos
        runkratosmpi MainKratosCoSim.py <number_of_processes>
        ```
    4. To see the results in ParaView:
        1.  ```
            paraFoam -case .&
            ```
        2. Load the KRATOS results from the folder > vtk_structure_output
    5. Displacement of the tip can be seen using:
        ```
        python plot_displacement.py
        ```
    6. To clean all generated files during simulation:
        ```
        ./clean.sh
        ```

**Note:** To compile the export and import data files, Dont forget to add Path of the **co_sim_io.hpp** ($HOME/CoSimIO) into your include path.

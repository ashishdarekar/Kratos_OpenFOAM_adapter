# Kratos-OpenFOAM adapter
For more information about this Master Thesis: [Abstract of this Master thesis](https://github.com/ashishdarekar/OpenFoam_Kratos_adapter/blob/main/Abstract_of_Master_Thesis_ashish_darekar.pdf)

# Requirements:
1. [OpenFOAM-7](https://openfoam.org/download/7-ubuntu/) - Fluid Solver
2. [KRATOS](https://github.com/KratosMultiphysics/Kratos) - Structural Solver and "CoSimulation" Coupling tool
3. [CoSimIO](https://github.com/KratosMultiphysics/CoSimIO) - Communication tool between the solvers
4. [Scripts](https://github.com/philbucher/bash_scripts) - Bash scripts for consistent & convenient use of Kratos

# Configuration of the Adapter function object:
1. Go to the folder > KratosOpenfoamAdapterFunctionObject
2. Make linker related modifications in the file KratosOpenfoamAdapterFunctionObject/Make/options
    1. ```
        -I$<Path_to_directory_CoSimIO>/co_sim_io \
        -I$<Path_to_directory_CoSimIO> \
        -L$<Path_to_directory_CoSimIO>/bin
        ```

    2. If MPI related error occurs, *MPI_ROOT* should direct to *openmpi/include* in your system:
    e.g. using following commad (For future use its better to add this in *bashrc*)
        ```
        export MPI_ROOT=/usr/lib/x86_64-linux-gnu/openmpi/include
        ```

3.  Give following Commands to compile:
    ```
    wclean
    wmake
    ```
4. It will generate *"libKratosOpenfoamAdapterFunctionObjectFunctionObjects.so"* shared library file in the same folder. One need to add this file in controlDict while using this function object.

# Tutorials:
1. **Case 1 - Flow-induced vibration of a flexible beam or flap (without MPI):**
    1. Go to the folder > Tutorial_case_1_flap
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
        2. Load the KRATOS results from the folder > vtk_output_structure
    5. Displacement of the tip can be seen using:
        ```
        python3 plot_displacement.py
        ```
    6. To clean all generated files during simulation:
        ```
        ./clean.sh
        ```

2. **Case 1 - Flow-induced vibration of a flexible beam or flap (with MPI):**
    1. Go to the folder > Tutorial_case_1_flap
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
            reconstructPar
            paraFoam -case .&
            ```
        2. Load the KRATOS results from the folder > vtk_output_structure
    5. Displacement of the tip can be seen using:
        ```
        python3 plot_displacement.py
        ```
    6. To clean all generated files during simulation:
        ```
        ./clean.sh
        ```

3. **Case 2 - FSI of the CAARC building simulation (with MPI):**
    1. Go to the folder > Tutorial_case_2_CAARC
    2. In one terminal run the OpneFOAM commands:
        ```
        blockMesh
        decomposePar
        mv 0/ 0_org/
        mkdir 0
        mpirun -np <number_of_processes> snappyHexMesh -overwrite -parallel
        reconstructParMesh -constant
        rm -r 0/ processor0/ processor1/ processor2/ ..... processor<number_of_processes-1>
        mv 0_org/ 0/
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
            reconstructPar
            paraFoam -case .&
            ```
        2. Load the KRATOS results from the folder > vtk_output_structure
    5. Displacement of the points specified in the *ProjectParametesCSD.json* can be seen using:
        ```
        python3 plot_displacement.py
        ```
    6. To clean all generated files during simulation:
        ```
        ./clean.sh
        ```

**Note:** To compile the export and import data files, Dont forget to add Path of the **co_sim_io.hpp** ($HOME/CoSimIO) into your include path.

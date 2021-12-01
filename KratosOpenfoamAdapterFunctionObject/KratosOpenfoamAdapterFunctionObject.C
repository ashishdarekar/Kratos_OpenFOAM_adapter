/*-----------------------------------------------------------------------*\

Master-Thesis Work
Ashish Darekar
Sourcefile for the KratosOpenfoamAdpterFunctionObject.H

\*-----------------------------------------------------------------------*/

#include "KratosOpenfoamAdapterFunctionObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(KratosOpenfoamAdapterFunctionObject, 0);
    addToRunTimeSelectionTable(functionObject, KratosOpenfoamAdapterFunctionObject, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::KratosOpenfoamAdapterFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
fvMeshFunctionObject(name, runTime, dict), runTime_(runTime), dict_(dict)//, CoSimulationAdapter_()
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::~KratosOpenfoamAdapterFunctionObject()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::read(const dictionary& dict)
{
    debugInfo("CoSimulation Adapter's function object : read()", debugLevel);

    if(readConfig(dict))
    {
        connectKratos();

        exportMeshToCosim();

        resizeDataVectors();
    }

    return true;
}

bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::execute()
{
    debugInfo( runTime_.timeName()  + " : CoSimulation Adapter's function object : execution()", debugLevel);

    exportDataToKratos();

    importDataFromKratos();

    return true;
}

bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::end()
{
    debugInfo("CoSimulation Adapter's function object : end()", debugLevel);

    // DisConnection between OpenFOAM and Kratos-CoSimulation using CoSimIO
    disconnectKratos();

    return true;
}

bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::write()
{
    debugInfo("CoSimulation Adapter's function object : write()", debugLevel);

    return true;
}



// ****************************************** Auxillary functions Read config *********************************************** //
bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::readConfig(const dictionary& dict)
{
    // Reading configuration parameters from ControlDict related to function object
    fvMeshFunctionObject::read(dict);

    debugLevel = dict.lookupOrDefault<word>("debugLevel", "info");

    dim = dict.lookupOrDefault<int>("dim", 3);
    debugInfo( "Dimension of a problem is: "  + std::to_string(dim), debugLevel);

    thick = dict.lookupOrDefault<double>("thick", 1.0);
    debugInfo( "Thickness of a domain is: "  + std::to_string(thick), debugLevel);

    if(runTime_.controlDict().lookupOrDefault("adjustTimeStep", false))
    {
        debugInfo( "Cannot support adjustable time step", debugLevel);

        //EXIT THE SIMULATION
        connectKratos();
        end();
    }

    // Check the solver type and determine it if needed
    solverType_ = dict.lookupOrDefault<word>("solvertype", "none");

    if (solverType_.compare("compressible") == 0 || solverType_.compare("incompressible") == 0)
    {
        debugInfo( "Known solver type: " + solverType_ , debugLevel);
    }
    else if (solverType_.compare("none") == 0)
    {
        debugInfo("Determining the solver type", debugLevel);
        solverType_ = determineSolverType();
    }
    else
    {
        debugInfo("Unknown solver type. Determining the solver type", debugLevel);
        solverType_ = determineSolverType();
    }

    // Every interface is a subdictionary of "interfaces"
    const dictionary * interfaceSubdictPtr = dict.subDictPtr("interfaces");

    debugInfo("Coupling Interfaces Reading: Start", debugLevel);

    if(!interfaceSubdictPtr)
    {
        debugInfo("No Coupling Interfaces found", debugLevel);

        //EXIT THE SIMULATION
        connectKratos();
        end();
    }
    else
    {
        for(const entry& interfaceSubdictEntry : *interfaceSubdictPtr)
        {
            if(interfaceSubdictEntry.isDict())
            {
                dictionary interfaceSubdict = interfaceSubdictEntry.dict();

                struct InterfaceData interfacedata;

                interfacedata.nameOfInterface = interfaceSubdict.lookupType<word>("name");
                debugInfo( "Name of the coupling interface is = " + interfacedata.nameOfInterface , debugLevel);

                wordList patches = interfaceSubdict.lookupType<wordList>("patches");
                for(auto patch : patches)
                {
                    interfacedata.patchNames.push_back(patch);
                }

                wordList importData = interfaceSubdict.lookupType<wordList>("importData");
                for(auto rData : importData)
                {
                    interfacedata.importData.push_back(rData);
                }

                wordList importDataIdentifier = interfaceSubdict.lookupType<wordList>("importDataIdentifier");
                for(auto rData : importDataIdentifier)
                {
                    interfacedata.importDataIdentifier.push_back(rData);
                }

                wordList exportData = interfaceSubdict.lookupType<wordList>("exportData");
                for(auto wData : exportData)
                {
                    interfacedata.exportData.push_back(wData);
                }

                wordList exportDataIdentifier = interfaceSubdict.lookupType<wordList>("exportDataIdentifier");
                for(auto rData : exportDataIdentifier)
                {
                    interfacedata.exportDataIdentifier.push_back(rData);
                }

                //Add this interface in the Array of all interfaces
                interfaces_.push_back(interfacedata);

                num_interfaces_++;
            }
        }
    }

    // Loop to save the PatchIds associated with each interface for later use
    for(std::size_t j = 0; j < num_interfaces_; j++)
    {
        // For every patch that participates in the coupling interface. We are keeping one patch for one interface
        for (std::size_t i = 0; i < interfaces_.at(j).patchNames.size(); i++)
        {
            // Get the patchID
            int patchID = mesh_.boundaryMesh().findPatchID((interfaces_.at(j).patchNames).at(i));

            // Throw an error if the patch was not found
            if (patchID == -1)
            {
                debugInfo( "ERROR: Patch " + (interfaces_.at(j).patchNames).at(i) + " does not exist.", debugLevel );
            }

            // Add the patch in the list
            interfaces_.at(j).patchIDs.push_back(patchID);
        }

    }

    debugInfo("Coupling Interfaces Reading: Done" , debugLevel);
    debugInfo("Number of coupling interfaces found: " + std::to_string(num_interfaces_) , debugLevel);

    return true;

}

void Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::connectKratos()
{
    // Connection between openFOAM and Kratos-CoSimulation using CoSimIO (ONLY ONE TIME for multiple interfaces)
    CoSimIO::Info settings;
    settings.Set("my_name", "Openfoam_Adapter");
    settings.Set("connect_to", "Openfoam_Kratos_Wrapper");
    settings.Set("echo_level", 0);
    settings.Set("version", "1.25");
    CoSimIO::Info connect_info;

    if(TotalNumOfProcesses == 1)
    {
        connect_info = CoSimIO::Connect(settings);
    }
    else{
        connect_info = CoSimIO::ConnectMPI(settings, MPI_COMM_WORLD);
    }

    COSIMIO_CHECK_EQUAL(connect_info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Connected);
    connection_name = connect_info.Get<std::string>("connection_name");
}

void Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::disconnectKratos()
{
    CoSimIO::Info connect_info;
    CoSimIO::Info disconnect_settings;
    disconnect_settings.Set("connection_name", connection_name);
    connect_info = CoSimIO::Disconnect(disconnect_settings);
    COSIMIO_CHECK_EQUAL(connect_info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Disconnected);

}

void Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::exportMeshToCosim()
{
    debugInfo( "Exporting All InterfaceMeshes to KRATOS using CoSimIO::ModelPart : Start", debugLevel);

    for(std::size_t j = 0; j < num_interfaces_; j++)
    {
        debugInfo( "[OF]Reading mesh data from OF for coupling interface : " + interfaces_.at(j).nameOfInterface, debugLevel );

        // Arrays to exchange the number of nodes with all ranks
        scalarListList sendDataNumNodeIndex(Pstream::nProcs());
        scalarListList recvDataNumNodeIndex(Pstream::nProcs());

        // Count the Nodes and Elements for all the patches in that interface
        for (std::size_t i = 0; i < interfaces_.at(j).patchIDs.size(); i++)
        {
            interfaces_.at(j).numNodes += mesh_.boundaryMesh()[interfaces_.at(j).patchIDs.at(i)].localPoints().size();
            interfaces_.at(j).numElements += mesh_.boundary()[interfaces_.at(j).patchIDs[i]].size();
        }

        // Save the Number of data to send to other Ranks
        for(int i=0; i< TotalNumOfProcesses ; i++)
        {
            sendDataNumNodeIndex[i].setSize(1);
            sendDataNumNodeIndex[i][0] = interfaces_.at(j).numNodes ;
        }

        // Make CoSimIO::ModelPart and push in the array of model_part_interfaces
        model_part_interfaces_.push_back(CoSimIO::make_unique<CoSimIO::ModelPart>(interfaces_.at(j).nameOfInterface));

        //MPI Exchange to know start Index for Node formation
        Pstream::exchange<scalarList, scalar>(sendDataNumNodeIndex, recvDataNumNodeIndex);

        int nodeIndex = 1 ; //Default value of Serial run
        int elemIndex = 1;

        for(int i = 0 ; i < MyRank; i++)
        {
            nodeIndex += recvDataNumNodeIndex[i][0];
        }
        int globalNodeIndexBegin = nodeIndex; //To presrve its value(use is later)

        // Accessing the coordinates of all nodes in the Inteface and collecting nodal and elemental data
        for(std::size_t i = 0; i < interfaces_.at(j).patchIDs.size(); i++)
        {
            label patchIndex1 = mesh_.boundaryMesh().findPatchID((interfaces_.at(j).patchNames).at(i)); //Patch Id of Interface Flap
            const UList<label> &bfaceCells1 = mesh_.boundaryMesh()[patchIndex1].faceCells();
            label patchIndex2 = 0; //Check the default value, What should I provide?

            // Travel through all the nodes in that MPI Rank
            forAll(bfaceCells1, bfacei1)
            {
                const label& faceID1 = mesh_.boundaryMesh()[interfaces_.at(j).patchIDs[i]].start() + bfacei1;

                Element new_element = Element(elemIndex++);//make new element and then increase the elemental index counter

                forAll(mesh_.faces()[faceID1], nodei1)
                {
                    const label& nodeID1 = mesh_.faces()[faceID1][nodei1]; //for OpenFOAM
                    auto pointX = mesh_.points()[nodeID1];
                    std::vector<bool> is_common_node(TotalNumOfProcesses , 0); //Total number of max neighbours  = NumOfProcesses-1

                    int result = compare_nodes(pointX); // return nodeIndex(starting from 1) if node is already present and (-1) if node is not present

                    if(result == (-1)) // For new node
                    {
                        //--------------------------------- FORALL Neighbouring Boundaries -----------------------------------//
                        forAll(mesh_.boundaryMesh() , ipatch)
                        {
                            word BCtype = mesh_.boundaryMesh().types()[ipatch];
                            if(BCtype == "processor")
                            {
                                patchIndex2 = ipatch; //patchIndex 2 can be replace with ipatch
                                const processorPolyPatch& pp = refCast<const processorPolyPatch>( mesh_.boundaryMesh()[patchIndex2] );
                                int nighbor_proc_id = pp.neighbProcNo();

                                const UList<label> &bfaceCells2 = mesh_.boundaryMesh()[patchIndex2].faceCells();

                                //While creation of node check if it has to be created in the way of ghost or local, CHECK it in ALL neighbour or sharing boundaries
                                forAll(bfaceCells2, bfacei2) //compare with all the nodes in the Processor common patch
                                {
                                    const label& faceID2 = mesh_.boundaryMesh()[patchIndex2].start() + bfacei2;
                                    forAll(mesh_.faces()[faceID2], nodei2)
                                    {
                                        const label& nodeID2 = mesh_.faces()[faceID2][nodei2]; //for OpenFOAM
                                        auto pointY = mesh_.points()[nodeID2];

                                        if(is_same_points(pointX,pointY) == 1 )//once this is found Go out from the For loop and create the node
                                        {
                                            is_common_node.at(nighbor_proc_id) = 1; //need to create this node as a gghost node
                                        }
                                    }
                                }
                            }
                        }

                        bool new_node_created = 0; //ON/OFF
                        for(int i = 0; i< TotalNumOfProcesses ; i++)
                        {
                            if(is_common_node.at(i) == 1)//Any one
                            {
                                if(new_node_created == 0)
                                {
                                    Node new_node = Node(pointX, nodeIndex, nodeIndex, i);//cosimIndex kept same as local index (by default), later need to change
                                    new_node_created = 1;
                                    interfaces_.at(j).Interface_nodes.push_back(new_node);
                                    new_element.addNodeIndexInList(nodeIndex);// connectivity to make that element
                                    new_element.addNodesInList(new_node);// Add new node in the element
                                    nodeIndex++;
                                    array_of_nodes.push_back(pointX);// Push new node in the list to compare
                                }
                                else{
                                    interfaces_.at(j).Interface_nodes.at(nodeIndex-2).setCommonWithRank(i); //With which Id it is common
                                }
                            }
                        }
                        is_common_node.clear(); //No need anyways it will be deleted once it goes out of scope

                        if(new_node_created == 0) //if Still node is not created then only create it normally
                        {
                            Node new_node = Node(pointX, nodeIndex, nodeIndex, MyRank);//cosimIndex kept same as local index (by default), later need to change
                            interfaces_.at(j).Interface_nodes.push_back(new_node);
                            new_element.addNodeIndexInList(nodeIndex);// connectivity to make that element
                            new_element.addNodesInList(new_node);// Add new node in the element
                            nodeIndex++;
                            array_of_nodes.push_back(pointX);// Push new node in the list to compare
                        }

                    }
                    else{
                        //No need to produce new Node
                        new_element.addNodeIndexInList(result);// connectivity to make that element
                        new_element.addNodesInList(interfaces_.at(j).Interface_nodes.at(result-1));// Add old node in the List of nodes for that element, result-1 as Vector index starts from 0
                    }
                }
                //Once Element is set push it to the List of Elements
                interfaces_.at(j).Interface_elements.push_back(new_element);
            }

            debugInfo( "[OF]Total number of Nodes in coupling interface "  + interfaces_.at(j).nameOfInterface + " are " +  std::to_string(interfaces_.at(j).Interface_nodes.size()), debugLevel);
            debugInfo( "[OF]Total number of Elements in coupling interface "  + interfaces_.at(j).nameOfInterface + " are " +  std::to_string(interfaces_.at(j).Interface_elements.size()), debugLevel);

        }

        //Save the info about neighbour rank and corresponding number of nodes common with
        int num_common_nodes = 0;
        for(int k = 0 ; k<TotalNumOfProcesses ;k++)
        {
            for(auto& nodei: interfaces_.at(j).Interface_nodes)
            {
                if(nodei.getCommonWithRank().at(k) == 1)
                {
                    num_common_nodes++;
                }
            }
            interfaces_.at(j).neighbour_ids_comm_num_of_nodes.push_back(k);
            interfaces_.at(j).neighbour_ids_comm_num_of_nodes.push_back(num_common_nodes);
            num_common_nodes = 0 ; //Reset to 0 , to count for next
        }

        // ------------------Common Nodal Data exchange with all ranks --------------------------------//
        scalarListList sendNodalData(Pstream::nProcs());
        scalarListList recvNodalData(Pstream::nProcs());

        // Decide the size of a send array
        for(int i = 0; i < TotalNumOfProcesses ; i++)
        {
            sendNodalData[i].setSize(4 * interfaces_.at(j).neighbour_ids_comm_num_of_nodes.at( 2*i + 1 )); //4 for each node, 1 for node index and 3 points coordinates
        }

        // Filiing of the send array
        std::vector<int>counter(TotalNumOfProcesses , 0);
        for (auto& nodei: interfaces_.at(j).Interface_nodes)
        {
            for(std::size_t i = 0 ;  i < nodei.getCommonWithRank().size() ; i++)
            {
                if(nodei.getCommonWithRank()[i] == 1)
                {
                    sendNodalData[i][(counter.at(i))++] = nodei.getLocalNodeIndex(); //Provide Local Id only
                    vector nodePosition  = nodei.getNodePosition();
                    sendNodalData[i][(counter.at(i))++] = nodePosition[0];
                    sendNodalData[i][(counter.at(i))++] = nodePosition[1];
                    sendNodalData[i][(counter.at(i))++] = nodePosition[2];
                }
            }
        }
        counter.clear();

        // MPI_Exchange
        Pstream::exchange<scalarList, scalar>(sendNodalData, recvNodalData);

        // Update the nodeIndexes for common nodes (before making CoSim nodes)
        for (auto& nodei: interfaces_.at(j).Interface_nodes)
        {
            for(std::size_t i = MyRank ;  i < nodei.getCommonWithRank().size() ; i++)
            {
                if(nodei.getCommonWithRank()[i] == 1)
                {
                    vector nodePosition  = nodei.getNodePosition();

                    for(int j = 0; j<recvNodalData[i].size()/4 ; j++)
                    {
                        vector trial_node(recvNodalData[i][j*4+1], recvNodalData[i][j*4+2], recvNodalData[i][j*4+3]);

                        if(is_same_points(nodePosition,trial_node))
                        {
                            nodei.setNodeIndexForCoSim(recvNodalData[i][j*4+0]); //Take node Id of that
                        }
                    }
                }
            }
        }

        // Make CoSim nodes
        for(auto& nodei : interfaces_.at(j).Interface_nodes)
        {
            if(nodei.getLocalNodeIndex() != nodei.getNodeIndexForCoSim() && (nodei.getHighestCommRank()<TotalNumOfProcesses)) //If both the ids are not same then they are ghost
            {
                vector nodePosition  = nodei.getNodePosition();
                model_part_interfaces_.at(j)->CreateNewGhostNode( nodei.getNodeIndexForCoSim(), nodePosition[0], nodePosition[1], nodePosition[2], nodei.getHighestCommRank());
            }
            else //Local nodes
            {
                vector nodePosition  = nodei.getNodePosition();
                model_part_interfaces_.at(j)->CreateNewNode( nodei.getNodeIndexForCoSim(), nodePosition[0], nodePosition[1], nodePosition[2]);
            }
        }
        debugInfo( "[COSIM]Total number of Nodes formed in coupling interface " + interfaces_.at(j).nameOfInterface +  " (local, Ghost, total) = (" + std::to_string(model_part_interfaces_.at(j)->NumberOfLocalNodes()) + " , " + std::to_string(model_part_interfaces_.at(j)->NumberOfGhostNodes()) +  " , " + std::to_string(model_part_interfaces_.at(j)->NumberOfNodes()) + " )" , debugLevel);

        // Make Cosim elements
        std::vector<CoSimIO::IdType> connectivity;
        int quad_count = 0;
        int pyramid_count= 0;
        int prism_count = 0;

        for(auto& elementi : interfaces_.at(j).Interface_elements)
        {
            for(auto& elementalNodeIndexi : elementi.getElementalNodes())
            {
                connectivity.push_back(interfaces_.at(j).Interface_nodes.at(elementalNodeIndexi.getLocalNodeIndex()-globalNodeIndexBegin).getNodeIndexForCoSim());
            }

            // To select the Element from the available list of elements in the CoSimIO
            switch (connectivity.size())
            {
            case 4:
                model_part_interfaces_.at(j)->CreateNewElement( elementi.getLocalElementIndex(), CoSimIO::ElementType::Quadrilateral2D4, connectivity );
                quad_count++;
                break;
            case 5: // When not written, CoSimIO error
                model_part_interfaces_.at(j)->CreateNewElement( elementi.getLocalElementIndex(), CoSimIO::ElementType::Pyramid3D5, connectivity );
                pyramid_count++;
                break;
            case 6: // When not written, these elements got skipped, No CoSimIO error
                model_part_interfaces_.at(j)->CreateNewElement( elementi.getLocalElementIndex(), CoSimIO::ElementType::Prism3D6, connectivity );
                prism_count++;
                break;
            default: // Add more cases for more elements
                break;
            }

            connectivity.clear();
        }
        debugInfo( "[COSIM]Total number of Quad Elements formed in coupling interface " + interfaces_.at(j).nameOfInterface + " = " + std::to_string(quad_count), debugLevel);
        debugInfo( "[COSIM]Total number of Pyramid Elements formed in coupling interface " + interfaces_.at(j).nameOfInterface + " = " + std::to_string(pyramid_count), debugLevel);
        debugInfo( "[COSIM]Total number of Prism Elements formed in coupling interface " + interfaces_.at(j).nameOfInterface + " = " + std::to_string(prism_count), debugLevel);

        debugInfo( "[COSIM]Total number of Elements formed in coupling interface " + interfaces_.at(j).nameOfInterface + " = " + std::to_string(model_part_interfaces_.at(j)->NumberOfElements()), debugLevel);

        // Export InterfaceMesh/ModelPart to CoSimulation using CoSimIO
        info.Clear();
        info.Set("identifier", interfaces_.at(j).nameOfInterface);
        info.Set("connection_name", connection_name);
        auto export_info = CoSimIO::ExportMesh(info, *model_part_interfaces_.at(j));
        debugInfo( "Finished Exporting interface Mesh " +  interfaces_.at(j).nameOfInterface + " to Kratos as a ModelPart "  , debugLevel);

        // Clear all the entries of temporary arrays for this interface, make it available for nex interface
        array_of_nodes.clear();
        sendNodalData.clear();
        recvNodalData.clear();
        sendDataNumNodeIndex.clear();
        recvDataNumNodeIndex.clear();
    }
    debugInfo( "Exporting All InterfaceMeshes to KRATOS using CoSimIO::ModelPart : End", debugLevel);

    Foam::wordList Objectnames_ = mesh_.names(); // calling names() method
    forAll(Objectnames_,i){
        std::cout << Objectnames_[i] << ", ";}
}

void Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::resizeDataVectors()
{
    // Resizing the Data Vectors for Import Export operations with CoSimIO
    // Import Data from CoSimulation (Displacement/Delta) present only on faceNodes/Nodes. No need to resize it
    for(std::size_t i=0; i < num_interfaces_; i++)
    {
        // Export Data to CoSimulation (force/Stress) present only on faceCenters/Elements
        for(std::size_t j=0; j< interfaces_.at(i).exportData.size(); j++)
        {
            std::string dataName = interfaces_.at(i).exportData.at(j);

            if(dataName.find("Force") == 0 || dataName.find("Stress") == 0) //If "force" or "stress" string is found it will return 0
            {
                interfaces_.at(i).data_to_send.resize((interfaces_.at(i).numElements) * dim);
            }
        }

    }
}

void Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::exportDataToKratos()
{
    for(std::size_t i=0; i < num_interfaces_; i++)
    {
        debugInfo( "Force calculation started for coupling interface : " + interfaces_.at(i).nameOfInterface , debugLevel);

        // For "Write Data" variables which need to send to CoSimulation
        for(std::size_t j=0; j< interfaces_.at(i).exportData.size(); j++)
        {
            std::string dataName = interfaces_.at(i).exportData.at(j);

            if(dataName.find("Force") == 0 )
            {
                calculateForces(i);
            }

            // Conversion Utilities for converting Elemental Load values to Nodal values (Variable is load/Force/Reaction now)
            conversionElementalToNodalValues(i);
        }

        // Export this force array to CoSimulation //Elemental Force Data
        CoSimIO::Info connect_info;
        connect_info.Clear();
        connect_info.Set("identifier", interfaces_.at(i).exportDataIdentifier[0]);
        connect_info.Set("connection_name", connection_name);
        connect_info = CoSimIO::ExportData(connect_info, interfaces_.at(i).data_to_send);

        //debugInfo("Data has been exported from OpenFOAM to CoSimulation (for coupling interface name = " + interfaces_.at(i).nameOfInterface + ") , Force values with array size = " + std::to_string(interfaces_.at(i).data_to_send.size()), debugLevel);
        debugInfo("Data has been exported from OpenFOAM to CoSimulation (for coupling interface name = " + interfaces_.at(i).nameOfInterface + ")", debugLevel);
    }
}

void Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::importDataFromKratos()
{
    for (std::size_t i = 0; i < interfaces_.size(); i++)
    {
        // Import the displacement array from the CoSimulation
        CoSimIO::Info connect_info;
        connect_info.Clear();
        connect_info.Set("identifier", interfaces_.at(i).importDataIdentifier[0]);
        connect_info.Set("connection_name", connection_name);
        connect_info = CoSimIO::ImportData(connect_info, interfaces_.at(i).data_to_recv);
        COSIMIO_CHECK_EQUAL(interfaces_.at(i).data_to_recv.size() , (interfaces_.at(i).numNodes) * dim ); //Check size of Receive data = Number of nodes*dim Is it require??

        //debugInfo( "Data has been imported from CoSimulation to OpenFOAM: (for coupling interface name = " + interfaces_.at(i).nameOfInterface + ") , Disp values with array size = " + std::to_string(interfaces_.at(i).data_to_recv.size()) , debugLevel);
        debugInfo( "Data has been imported from CoSimulation to OpenFOAM: (for coupling interface name = " + interfaces_.at(i).nameOfInterface + ")" , debugLevel);

        // Get the displacement on the patch(for every patch in the interface) and assign it those values received from CoSimulation
        for (std::size_t j = 0; j < interfaces_.at(i).patchNames.size(); j++)
        {
            Foam::pointVectorField* point_disp;

            point_disp = const_cast<pointVectorField*>( &mesh_.lookupObject<pointVectorField>("pointDisplacement") );
            label patchIndex = mesh_.boundaryMesh().findPatchID((interfaces_.at(i).patchNames).at(j));
            fixedValuePointPatchVectorField& pointDisplacementFluidPatch = refCast<fixedValuePointPatchVectorField>(point_disp->boundaryFieldRef()[patchIndex]);

            int iterator = 0;
            forAll(point_disp->boundaryFieldRef()[patchIndex] ,k)
            {
                pointDisplacementFluidPatch[k][0] = interfaces_.at(i).data_to_recv[iterator++];
                pointDisplacementFluidPatch[k][1] = interfaces_.at(i).data_to_recv[iterator++];
                if (dim ==3)
                    pointDisplacementFluidPatch[k][2] = interfaces_.at(i).data_to_recv[iterator++];
            }
        }

        debugInfo( "Displacement replacement ended for coupling interface : " + interfaces_.at(i).nameOfInterface, debugLevel);
    }
}



// *********************************************** Some Auxillary Functions **************************************************//
// Calculate the Solver Type - according to the pressure dimension
std::string Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::determineSolverType()
{
    dimensionSet pressureDimensionsCompressible(1, -1, -2, 0, 0, 0, 0);
    dimensionSet pressureDimensionsIncompressible(0, 2, -2, 0, 0, 0, 0);

    if (mesh_.foundObject<volScalarField>("p"))
    {
        volScalarField p_ = mesh_.lookupObject<volScalarField>("p");

        if (p_.dimensions() == pressureDimensionsCompressible)
        {
            solverType_ = "compressible";
            debugInfo("Solver Type : Compressible ", debugLevel);
        }
        else if (p_.dimensions() == pressureDimensionsIncompressible)
        {
            solverType_ = "incompressible";
            debugInfo("Solver Type : InCompressible ", debugLevel);
        }
    }

    if(solverType_  == "unknown")
    {
        debugInfo("Solver Type: Neither Compressible nor Incompresible", debugLevel);
    }

    return solverType_;
}

// To compare the Foam::Vector. Check whether it is really required?
bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::is_same_points(Foam::vector& pointX, Foam::vector& pointY)
{
    if(pointX[0] == pointY[0] && pointX[1] == pointY[1] && pointX[2] == pointY[2])
        return true;
    else
        return false;
}

// To Compare the new node with all previous nodes before creating the new CoSimIO::Node
int Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::compare_nodes(Foam::vector& pointX)
{
    int answer = 0;

    if(array_of_nodes.size()==0)
    {
        answer = -1;
    }
    else
    {
        for(Foam::vector& nodei : array_of_nodes)
        {
            if(is_same_points(pointX , nodei)) //current node == previous all nodes
            {
                auto itr = std::find( array_of_nodes.begin(), array_of_nodes.end(), nodei ) ;
                answer = std::distance(array_of_nodes.begin(), itr) + 1 ; // nodeindex starts from 1 in CoSimIO
                break; // Once repeated node is found immediate break it. No need to check it later
            }
            else
            {
                answer = -1; // Compare with all the available nodes and then return (-1) if pointX is not found. Then new node is created in the main()
            }
        }
    }

    return answer;
}

// ************************** Force / load calculation ************************************//
// For viscous Force
Foam::tmp<Foam::volSymmTensorField> Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::devRhoReff() const
{
    //For turbulent flows
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (mesh_.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const cmpTurbModel & turb
        (
            mesh_.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName)
        );

        return turb.devRhoReff();

    }
    else if (mesh_.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const incompressible::turbulenceModel& turb
        (
            mesh_.lookupObject<icoTurbModel>(icoTurbModel::propertiesName)
        );

        return rho()*turb.devReff();
    }
    else
    {
        // For laminar flows get the velocity
        const volVectorField & U
        (
            mesh_.lookupObject<volVectorField>("U")
        );

        return -mu()*dev(twoSymm(fvc::grad(U)));
    }
}

// Finding correct rho
Foam::tmp<Foam::volScalarField> Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::rho() const
{
    // If volScalarField exists, read it from registry (for compressible cases)
    // interFoam is incompressible but has volScalarField rho
    if (mesh_.foundObject<volScalarField>("rho"))
    {
        return mesh_.lookupObject<volScalarField>("rho");
    }
    else if (solverType_.compare("incompressible") == 0)
    {
        const dictionary& FSIDict = dict_.subOrEmptyDict("Parameters");

        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(FSIDict.lookup("rho"))
            )
        );
    }
    else
    {
        FatalErrorInFunction << "Exiting the simulation: correct rho not found" << exit(FatalError);

        return volScalarField::null();
    }
}

// Finding correct mu
Foam::tmp<Foam::volScalarField> Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::mu() const
{
    if (solverType_.compare("incompressible") == 0)
    {
        typedef immiscibleIncompressibleTwoPhaseMixture iitpMixture;
        if (mesh_.foundObject<iitpMixture>("mixture"))
        {
            const iitpMixture& mixture
            (
                mesh_.lookupObject<iitpMixture>("mixture")
            );

            return mixture.mu();
        }
        else
        {
            const dictionary& FSIDict = dict_.subOrEmptyDict("Parameters");

            dimensionedScalar nu(FSIDict.lookup("nu"));

            return tmp<volScalarField>
            (
                new volScalarField
                (
                    nu*rho()
                )
            );
        }

    }
    else if (solverType_.compare("compressible") == 0)
    {
        return mesh_.lookupObject<volScalarField>("thermo:mu");
    }
    else
    {
        FatalErrorInFunction << "Exiting the simulation: correct mu not found" << exit(FatalError);

        return volScalarField::null();
    }
}

// Normal vectors multiplied by face area
Foam::vectorField Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::getFaceVectors(const unsigned int patchID) const
{
    // Normal vectors multiplied by face area
    return mesh_.boundary()[patchID].Sf();
}

// Calculate Total Force
bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::calculateForces(std::size_t interface_index)
{
    // Get different force fields from OpenFOAM, Refering Force Function Object (OpenFOAM User guide)
    // 1.Stress tensor boundary field
    tmp<volSymmTensorField> tdevRhoReff(devRhoReff());
    const volSymmTensorField::Boundary& devRhoReffb
    (
        tdevRhoReff().boundaryField()
    );

    // 2.Density boundary field
    tmp<volScalarField> trho(rho());
    const volScalarField::Boundary& rhob = trho().boundaryField();

    // 3.Pressure boundary field
    tmp<volScalarField> tp = mesh_.lookupObject<volScalarField>("p");
    const volScalarField::Boundary& pb
    (
        tp().boundaryField()
    );

    // For every boundary patch of the interface
    for(std::size_t j=0; j< interfaces_.at(interface_index).patchIDs.size(); j++)
    {
        int patchID = interfaces_.at(interface_index).patchIDs.at(j);

        const auto& surface = getFaceVectors(patchID);

        // Pressure forces
        if(solverType_.compare("incompressible") == 0)
        {
            Force_->boundaryFieldRef()[patchID] = (surface/thick) * pb[patchID] * rhob[patchID];
        }
        else if(solverType_.compare("compressible") == 0)
        {
            Force_->boundaryFieldRef()[patchID] = (surface/thick) * pb[patchID];
        }
        else
        {
            FatalErrorInFunction << "Forces calculation does only support compressible or incompressible solver type." << exit(FatalError);
        }

        // Viscous forces
        Force_->boundaryFieldRef()[patchID] += (surface/thick) & devRhoReffb[patchID];

        // Writing this forces into Buffer to export to CoSimulation
        int bufferIndex = 0;
        // For every cell in the patch
        forAll(Force_->boundaryField()[patchID], i)
        {
            // Copy the force into the buffer
            // x-dimension
            interfaces_.at(interface_index).data_to_send[bufferIndex++] = Force_->boundaryField()[patchID][i].x();

            // y-dimension
            interfaces_.at(interface_index).data_to_send[bufferIndex++] = Force_->boundaryField()[patchID][i].y();

            if(dim == 3)
            {
                // z-dimension
                interfaces_.at(interface_index).data_to_send[bufferIndex++] = Force_->boundaryField()[patchID][i].z();
            }
        }

    }

    return true;
}

void Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::conversionElementalToNodalValues(std::size_t interface_index)
{
    int bufferIndex = 0;
    int number_of_nodes_in_element = 0;

    // Travel all Elements and Distribute the forces on the nodes
    for(auto& elementi : interfaces_.at(interface_index).Interface_elements)
    {
        number_of_nodes_in_element = elementi.getNumberOfNodes();

        for(auto& elementalNodeIndexi : elementi.getElementalNodes())
        {
            elementalNodeIndexi.getLoadValues()[0] +=  ( (interfaces_.at(interface_index).data_to_send[bufferIndex++]) / double(number_of_nodes_in_element) );
            elementalNodeIndexi.getLoadValues()[1] +=  ( (interfaces_.at(interface_index).data_to_send[bufferIndex++]) / double(number_of_nodes_in_element) );
            elementalNodeIndexi.getLoadValues()[2] +=  ( (interfaces_.at(interface_index).data_to_send[bufferIndex++]) / double(number_of_nodes_in_element) );
        }
    }

    // Clear all the entried of data_to_send array
    interfaces_.at(interface_index).data_to_send.clear();

    // Resize the data_to_send to keep nodal data
    interfaces_.at(interface_index).data_to_send.resize( (interfaces_.at(interface_index).numNodes) * dim);

    // Fill the data_to_send to keep nodal data
    // Data_send will be in the OF order
    bufferIndex = 0;
    for(auto& nodei : interfaces_.at(interface_index).Interface_nodes)
    {
        interfaces_.at(interface_index).data_to_send[bufferIndex++] = nodei.getLoadValues()[0];
        interfaces_.at(interface_index).data_to_send[bufferIndex++] = nodei.getLoadValues()[1];
        interfaces_.at(interface_index).data_to_send[bufferIndex++] = nodei.getLoadValues()[2];
    }
}

// *********************************************** Utilities **************************************************//
void Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::debugInfo(const std::string message, const std::string level )
{
    if(level.compare("debug") == 0)
    {
        Pout << "[AdapterInfo] " << message.c_str() << endl;
    }
    else // by default "info"
    {
        Info << "[AdapterInfo] " << message.c_str() << nl;
    }
}
}//namespace Foam

// ************************************************************************* //
/*--------------------------------------------------------------------------------------------------------------*\

Master-Thesis Work
Ashish Darekar
Sourcefile for the KratosOpenfoamAdpterFunctionObject.H

\*--------------------------------------------------------------------------------------------------------------*/

#include "KratosOpenfoamAdapterFunctionObject.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(KratosOpenfoamAdapterFunctionObject, 0);
    addToRunTimeSelectionTable(functionObject, KratosOpenfoamAdapterFunctionObject, dictionary);
}


// ****************************************** Constructors ******************************************************//
Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::KratosOpenfoamAdapterFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
fvMeshFunctionObject(name, runTime, dict), runTime_(runTime), dict_(dict)
{
    read(dict);
}

// ******************************************* Destructor ******************************************************//
Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::~KratosOpenfoamAdapterFunctionObject()
{
}


// ********************************* FunctionObject's Member Functions ***************************************** //
bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::read(const dictionary& dict)
{
    debugInfo("CoSimulation Adapter's function object : read()", debugLevel);

    if(readConfig(dict))
    {
        // Coneect to CoSimulation
        connectKratos();

        // Export Mesh related data to CoSimulation in the form of ModelPart (only once)
        exportMeshToCosim();

        // Resize the data vectors required for communication
        resizeDataVectors();
    }

    return true;
}

bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::execute()
{
    debugInfo( runTime_.timeName()  + " : CoSimulation Adapter's function object : execution()", debugLevel);

    // Export Load data to CoSimulation
    exportDataToKratos();

    // Import Displacement data from CoSimulation
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
    return true;
}


// ********************************* Auxiliar functions for configuration ****************************************//
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
        //Connect and Dicsconnect
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
        //Connect and Dicsconnect
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

                wordList importDataIdentifier = interfaceSubdict.lookupType<wordList>("importDataIdentifier");
                for(auto rData : importDataIdentifier)
                {
                    interfacedata.importDataIdentifier.push_back(rData);
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

bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::is_same_points(Foam::vector& pointX, Foam::vector& pointY)
{
    if(pointX[0] == pointY[0] && pointX[1] == pointY[1] && pointX[2] == pointY[2])
        return true;
    else
        return false;
}

int Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::compare_nodes(Foam::vector& pointX, std::size_t interface_index)
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
            if(is_same_points(pointX , nodei))
            {
                auto itr = std::find( array_of_nodes.begin(), array_of_nodes.end(), nodei ) ;
                answer = std::distance(array_of_nodes.begin(), itr) + interfaces_.at(interface_index).globalNodeIndexBegin ;
                break;
            }
            else
            {
                answer = -1;
            }
        }
    }

    // Return position if node found, otherwise (-1)
    return answer;
}


// *********************************** Auxiliar functions for load calculation ************************************//
bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::calculateForces(std::size_t interface_index)
{
    // Get different force fields from OpenFOAM, Refering Force-Function Object (OpenFOAM User guide)
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

        tmp<vectorField> tempsurface = getFaceVectors(patchID);
        const auto& surface = tempsurface();

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

        // Writing this forces into vecotrs to export to CoSimulation
        int bufferIndex = 0;

        // For every cell in the interface-patch
        forAll(Force_->boundaryField()[patchID], i)
        {
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

Foam::tmp<Foam::volScalarField> Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::rho() const
{
    // If volScalarField exists, read it from registry (for compressible cases)
    if (mesh_.foundObject<volScalarField>("rho"))
    {
        return mesh_.lookupObject<volScalarField>("rho");
    }
    else if (solverType_.compare("incompressible") == 0)
    {
        const dictionary& FSIDict = dict_.subOrEmptyDict("parameters");

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
            const dictionary& FSIDict = dict_.subOrEmptyDict("parameters");

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

Foam::tmp<Foam::vectorField> Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::getFaceVectors(const unsigned int patchID) const
{
    // Normal vectors multiplied by face area
    return mesh_.boundary()[patchID].Sf();
}

void Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::conversionElementalToNodalValues(std::size_t interface_index)
{
    int number_of_nodes_in_element = 0;

    // Travel all Elements and Distribute the loads on the nodes
    for(auto& elementi : interfaces_.at(interface_index).Interface_elements)
    {
        number_of_nodes_in_element = elementi.getNumberOfNodes();

        for(auto& elementalNodeIndexi : elementi.getElementalNodeIndexes())
        {
            Node& temp_node = interfaces_.at(interface_index).Interface_nodes.at( elementalNodeIndexi - interfaces_.at(interface_index).globalNodeIndexBegin );
            std::vector<double>& temp_force = temp_node.getLoadValues();

            temp_force[0] += ( (interfaces_.at(interface_index).data_to_send[ (( elementi.getLocalElementIndex())- interfaces_.at(interface_index).globalElementIndexBegin ) * 3 + 0]) / double(number_of_nodes_in_element) );
            temp_force[1] += ( (interfaces_.at(interface_index).data_to_send[ (( elementi.getLocalElementIndex())- interfaces_.at(interface_index).globalElementIndexBegin ) * 3 + 1]) / double(number_of_nodes_in_element) );
            temp_force[2] += ( (interfaces_.at(interface_index).data_to_send[ (( elementi.getLocalElementIndex())- interfaces_.at(interface_index).globalElementIndexBegin ) * 3 + 2]) / double(number_of_nodes_in_element) );
        }

    }

    // Clear all the entries of "data_to_send" array
    interfaces_.at(interface_index).data_to_send.clear();

    // Resize the "data_to_send" to keep nodal data (number of nodes = number of local nodes for this rank)
    interfaces_.at(interface_index).data_to_send.resize( (( model_part_interfaces_.at(interface_index)->NumberOfNodes() ) * dim ), 0.0);

    // Fill the vector "data_to_send" with nodal load data. Filling data in Local ordering
    for(auto& nodei : interfaces_.at(interface_index).Interface_nodes)
    {
        interfaces_.at(interface_index).data_to_send[( ( nodei.getLocalNodeIndex() - interfaces_.at(interface_index).globalNodeIndexBegin ) * 3) + 0] = (nodei.getLoadValues()[0]) ;

        interfaces_.at(interface_index).data_to_send[( ( nodei.getLocalNodeIndex() - interfaces_.at(interface_index).globalNodeIndexBegin ) * 3) + 1] = (nodei.getLoadValues()[1]) ;

        interfaces_.at(interface_index).data_to_send[( ( nodei.getLocalNodeIndex() - interfaces_.at(interface_index).globalNodeIndexBegin ) * 3) + 2] = (nodei.getLoadValues()[2]) ;
    }

    // Set local load value of all the nodes to zero, else it will keep on accumulating
    for(auto& nodei : interfaces_.at(interface_index).Interface_nodes)
    {
        nodei.setLoadValuesToZeros();
    }

}


// *********************************** Auxiliar functions for interaction with CoSimulation *************************//
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
        debugInfo( "Reading mesh data from OF for coupling interface : " + interfaces_.at(j).nameOfInterface, debugLevel );

        // Arrays to exchange the number of nodes and Elements with all ranks
        scalarListList sendDataNumNodeIndex(Pstream::nProcs());
        scalarListList recvDataNumNodeIndex(Pstream::nProcs());
        scalarListList sendDataNumElementIndex(Pstream::nProcs());
        scalarListList recvDataNumElementIndex(Pstream::nProcs());

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
            sendDataNumElementIndex[i].setSize(1);
            sendDataNumElementIndex[i][0] = interfaces_.at(j).numElements ;
        }

        // Make CoSimIO::ModelPart and push in the array of model_part_interfaces
        model_part_interfaces_.push_back(CoSimIO::make_unique<CoSimIO::ModelPart>(interfaces_.at(j).nameOfInterface));

        // MPI Exchange to know "Node-start Index" for all the ranks
        Pstream::exchange<scalarList, scalar>(sendDataNumNodeIndex, recvDataNumNodeIndex);
        Pstream::exchange<scalarList, scalar>(sendDataNumElementIndex, recvDataNumElementIndex);

        // Default values of Serial run, playes important role in parallel simulations
        int nodeIndex = 1;
        int elemIndex = 1;

        for(int i = 0 ; i < MyRank; i++)
        {
            nodeIndex += recvDataNumNodeIndex[i][0];
            elemIndex += recvDataNumElementIndex[i][0];
        }
        interfaces_.at(j).globalNodeIndexBegin = nodeIndex;
        interfaces_.at(j).globalElementIndexBegin = elemIndex;

        // Accessing all the nodes in the Inteface and collect the nodal data to make the Mesh(Model Part),
        // which is then Export to CoSimulation
        for(std::size_t i = 0; i < interfaces_.at(j).patchIDs.size(); i++)
        {
            label patchIndex1 = mesh_.boundaryMesh().findPatchID((interfaces_.at(j).patchNames).at(i));
            const UList<label> &bfaceCells1 = mesh_.boundaryMesh()[patchIndex1].faceCells();
            label patchIndex2 = 0;

            // Travel through all the nodes in that MPI Rank
            forAll(bfaceCells1, bfacei1)
            {
                const label& faceID1 = mesh_.boundaryMesh()[interfaces_.at(j).patchIDs[i]].start() + bfacei1;

                Element new_element = Element(elemIndex++);

                forAll(mesh_.faces()[faceID1], nodei1)
                {
                    const label& nodeID1 = mesh_.faces()[faceID1][nodei1];
                    auto pointX = mesh_.points()[nodeID1];
                    std::vector<bool> is_common_node(TotalNumOfProcesses , 0);

                    // Return position if node found, otherwise (-1)
                    int result = compare_nodes(pointX , j);

                    if(result == (-1))
                    {
                        // For all Neighbouring Boundaries
                        forAll(mesh_.boundaryMesh() , ipatch)
                        {
                            word BCtype = mesh_.boundaryMesh().types()[ipatch];
                            if(BCtype == "processor")
                            {
                                patchIndex2 = ipatch;
                                const processorPolyPatch& pp = refCast<const processorPolyPatch>( mesh_.boundaryMesh()[patchIndex2] );
                                int nighbor_proc_id = pp.neighbProcNo();

                                const UList<label> &bfaceCells2 = mesh_.boundaryMesh()[patchIndex2].faceCells();

                                forAll(bfaceCells2, bfacei2)
                                {
                                    const label& faceID2 = mesh_.boundaryMesh()[patchIndex2].start() + bfacei2;
                                    forAll(mesh_.faces()[faceID2], nodei2)
                                    {
                                        const label& nodeID2 = mesh_.faces()[faceID2][nodei2];
                                        auto pointY = mesh_.points()[nodeID2];

                                        if(is_same_points(pointX,pointY) == 1 )
                                        {
                                            is_common_node.at(nighbor_proc_id) = 1;
                                        }
                                    }
                                }
                            }
                        }

                        bool new_node_created = 0;

                        for(int i = 0; i< TotalNumOfProcesses ; i++)
                        {
                            if(is_common_node.at(i) == 1)
                            {
                                if(new_node_created == 0)
                                {
                                    Node new_node = Node(pointX, nodeIndex, nodeIndex, i);
                                    new_node_created = 1;
                                    interfaces_.at(j).Interface_nodes.push_back(new_node);
                                    new_element.addNodeIndexInList(nodeIndex);
                                    new_element.addNodesInList(new_node);
                                    nodeIndex++;
                                    array_of_nodes.push_back(pointX);
                                }
                                else{
                                    interfaces_.at(j).Interface_nodes.at(nodeIndex-2).setCommonWithRank(i);
                                }
                            }
                        }
                        is_common_node.clear();

                        if(new_node_created == 0)
                        {
                            Node new_node = Node(pointX, nodeIndex, nodeIndex, MyRank);
                            interfaces_.at(j).Interface_nodes.push_back(new_node);
                            new_element.addNodeIndexInList(nodeIndex);
                            new_element.addNodesInList(new_node);
                            nodeIndex++;
                            array_of_nodes.push_back(pointX);
                        }

                    }
                    else{
                        // No need to produce new Node, just add in the list of nodes
                        new_element.addNodeIndexInList(result);
                        new_element.addNodesInList(interfaces_.at(j).Interface_nodes.at( result- interfaces_.at(j).globalNodeIndexBegin));
                    }
                }
                // Once Element is set, push it to the List of Elements
                interfaces_.at(j).Interface_elements.push_back(new_element);
            }

            //debugInfo( "[OF]Total number of Nodes in coupling interface "  + interfaces_.at(j).nameOfInterface + " are " +  std::to_string(interfaces_.at(j).numNodes), debugLevel);
            //debugInfo( "[OF]Total number of Elements in coupling interface "  + interfaces_.at(j).nameOfInterface + " are " +  std::to_string(interfaces_.at(j).numElements), debugLevel);

        }

        // Save the info about neighbour rank and corresponding number of nodes common with it
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
            num_common_nodes = 0 ;
        }

        // -----------------------Common Nodal Data exchange with all ranks --------------------------------//
        scalarListList sendNodalData(Pstream::nProcs());
        scalarListList recvNodalData(Pstream::nProcs());

        // Decide the size of a send array
        for(int i = 0; i < TotalNumOfProcesses ; i++)
        {
            sendNodalData[i].setSize(4 * interfaces_.at(j).neighbour_ids_comm_num_of_nodes.at( 2*i + 1 ));
        }

        // Filiing of the send array
        std::vector<int>counter(TotalNumOfProcesses , 0);
        for (auto& nodei: interfaces_.at(j).Interface_nodes)
        {
            for(std::size_t i = 0 ;  i < nodei.getCommonWithRank().size() ; i++)
            {
                if(nodei.getCommonWithRank()[i] == 1)
                {
                    sendNodalData[i][(counter.at(i))++] = nodei.getLocalNodeIndex();
                    vector nodePosition  = nodei.getNodePosition();
                    sendNodalData[i][(counter.at(i))++] = nodePosition[0];
                    sendNodalData[i][(counter.at(i))++] = nodePosition[1];
                    sendNodalData[i][(counter.at(i))++] = nodePosition[2];
                }
            }
        }
        counter.clear();

        // MPI_Exchange (Common Nodal data between all the ranks)
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
                            nodei.setNodeIndexForCoSim(recvNodalData[i][j*4+0]);
                        }
                    }
                }
            }
        }

        // Make CoSim nodes (Either Local of Ghost)
        for(auto& nodei : interfaces_.at(j).Interface_nodes)
        {
            if(nodei.getLocalNodeIndex() != nodei.getNodeIndexForCoSim() && (nodei.getHighestCommRank()<TotalNumOfProcesses)) // ghost nodes
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
        //debugInfo( "[COSIM]Total number of Nodes formed in coupling interface " + interfaces_.at(j).nameOfInterface +  " (local, Ghost, total) = (" + std::to_string(model_part_interfaces_.at(j)->NumberOfLocalNodes()) + " , " + std::to_string(model_part_interfaces_.at(j)->NumberOfGhostNodes()) +  " , " + std::to_string(model_part_interfaces_.at(j)->NumberOfNodes()) + " )" , debugLevel);

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

}

void Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::resizeDataVectors()
{
    for(std::size_t i=0; i < num_interfaces_; i++)
    {
        interfaces_.at(i).data_to_send.resize((interfaces_.at(i).numElements) * dim);
    }
}

void Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::exportDataToKratos()
{
    for(std::size_t i=0; i < num_interfaces_; i++)
    {
        debugInfo( "Force calculation started for coupling interface : " + interfaces_.at(i).nameOfInterface , debugLevel);

        // Resize to Elemental size, clear data before resizing
        interfaces_.at(i).data_to_send.clear();
        resizeDataVectors();

        // Force Calculations on Elements or cell centers
        calculateForces(i);

        // Conversion Utilities for converting Elemental Load values to Nodal values
        debugInfo( "Conversion of Elemental Force to Nodal Force : " + interfaces_.at(i).nameOfInterface , debugLevel);
        conversionElementalToNodalValues(i);

        // Export this force array to CoSimulation
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
        COSIMIO_CHECK_EQUAL(interfaces_.at(i).data_to_recv.size() , (interfaces_.at(i).numNodes) * dim );

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

}//namespace Foam

// ************************************************************************* //
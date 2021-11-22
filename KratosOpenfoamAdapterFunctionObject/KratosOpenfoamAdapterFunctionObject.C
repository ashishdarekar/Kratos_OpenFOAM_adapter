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

    // Add Try and catch block, Catch all possible error here, because read() is the methods
    // Which is going to executed before first time step and hence we can exit the simulation
    // if something is wrong here.

    Pout << "CoSimulation Adapter's function object : read()" << endl;

    try
    {
        // Reading configuration parameters from ControlDict related to function object
        fvMeshFunctionObject::read(dict);

        my_name = dict.lookupOrDefault<word>("participant", "fluid");
        Pout << "Name of the participant is: " << my_name <<endl;

        dim = dict.lookupOrDefault<int>("dim", 3);
        Pout << "Dimension of a problem is: " << dim <<endl;

        thick = dict.lookupOrDefault<double>("thick", 1.0);
        Pout << "Thickness of a domain is: " << thick <<endl;

        word name_of_interface;

        if(dict.lookupOrDefault("adjustTimeStep", false))
        {
            Pout << "Cannot support adjustable time step" << endl;
            //EXIT THE SIMULATION
        }

        // Check the solver type and determine it if needed
        solverType_ = dict.lookupOrDefault<word>("solvertype", "none");
        if (solverType_.compare("compressible") == 0 || solverType_.compare("incompressible") == 0)
        {
            Pout << "Known solver type: " << solverType_ << endl;
        }
        else if (solverType_.compare("none") == 0)
        {
            Pout << "Determining the solver type..." << endl;
            solverType_ = determineSolverType();
        }
        else
        {
            Pout << "Unknown solver type. Determining the solver type..." << endl;
            solverType_ = determineSolverType();
        }

        // Every interface is a subdictionary of "interfaces"
        const dictionary * interfaceSubdictPtr = dict.subDictPtr("interfaces");

        Pout << "Interfaces Reading: Start" << endl;

        if(!interfaceSubdictPtr)
        {
            Pout << "No Interfaces found" << endl;
            //EXIT THE SIMULATION
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
                    Pout << "Name of the interface is = " << interfacedata.nameOfInterface << endl;

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

        Pout << "Interfaces Reading: Done" << endl;
        Pout << "Number of interfaces found: " << num_interfaces_<< endl;

        Pout << "**************** Exporting InterfaceMesh using ModelPart: Start ******************" << endl;
        for(std::size_t j = 0; j < num_interfaces_; j++)
        {
            Pout << "Name of the interface under progress : " << interfaces_.at(j).nameOfInterface << endl;

            // *******************************Create a mesh as a ModelPart************************************ //
            std::vector<int> patchIDs;

            // For every patch that participates in the coupling interface. We are keeping one patch for one interface
            for (std::size_t i = 0; i < interfaces_.at(j).patchNames.size(); i++)
            {
                // Get the patchID
                int patchID = mesh_.boundaryMesh().findPatchID((interfaces_.at(j).patchNames).at(i));


                // Throw an error if the patch was not found
                if (patchID == -1)
                {
                    Pout << "ERROR: Patch " << (interfaces_.at(j).patchNames).at(i) << " does not exist." << endl;
                }

                // Add the patch in the list
                patchIDs.push_back(patchID);
            }

            // Count the Nodes for all the patches in that interface
            for (std::size_t i = 0; i < patchIDs.size(); i++)
            {
                interfaces_.at(j).numNodes += mesh_.boundaryMesh()[patchIDs.at(i)].localPoints().size();
            }
            Pout << "Total Number of Nodes in this interface: " << interfaces_.at(j).numNodes  << endl;

            // Count the number of elements/faces for all the patches in that interface
            for (std::size_t i = 0; i < patchIDs.size(); i++)
            {
                interfaces_.at(j).numElements += mesh_.boundary()[patchIDs[i]].size();
            }

            Pout << "Total Number of Elements/faces in this interface: " << interfaces_.at(j).numElements << endl;

            // Make CoSimIO::ModelPart and push in the array of model_part_interfaces
            model_part_interfaces_.push_back(CoSimIO::make_unique<CoSimIO::ModelPart>(interfaces_.at(j).nameOfInterface));

            // For Nodes and Element IDs for CoSimIO
            Pout << "Creating Model Part (Nodes and Elements) for CoSimIO : start" << endl;
            if(0)
            {
                int nodeIndex;
                if(MyRank ==0)
                {
                    nodeIndex = 1; //As Node indexing starts with 1 in CoSimIO
                }
                else
                {
                    nodeIndex = 273;

                }

                int elemIndex = 1; //As element indexing starts with 1 in CoSimIO

                //Customization
                vector point135(8.85, 5.97,  0.5);
                vector point136(8.85, 5.97, -0.5);
                vector point252(8.85, 6.03,  0.5);
                vector point251(8.85, 6.03, -0.5);

                // Accessing the coordinates of nodes in the Inteface and making CoSimIO nodes and elements
                for(std::size_t i = 0; i < patchIDs.size(); i++)
                {
                    label patchIndex1 = mesh_.boundaryMesh().findPatchID(interfaces_.at(j).patchNames[i]);
                    label patchIndex2 = 0; //Check the default value, What should I provide?

                    forAll(mesh_.boundaryMesh() , ipatch)
                    {
                        word BCtype = mesh_.boundaryMesh().types()[ipatch];
                        if( BCtype == "processor" )
                        {
                            patchIndex2 = ipatch;
                        }
                    }
                    const UList<label> &bfaceCells1 = mesh_.boundaryMesh()[patchIndex1].faceCells();
                    const UList<label> &bfaceCells2 = mesh_.boundaryMesh()[patchIndex2].faceCells();

                    int is_ghost_node = 0;
                    forAll(bfaceCells1, bfacei1)
                    {
                        const label& faceID1 = mesh_.boundaryMesh()[patchIDs[i]].start() + bfacei1;

                        std::vector<CoSimIO::IdType> connectivity;
                        forAll(mesh_.faces()[faceID1], nodei1)
                        {
                            const label& nodeID1 = mesh_.faces()[faceID1][nodei1]; //for OpenFOAM
                            auto pointX = mesh_.points()[nodeID1];

                            int result = compare_nodes(pointX); // return nodeIndex if node is already present and (-1) if node is not present

                            if(result == (-1)) // For new node
                            {
                                //While creation of node check if it has to be created in the way of ghost or local
                                forAll(bfaceCells2, bfacei2) //compare with all the nodes in the Processor common patch
                                {
                                    const label& faceID2 = mesh_.boundaryMesh()[patchIndex2].start() + bfacei2;
                                    forAll(mesh_.faces()[faceID2], nodei2)
                                    {
                                        const label& nodeID2 = mesh_.faces()[faceID2][nodei2]; //for OpenFOAM
                                        auto pointY = mesh_.points()[nodeID2];

                                        if(is_same_points(pointX,pointY) == 1 )//once this is found Go out from the For loop and create the node
                                        {
                                            Pout << "Found match in the nodeIDs with NodeID = "  << nodeID1 << " , " << nodeID2 << endl;
                                            Pout<<"{X,Y,Z = " << pointY[0] << " , " << pointY[1] << " , " <<pointY[2] << " }"<< endl;
                                            is_ghost_node = 1; //need to create this node as a gghost node
                                        }
                                    }
                                }

                                //Creation of the node either normal of Ghost
                                if(is_ghost_node && MyRank==1)
                                {
                                    if(is_same_points(pointX, point136))
                                    {
                                        Pout << "I m in point136 : Actaul node index going on = "<< nodeIndex <<endl;
                                        const processorPolyPatch& pp = refCast<const processorPolyPatch>( mesh_.boundaryMesh()[patchIndex2] );
                                        model_part_interfaces_.at(j)->CreateNewGhostNode( 136, pointX[0], pointX[1], pointX[2], pp.neighbProcNo());
                                        {Pout << "Making ghost Node with Index = "<< 136 <<endl;}
                                        is_ghost_node = 0; //make it zero again

                                        //array_of_nodes.push_back(pointX);// Push new node in the list to compare
                                        connectivity.push_back(136); // connectivity to make that element
                                        //nodeIndex ++;
                                    }else if(is_same_points(pointX, point135))
                                    {
                                        Pout << "I m in point135 : Actaul node index going on = "<< nodeIndex <<endl;
                                        const processorPolyPatch& pp = refCast<const processorPolyPatch>( mesh_.boundaryMesh()[patchIndex2] );
                                        model_part_interfaces_.at(j)->CreateNewGhostNode( 135, pointX[0], pointX[1], pointX[2], pp.neighbProcNo());
                                        {Pout << "Making ghost Node with Index = "<< 135 <<endl;}
                                        is_ghost_node = 0; //make it zero again

                                        //array_of_nodes.push_back(pointX);// Push new node in the list to compare
                                        connectivity.push_back(135); // connectivity to make that element
                                        //nodeIndex ++;
                                    }else if(is_same_points(pointX, point251))
                                    {
                                        Pout << "I m in point251 : Actaul node index going on = "<< nodeIndex <<endl;
                                        const processorPolyPatch& pp = refCast<const processorPolyPatch>( mesh_.boundaryMesh()[patchIndex2] );
                                        model_part_interfaces_.at(j)->CreateNewGhostNode( 251, pointX[0], pointX[1], pointX[2], pp.neighbProcNo());
                                        {Pout << "Making ghost Node with Index = "<< 251 <<endl;}
                                        is_ghost_node = 0; //make it zero again

                                        //array_of_nodes.push_back(pointX);// Push new node in the list to compare
                                        connectivity.push_back(251); // connectivity to make that element
                                        //nodeIndex ++;
                                    }else if(is_same_points(pointX, point252))
                                    {
                                        Pout << "I m in point252 : Actaul node index going on = "<< nodeIndex <<endl;
                                        const processorPolyPatch& pp = refCast<const processorPolyPatch>( mesh_.boundaryMesh()[patchIndex2] );
                                        model_part_interfaces_.at(j)->CreateNewGhostNode( 252, pointX[0], pointX[1], pointX[2], pp.neighbProcNo());
                                        {Pout << "Making ghost Node with Index = "<< 252 <<endl;}
                                        is_ghost_node = 0; //make it zero again

                                        //array_of_nodes.push_back(pointX);// Push new node in the list to compare
                                        connectivity.push_back(252); // connectivity to make that element
                                        //nodeIndex ++;
                                    }else
                                    {
                                        Pout << "I m in default : Actaul node index going on = "<< nodeIndex <<endl;
                                        const processorPolyPatch& pp = refCast<const processorPolyPatch>( mesh_.boundaryMesh()[patchIndex2] );
                                        model_part_interfaces_.at(j)->CreateNewGhostNode( nodeIndex, pointX[0], pointX[1], pointX[2], pp.neighbProcNo());

                                        is_ghost_node = 0; //make it zero again

                                        array_of_nodes.push_back(pointX);// Push new node in the list to compare
                                        connectivity.push_back(nodeIndex); // connectivity to make that element
                                        nodeIndex ++;
                                    }

                                    /* const processorPolyPatch& pp = refCast<const processorPolyPatch>( mesh_.boundaryMesh()[patchIndex2] );
                                    model_part_interfaces_.at(j)->CreateNewGhostNode( nodeIndex, pointX[0], pointX[1], pointX[2], pp.neighbProcNo());//ASk Ghost nodes ID should start from 1???
                                    if(1)
                                        {Pout << "Making Ghost Node with Index = "<< nodeIndex << " and original ID = " << nodeID1 << endl;}
                                    is_ghost_node = 0; //make it zero again

                                    array_of_nodes.push_back(pointX);// Push new node in the list to compare
                                    connectivity.push_back(nodeIndex); // connectivity to make that element
                                    nodeIndex++; */
                                }
                                else
                                {
                                    model_part_interfaces_.at(j)->CreateNewNode( nodeIndex, pointX[0], pointX[1], pointX[2]);
                                    if(MyRank ==1)
                                        {Pout << "Making Normal Node with Index = "<< nodeIndex << " and original ID = " << nodeID1 <<endl;}

                                    array_of_nodes.push_back(pointX);// Push new node in the list to compare
                                    connectivity.push_back(nodeIndex); // connectivity to make that element
                                    nodeIndex++;
                                }

                            }
                            else // Old node index just push in connectivity to make new element
                            {
                                if(MyRank ==1)
                                    {Pout << "Found Repetitive node = "<< result <<endl;}
                                if(MyRank ==1)
                                {
                                    result =result+272;
                                }
                                connectivity.push_back(result); // connectivity to make that element
                            }
                        }

                        model_part_interfaces_.at(j)->CreateNewElement( elemIndex, CoSimIO::ElementType::Quadrilateral2D4, connectivity );
                        if(MyRank ==1)
                            {Pout << "Making Element with Index = "<< elemIndex << ": using nodeindexes ={" << connectivity.at(0) << "," << connectivity.at(1) << "," <<connectivity.at(2) << "," <<connectivity.at(3) << "}" <<endl;}
                        elemIndex++;
                    }

                }

                Pout << "Name of the interface done : " << interfaces_.at(j).nameOfInterface << endl;

                Pout << "Total number of (local, Ghost, total) = (" << model_part_interfaces_.at(j)->NumberOfLocalNodes() << " , " <<model_part_interfaces_.at(j)->NumberOfGhostNodes() <<  " , "<<model_part_interfaces_.at(j)->NumberOfNodes() << " )"<<endl;
                Pout << "Total number of Elements formed: " << model_part_interfaces_.at(j)->NumberOfElements() <<endl;


                // iterate elements (with range based loop)
                /* for (auto& element : model_part_interfaces_.at(j)->Elements()) {
                    Pout<< "Element Id = " << element.Id() << " : Made up of following nodes: " <<endl;
                    for (auto node_it=element.NodesBegin(); node_it!=element.NodesEnd(); ++node_it) {
                        CoSimIO::Node& node = **node_it;
                        Pout<< "Node Id = " << node.Id() << " : with Co-ordinates (x,y,z) = ("<< node.X() << " , " << node.Y() << " , " << node.Z() << ")" << endl;
                    }
                }

                // iterate ghost nodes
                for (auto& ghost_node : model_part_interfaces_.at(j)->GhostNodes()) {
                    // do sth with node, e.g. print the id:
                    Pout<< "Ghost Node Id = " << ghost_node.Id() << " : with Co-ordinates (x,y,z) = ("<< ghost_node.X() << " , " << ghost_node.Y() << " , " << ghost_node.Z() << ")" << endl;
                } */

            }

            //New Code 14.11.2021
            if(0)
            {
                int nodeIndex = 1;
                int elemIndex = 1;
                int num_ghost_nodes=0;
                int nighbor_proc_id = 0;
                // Accessing the coordinates of nodes in the Inteface and making CoSimIO nodes and elements
                for(std::size_t i = 0; i < patchIDs.size(); i++)
                {
                    label patchIndex1 = mesh_.boundaryMesh().findPatchID(interfaces_.at(j).patchNames[i]);
                    label patchIndex2 = 0; //Check the default value, What should I provide?

                    forAll(mesh_.boundaryMesh() , ipatch)
                    {
                        word BCtype = mesh_.boundaryMesh().types()[ipatch];
                        if( BCtype == "processor" )
                        {
                            patchIndex2 = ipatch;
                        }
                    }
                    const UList<label> &bfaceCells1 = mesh_.boundaryMesh()[patchIndex1].faceCells();
                    const UList<label> &bfaceCells2 = mesh_.boundaryMesh()[patchIndex2].faceCells();
                    const processorPolyPatch& pp = refCast<const processorPolyPatch>( mesh_.boundaryMesh()[patchIndex2] );
                    nighbor_proc_id = pp.neighbProcNo();

                    int is_ghost_node = 0;

                    forAll(bfaceCells1, bfacei1)
                    {
                        const label& faceID1 = mesh_.boundaryMesh()[patchIDs[i]].start() + bfacei1;
                        Element new_element = Element(elemIndex++);//make new element and then increase the elemental index counter

                        forAll(mesh_.faces()[faceID1], nodei1)
                        {
                            const label& nodeID1 = mesh_.faces()[faceID1][nodei1]; //for OpenFOAM
                            auto pointX = mesh_.points()[nodeID1];

                            int result = compare_nodes(pointX); // return nodeIndex(starting from 1) if node is already present and (-1) if node is not present

                            if(result == (-1)) // For new node
                            {
                                //While creation of node check if it has to be created in the way of ghost or local
                                forAll(bfaceCells2, bfacei2) //compare with all the nodes in the Processor common patch
                                {
                                    const label& faceID2 = mesh_.boundaryMesh()[patchIndex2].start() + bfacei2;
                                    forAll(mesh_.faces()[faceID2], nodei2)
                                    {
                                        const label& nodeID2 = mesh_.faces()[faceID2][nodei2]; //for OpenFOAM
                                        auto pointY = mesh_.points()[nodeID2];

                                        if(is_same_points(pointX,pointY) == 1 )//once this is found Go out from the For loop and create the node
                                        {
                                            //Pout << "Found match in the nodeIDs with NodeID = "  << nodeID1 << " , " << nodeID2 << endl;
                                            //Pout<<"{X,Y,Z = " << pointY[0] << " , " << pointY[1] << " , " <<pointY[2] << " }"<< endl;
                                            is_ghost_node = 1; //need to create this node as a gghost node
                                        }
                                    }
                                }

                                if(is_ghost_node == 1){
                                    //Node new_node = Node(pointX, nodeID1, nodeIndex, 1);
                                    Node new_node = Node(pointX, nodeIndex, nodeIndex, 1);//cosimIndex kept same as local index (by default), later need to change
                                    Interface_nodes.push_back(new_node);
                                    new_element.addNodeIndexInList(nodeIndex);// connectivity to make that element
                                    new_element.addNodesInList(new_node);// Add new node in the element
                                    nodeIndex++;
                                    array_of_nodes.push_back(pointX);// Push new node in the list to compare
                                    is_ghost_node = 0; //For next Node
                                }
                                else{
                                    Node new_node = Node(pointX, nodeIndex, nodeIndex, 0);//cosimIndex kept same as local index (by default), later need to change
                                    Interface_nodes.push_back(new_node);
                                    new_element.addNodeIndexInList(nodeIndex);// connectivity to make that element
                                    new_element.addNodesInList(new_node);// Add new node in the element
                                    nodeIndex++;
                                    array_of_nodes.push_back(pointX);// Push new node in the list to compare
                                }
                            }
                            else{
                                //No need to produce new Node
                                new_element.addNodeIndexInList(result);// connectivity to make that element
                                new_element.addNodesInList(Interface_nodes.at(result-1));// Add old node in the List of nodes for that element, result-1 as Vector index starts from 0
                            }
                        }
                        //Once Element is set push it to the List of Elements
                        Interface_elements.push_back(new_element);
                    }

                }

                Pout << "Name of the interface done : " << interfaces_.at(j).nameOfInterface << endl;

                //Pout << "Total number of (local, Ghost, total) = (" << model_part_interfaces_.at(j)->NumberOfLocalNodes() << " , " <<model_part_interfaces_.at(j)->NumberOfGhostNodes() <<  " , "<<model_part_interfaces_.at(j)->NumberOfNodes() << " )"<<endl;
                //Pout << "Total number of Elements formed: " << model_part_interfaces_.at(j)->NumberOfElements() <<endl;
                Pout << "[OF]Total number of Nodes formed: " << Interface_nodes.size() <<endl;
                Pout << "[OF]Total number of Elements formed: " << Interface_elements.size() <<endl;

                for(auto& nodei: Interface_nodes)
                {
                    if(nodei.getGhostStatus())
                    {
                        num_ghost_nodes++;
                    }
                }

                Pout << "[OF]Total number of Ghost Nodes formed: " << num_ghost_nodes<<endl;

                int shift = Interface_nodes.size() - num_ghost_nodes ; //272-4 = 268 (for proc 0), //56-4 = 52(for proc 1)
                int pass_ghost_counter = 0;

                //Trying MPI_Send and MPI_Recv
                scalarListList sendData(Pstream::nProcs());
                scalarListList recvData(Pstream::nProcs());
                //fill in the Send Data
                sendData[0].setSize(16);
                sendData[1].setSize(16);

                int counter = 0;

                //Making CoSim Nodes
                if(MyRank == 0) //Ghost only in 0
                {
                    for(auto& nodei : Interface_nodes)
                    {
                        if(nodei.getGhostStatus())
                        {
                            vector nodePosition  = nodei.getNodePosition();
                            nodei.setNodeIndexForCoSim(++shift);

                            model_part_interfaces_.at(j)->CreateNewGhostNode( nodei.getNodeIndexForCoSim(), nodePosition[0], nodePosition[1], nodePosition[2], nighbor_proc_id);

                            //if(MyRank == 0){Pout << "Making ghost Node with LocalIndex = "<<  nodei.getLocalNodeIndex() << " and CoSimNodeIndex = " << nodei.getNodeIndexForCoSim() <<endl;}
                            pass_ghost_counter ++;
                        }
                        else
                        {
                            vector nodePosition  = nodei.getNodePosition();
                            nodei.setNodeIndexForCoSim(nodei.getLocalNodeIndex() - pass_ghost_counter);
                            model_part_interfaces_.at(j)->CreateNewNode( nodei.getNodeIndexForCoSim(), nodePosition[0], nodePosition[1], nodePosition[2]);
                            //if(MyRank == 0) {Pout << "Making Normal Node with LocalIndex = "<<  nodei.getLocalNodeIndex() << " and CoSimNodeIndex = " << nodei.getNodeIndexForCoSim() <<endl;}
                        }

                    }

                    Pout << "[COSIM]Total number of Nodes(local, Ghost, total) = (" << model_part_interfaces_.at(j)->NumberOfLocalNodes() << " , " <<model_part_interfaces_.at(j)->NumberOfGhostNodes() <<  " , "<<model_part_interfaces_.at(j)->NumberOfNodes() << " )"<<endl;

                    //Making CoSim elements
                    std::vector<CoSimIO::IdType> connectivity;
                    for(auto& elementi : Interface_elements)
                    {
                        /* for(auto elementalNodeIndexi : elementi.getElementalNodeIndexes())
                        {
                            connectivity.push_back(elementalNodeIndexi);
                        } */
                        for(auto& elementalNodeIndexi : elementi.getElementalNodes())
                        {
                            connectivity.push_back(Interface_nodes.at(elementalNodeIndexi.getNodeIndexForCoSim()-1).getNodeIndexForCoSim());
                        }

                        model_part_interfaces_.at(j)->CreateNewElement( elementi.getLocalElementIndex(), CoSimIO::ElementType::Quadrilateral2D4, connectivity );
                        connectivity.clear();
                    }

                    Pout << "[COSIM]Total number of Elements = " << model_part_interfaces_.at(j)->NumberOfElements() <<endl;

                    // iterate ghost nodes
                    for (auto& ghost_node : model_part_interfaces_.at(j)->GhostNodes()) {
                        // do sth with node, e.g. print the id:
                        //Pout<< "Ghost Node Id = " << ghost_node.Id() << " : with Co-ordinates (x,y,z) = ("<< ghost_node.X() << " , " << ghost_node.Y() << " , " << ghost_node.Z() << ")" << endl;

                        sendData[1][counter++] = ghost_node.Id();
                        sendData[1][counter++] = ghost_node.X();
                        sendData[1][counter++] = ghost_node.Y();
                        sendData[1][counter++] = ghost_node.Z();
                    }


                }

                //MPI related stuff
                Pstream::exchange<scalarList, scalar>(sendData, recvData);

                if(MyRank == 1) //Proc = 1
                {
                    // Received data from proc 0
                    for(int i = 0; i<recvData[0].size();i+=4) {
                        // do sth with node, e.g. print the id:
                        //Pout<< "RECV Ghost Node Id = " << recvData[0][i] << " : with Co-ordinates (x,y,z) = ("<< recvData[0][i+1]  << " , " << recvData[0][i+2]  << " , " << recvData[0][i+3]  << ")" << endl;
                    }


                    for(auto& nodei : Interface_nodes)
                    {
                        //vector trial_node;
                        if(nodei.getGhostStatus())
                        {
                            for(int i = 0; i<recvData[0].size()/4 ; i++)
                            {
                                vector trial_node(recvData[0][i*4+1], recvData[0][i*4+2], recvData[0][i*4+3]);

                                vector nodePosition  = nodei.getNodePosition();

                                if (is_same_points(trial_node, nodePosition))
                                {
                                    nodei.setNodeIndexForCoSim(recvData[0][i*4+0]); //Take node Id of that

                                    model_part_interfaces_.at(j)->CreateNewNode( nodei.getNodeIndexForCoSim(), nodePosition[0], nodePosition[1], nodePosition[2]);

                                    //if(MyRank == 1){Pout << "Making local_ghost Node with LocalIndex = "<<  nodei.getLocalNodeIndex() << " and CoSimNodeIndex = " << nodei.getNodeIndexForCoSim() <<endl;}
                                    pass_ghost_counter ++;

                                }
                            }

                        }
                        else
                        {
                            vector nodePosition  = nodei.getNodePosition();
                            nodei.setNodeIndexForCoSim(272 + nodei.getLocalNodeIndex() - pass_ghost_counter);
                            model_part_interfaces_.at(j)->CreateNewNode( nodei.getNodeIndexForCoSim(), nodePosition[0], nodePosition[1], nodePosition[2]);
                            //if(MyRank == 1) {Pout << "Making Normal Node with LocalIndex = "<<  nodei.getLocalNodeIndex() << " and CoSimNodeIndex = " << nodei.getNodeIndexForCoSim() <<endl;}
                        }

                    }

                    Pout << "[COSIM]Total number of Nodes(local, Ghost, total) = (" << model_part_interfaces_.at(j)->NumberOfLocalNodes() << " , " <<model_part_interfaces_.at(j)->NumberOfGhostNodes() <<  " , "<<model_part_interfaces_.at(j)->NumberOfNodes() << " )"<<endl;

                    //Making CoSim elements
                    std::vector<CoSimIO::IdType> connectivity;
                    for(auto& elementi : Interface_elements)
                    {
                        for(auto& elementalNodeIndexi : elementi.getElementalNodes())
                        {
                            //Pout<< "Value of a index = " << elementalNodeIndexi.getNodeIndexForCoSim()<< " and "<<  Interface_nodes.at(elementalNodeIndexi.getNodeIndexForCoSim()-1).getNodeIndexForCoSim() << endl;
                            connectivity.push_back(Interface_nodes.at(elementalNodeIndexi.getNodeIndexForCoSim()-1).getNodeIndexForCoSim());
                        }
                        model_part_interfaces_.at(j)->CreateNewElement( elementi.getLocalElementIndex(), CoSimIO::ElementType::Quadrilateral2D4, connectivity );
                        connectivity.clear();
                    }

                    Pout << "[COSIM]Total number of Elements = " << model_part_interfaces_.at(j)->NumberOfElements() <<endl;

                }

                //Printing some things
                if(MyRank == 0 && 0)
                {
                    // Iterate over all the nodes
                    for(auto& node: model_part_interfaces_.at(j)->Nodes())
                    {
                        Pout<< "Node Id = " << node.Id() << " : with Co-ordinates (x,y,z) = ("<< node.X() << " , " << node.Y() << " , " << node.Z() << ")" << endl;
                    }
                    /*
                    // iterate elements (with range based loop)
                    for (auto& element : model_part_interfaces_.at(j)->Elements()) {
                        Pout<< "Element Id = " << element.Id() << " : Made up of following nodes: " <<endl;
                        for (auto node_it=element.NodesBegin(); node_it!=element.NodesEnd(); ++node_it) {
                            CoSimIO::Node& node = **node_it;
                            Pout<< "Node Id = " << node.Id() << " : with Co-ordinates (x,y,z) = ("<< node.X() << " , " << node.Y() << " , " << node.Z() << ")" << endl;
                        }
                    }

                    // iterate ghost nodes
                    for (auto& ghost_node : model_part_interfaces_.at(j)->GhostNodes()) {
                        // do sth with node, e.g. print the id:
                        Pout<< "Ghost Node Id = " << ghost_node.Id() << " : with Co-ordinates (x,y,z) = ("<< ghost_node.X() << " , " << ghost_node.Y() << " , " << ghost_node.Z() << ")" << endl;
                    } */
                }
            }

            if(0)
            {
                //Trying MPI_Send and MPI_Recv
                scalarListList sendData(Pstream::nProcs());
                scalarListList recvData(Pstream::nProcs());
                //fill in the Send Data
                sendData[0].setSize(4);
                sendData[1].setSize(4);


                if(MyRank == 0)
                {
                    sendData[0][0] = 1;
                    sendData[0][1] = 2;
                    sendData[0][2] = 3;
                    sendData[0][3] = 4;

                    sendData[1][0] = 8;
                    sendData[1][1] = 9;
                    sendData[1][2] = 10;
                    sendData[1][3] = 11;

                }

                if(MyRank == 1)
                {
                    sendData[0][0] = 11;
                    sendData[0][1] = 21;
                    sendData[0][2] = 31;
                    sendData[0][3] = 41;

                    sendData[1][0] = 81;
                    sendData[1][1] = 91;
                    sendData[1][2] = 101;
                    sendData[1][3] = 111;

                }


                Pstream::exchange<scalarList, scalar>(sendData, recvData);

                //Check the receive data now
                if(MyRank == 1)
                {
                    Pout<<"Value of the element = " <<recvData[0][0] << endl;
                    Pout<<"Value of the element = " <<recvData[0][1] << endl;
                    Pout<<"Value of the element = " <<recvData[0][2] << endl;
                    Pout<<"Value of the element = " <<recvData[0][3] << endl;

                    Pout<<"Value of the element = " <<recvData[1][0] << endl;
                    Pout<<"Value of the element = " <<recvData[1][1] << endl;
                    Pout<<"Value of the element = " <<recvData[1][2] << endl;
                    Pout<<"Value of the element = " <<recvData[1][3] << endl;
                }
            }


            //New code 19.11.2021
            if(1)
            {
                //Trying MPI_Send and MPI_Recv
                scalarListList sendData(Pstream::nProcs());
                scalarListList recvData(Pstream::nProcs());
                //fill in the Send Data
                sendData[0].setSize(16);
                sendData[1].setSize(16);
                int nodeIndex = 1 ; //Default value of Serial run
                if(MyRank == 0)
                {
                    nodeIndex = 1;
                }
                else//MyRank ==1
                {
                    nodeIndex = 273;
                }
                int elemIndex = 1;
                int num_ghost_nodes=0;
                int nighbor_proc_id = 0;
                // Accessing the coordinates of nodes in the Inteface and making CoSimIO nodes and elements
                for(std::size_t i = 0; i < patchIDs.size(); i++)
                {
                    label patchIndex1 = mesh_.boundaryMesh().findPatchID(interfaces_.at(j).patchNames[i]);
                    label patchIndex2 = 0; //Check the default value, What should I provide?

                    forAll(mesh_.boundaryMesh() , ipatch)
                    {
                        word BCtype = mesh_.boundaryMesh().types()[ipatch];
                        if( BCtype == "processor" )
                        {
                            patchIndex2 = ipatch;
                        }
                    }
                    const UList<label> &bfaceCells1 = mesh_.boundaryMesh()[patchIndex1].faceCells();
                    const UList<label> &bfaceCells2 = mesh_.boundaryMesh()[patchIndex2].faceCells();
                    const processorPolyPatch& pp = refCast<const processorPolyPatch>( mesh_.boundaryMesh()[patchIndex2] );
                    nighbor_proc_id = pp.neighbProcNo();

                    int is_ghost_node = 0;

                    forAll(bfaceCells1, bfacei1)
                    {
                        const label& faceID1 = mesh_.boundaryMesh()[patchIDs[i]].start() + bfacei1;
                        Element new_element = Element(elemIndex++);//make new element and then increase the elemental index counter

                        forAll(mesh_.faces()[faceID1], nodei1)
                        {
                            const label& nodeID1 = mesh_.faces()[faceID1][nodei1]; //for OpenFOAM
                            auto pointX = mesh_.points()[nodeID1];

                            int result = compare_nodes(pointX); // return nodeIndex(starting from 1) if node is already present and (-1) if node is not present

                            if(result == (-1)) // For new node
                            {
                                //While creation of node check if it has to be created in the way of ghost or local
                                forAll(bfaceCells2, bfacei2) //compare with all the nodes in the Processor common patch
                                {
                                    const label& faceID2 = mesh_.boundaryMesh()[patchIndex2].start() + bfacei2;
                                    forAll(mesh_.faces()[faceID2], nodei2)
                                    {
                                        const label& nodeID2 = mesh_.faces()[faceID2][nodei2]; //for OpenFOAM
                                        auto pointY = mesh_.points()[nodeID2];

                                        if(is_same_points(pointX,pointY) == 1 )//once this is found Go out from the For loop and create the node
                                        {
                                            //Pout << "Found match in the nodeIDs with NodeID = "  << nodeID1 << " , " << nodeID2 << endl;
                                            //Pout<<"{X,Y,Z = " << pointY[0] << " , " << pointY[1] << " , " <<pointY[2] << " }"<< endl;
                                            is_ghost_node = 1; //need to create this node as a gghost node
                                        }
                                    }
                                }

                                if(is_ghost_node == 1){
                                    //Node new_node = Node(pointX, nodeID1, nodeIndex, 1);
                                    Node new_node = Node(pointX, nodeIndex, nodeIndex, 1);//cosimIndex kept same as local index (by default), later need to change
                                    Interface_nodes.push_back(new_node);
                                    new_element.addNodeIndexInList(nodeIndex);// connectivity to make that element
                                    new_element.addNodesInList(new_node);// Add new node in the element
                                    nodeIndex++;
                                    array_of_nodes.push_back(pointX);// Push new node in the list to compare
                                    is_ghost_node = 0; //For next Node
                                }
                                else{
                                    Node new_node = Node(pointX, nodeIndex, nodeIndex, 0);//cosimIndex kept same as local index (by default), later need to change
                                    Interface_nodes.push_back(new_node);
                                    new_element.addNodeIndexInList(nodeIndex);// connectivity to make that element
                                    new_element.addNodesInList(new_node);// Add new node in the element
                                    nodeIndex++;
                                    array_of_nodes.push_back(pointX);// Push new node in the list to compare
                                }
                            }
                            else{
                                //No need to produce new Node
                                new_element.addNodeIndexInList(result);// connectivity to make that element
                                new_element.addNodesInList(Interface_nodes.at(result-1));// Add old node in the List of nodes for that element, result-1 as Vector index starts from 0
                            }
                        }
                        //Once Element is set push it to the List of Elements
                        Interface_elements.push_back(new_element);
                    }

                }

                Pout << "Name of the interface done : " << interfaces_.at(j).nameOfInterface << endl;

                //Pout << "Total number of (local, Ghost, total) = (" << model_part_interfaces_.at(j)->NumberOfLocalNodes() << " , " <<model_part_interfaces_.at(j)->NumberOfGhostNodes() <<  " , "<<model_part_interfaces_.at(j)->NumberOfNodes() << " )"<<endl;
                //Pout << "Total number of Elements formed: " << model_part_interfaces_.at(j)->NumberOfElements() <<endl;
                Pout << "[OF]Total number of Nodes formed: " << Interface_nodes.size() <<endl;
                Pout << "[OF]Total number of Elements formed: " << Interface_elements.size() <<endl;

                for(auto& nodei: Interface_nodes)
                {
                    if(nodei.getGhostStatus())
                    {
                        //Pout <<"[OF]Ghost Node with the Id = " << nodei.getLocalNodeIndex() <<endl;
                        num_ghost_nodes++;
                    }
                }

                Pout << "[OF]Total number of Ghost Nodes formed: " << num_ghost_nodes<<endl;

                // Get the Ghost data from Proc 1 to Proc 0
                int counter = 0;
                if(MyRank == 1)
                {
                    // iterate ghost nodes
                    for (auto& nodei: Interface_nodes)
                    {
                        if(nodei.getGhostStatus() == 1)
                        {
                            sendData[0][counter++] = nodei.getLocalNodeIndex();
                            vector nodePosition  = nodei.getNodePosition();
                            sendData[0][counter++] = nodePosition[0];
                            sendData[0][counter++] = nodePosition[1];
                            sendData[0][counter++] = nodePosition[2];
                        }
                    }
                }

                //MPI related stuff
                Pstream::exchange<scalarList, scalar>(sendData, recvData);

                if(MyRank == 0)
                {
                    // Received data from proc 1
                    /* for(int i = 0; i<recvData[1].size();i+=4) {
                        Pout<< "RECV Ghost Node Id = " << recvData[1][i] << " : with Co-ordinates (x,y,z) = ("<< recvData[1][i+1]  << " , " << recvData[1][i+2]  << " , " << recvData[1][i+3]  << ")" << endl;
                    } */

                    // Change the NodeIndexes for CoSimIO for the Ghost nodes in Proc 0 once with the received Data
                    for (auto& nodei: Interface_nodes)
                    {
                        if(nodei.getGhostStatus() == 1)
                        {
                            vector nodePosition  = nodei.getNodePosition();

                            for(int i = 0; i<recvData[0].size()/4 ; i++)
                            {
                                vector trial_node(recvData[1][i*4+1], recvData[1][i*4+2], recvData[1][i*4+3]);

                                if(is_same_points(nodePosition,trial_node))
                                {
                                    nodei.setNodeIndexForCoSim(recvData[1][i*4+0]); //Take node Id of that
                                }
                            }
                        }
                    }

                    for(auto& nodei: Interface_nodes)
                    {
                        if(nodei.getGhostStatus())
                        {
                            //Pout <<"[OF]Ghost Node with the Original Id = " << nodei.getLocalNodeIndex() << " and ID for CoSim = " << nodei.getNodeIndexForCoSim() <<endl;
                            num_ghost_nodes++;
                        }
                    }

                }

                //Now make the CoSim Nodes and Elements

                //Making CoSim Nodes
                if(MyRank == 0) //Ghost only in 0
                {
                    for(auto& nodei : Interface_nodes)
                    {
                        if(nodei.getGhostStatus())
                        {
                            vector nodePosition  = nodei.getNodePosition();
                            model_part_interfaces_.at(j)->CreateNewGhostNode( nodei.getNodeIndexForCoSim(), nodePosition[0], nodePosition[1], nodePosition[2], nighbor_proc_id);
                            //Pout << "Making ghost Node with LocalIndex = "<<  nodei.getLocalNodeIndex() << " and CoSimNodeIndex = " << nodei.getNodeIndexForCoSim() <<endl;
                        }
                        else
                        {
                            vector nodePosition  = nodei.getNodePosition();
                            model_part_interfaces_.at(j)->CreateNewNode( nodei.getNodeIndexForCoSim(), nodePosition[0], nodePosition[1], nodePosition[2]);
                            //Pout << "Making Normal Node with LocalIndex = "<<  nodei.getLocalNodeIndex() << " and CoSimNodeIndex = " << nodei.getNodeIndexForCoSim() <<endl;
                        }

                    }

                    Pout << "[COSIM]Total number of Nodes(local, Ghost, total) = (" << model_part_interfaces_.at(j)->NumberOfLocalNodes() << " , " <<model_part_interfaces_.at(j)->NumberOfGhostNodes() <<  " , "<<model_part_interfaces_.at(j)->NumberOfNodes() << " )"<<endl;

                    //Making CoSim elements
                    std::vector<CoSimIO::IdType> connectivity;
                    for(auto& elementi : Interface_elements)
                    {
                        for(auto& elementalNodeIndexi : elementi.getElementalNodes())
                        {
                            connectivity.push_back(Interface_nodes.at(elementalNodeIndexi.getLocalNodeIndex()-1).getNodeIndexForCoSim());
                        }

                        model_part_interfaces_.at(j)->CreateNewElement( elementi.getLocalElementIndex(), CoSimIO::ElementType::Quadrilateral2D4, connectivity );
                        connectivity.clear();
                    }

                    Pout << "[COSIM]Total number of Elements = " << model_part_interfaces_.at(j)->NumberOfElements() <<endl;

                }
                else //For rank =1
                {
                    for(auto& nodei : Interface_nodes)
                    {
                        vector nodePosition  = nodei.getNodePosition();
                        model_part_interfaces_.at(j)->CreateNewNode( nodei.getNodeIndexForCoSim(), nodePosition[0], nodePosition[1], nodePosition[2]);
                        //Pout << "Making Normal Node with LocalIndex = "<<  nodei.getLocalNodeIndex() << " and CoSimNodeIndex = " << nodei.getNodeIndexForCoSim() <<endl;
                    }

                    Pout << "[COSIM]Total number of Nodes(local, Ghost, total) = (" << model_part_interfaces_.at(j)->NumberOfLocalNodes() << " , " <<model_part_interfaces_.at(j)->NumberOfGhostNodes() <<  " , "<<model_part_interfaces_.at(j)->NumberOfNodes() << " )"<<endl;

                    //Making CoSim elements
                    std::vector<CoSimIO::IdType> connectivity;
                    for(auto& elementi : Interface_elements)
                    {
                        for(auto& elementalNodeIndexi : elementi.getElementalNodes())
                        {
                            connectivity.push_back(Interface_nodes.at(elementalNodeIndexi.getLocalNodeIndex()-272-1).getNodeIndexForCoSim());
                        }

                        model_part_interfaces_.at(j)->CreateNewElement( elementi.getLocalElementIndex(), CoSimIO::ElementType::Quadrilateral2D4, connectivity );
                        connectivity.clear();
                    }

                    Pout << "[COSIM]Total number of Elements = " << model_part_interfaces_.at(j)->NumberOfElements() <<endl;

                }

                //Printing some things
                if(MyRank == 1 && 0)
                {
                    // Iterate over all the nodes
                    for(auto& node: model_part_interfaces_.at(j)->Nodes())
                    {
                        Pout<< "Node Id = " << node.Id() << " : with Co-ordinates (x,y,z) = ("<< node.X() << " , " << node.Y() << " , " << node.Z() << ")" << endl;
                    }

                    // iterate elements (with range based loop)
                    for (auto& element : model_part_interfaces_.at(j)->Elements()) {
                        Pout<< "Element Id = " << element.Id() << " : Made up of following nodes: " <<endl;
                        for (auto node_it=element.NodesBegin(); node_it!=element.NodesEnd(); ++node_it) {
                            CoSimIO::Node& node = **node_it;
                            Pout<< "Node Id = " << node.Id() << " : with Co-ordinates (x,y,z) = ("<< node.X() << " , " << node.Y() << " , " << node.Z() << ")" << endl;
                        }
                    }

                    // iterate ghost nodes
                    for (auto& ghost_node : model_part_interfaces_.at(j)->GhostNodes()) {
                        // do sth with node, e.g. print the id:
                        Pout<< "Ghost Node Id = " << ghost_node.Id() << " : with Co-ordinates (x,y,z) = ("<< ghost_node.X() << " , " << ghost_node.Y() << " , " << ghost_node.Z() << ")" << endl;
                    }
                }

            }

            // Connection between openFOAM and Kratos-CoSimulation using CoSimIO
            CoSimIO::Info settings;
            settings.Set("my_name", "Openfoam_Adapter");
            settings.Set("connect_to", "Openfoam_Kratos_Wrapper");
            settings.Set("echo_level", 0);
            settings.Set("version", "1.25");

            auto connect_info = CoSimIO::ConnectMPI(settings, MPI_COMM_WORLD);
            COSIMIO_CHECK_EQUAL(connect_info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Connected);
            connection_name = connect_info.Get<std::string>("connection_name");

            // Export InterfaceMesh/ModelPart to CoSimulation using CoSimIO
            Pout << "Exporting Mesh as a ModelPart for an interface: " << interfaces_.at(j).nameOfInterface << endl;
            info.Clear();
            info.Set("identifier", interfaces_.at(j).nameOfInterface);
            info.Set("connection_name", connection_name);
            auto export_info = CoSimIO::ExportMesh(info, *model_part_interfaces_.at(j));

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
            Pout << "Name of the interface Finished processing : " << interfaces_.at(j).nameOfInterface << endl;
            array_of_nodes.clear(); //Clear all the entries of array_nodes, so that new empty vector will be available for the comparision

        }
        Pout << "**************** Exporting InterfaceMesh using ModelPart: End ********************" << endl;

    }

    catch(const std::exception& e)
    {
        std::cerr << "Exception" << e.what() << '\n';
    }

    return true;
}

bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::execute()
{
    Pout << "CoSimulation Adapter's function object : execution()" << endl;

    // *************************************** Force/Load Related ****************************************** //
    for(std::size_t i=0; i < num_interfaces_; i++)
    {
        Pout << "Force calculation started for the interface : " << interfaces_.at(i).nameOfInterface << endl;

        // For "Write Data" variables which need to send to CoSimulation
        for(std::size_t j=0; j< interfaces_.at(i).exportData.size(); j++)
        {
            std::string dataName = interfaces_.at(i).exportData.at(j);

            if(dataName.find("Force") == 0 )
            {
                calculateForces(i);
            }
        }
        Pout << "Force calculation ended for the interface : " << interfaces_.at(i).nameOfInterface << endl;

        // Export this force array to CoSimulation //Elemental Force Data
        CoSimIO::Info connect_info;
        connect_info.Clear();
        connect_info.Set("identifier", interfaces_.at(i).exportDataIdentifier[0]);
        connect_info.Set("connection_name", connection_name);
        connect_info = CoSimIO::ExportData(connect_info, interfaces_.at(i).data_to_send);

        for(std::size_t p=0; p<interfaces_.at(i).data_to_send.size(); p++)
        {
            Pout<<"force output : " <<interfaces_.at(i).data_to_send.at(p)<<endl;
        }

        Pout << runTime_.timeName() << " : Data has been exported from OpenFOAM to CoSimulation (interface name = " << interfaces_.at(i).nameOfInterface << ") , Force values with array size = " << interfaces_.at(i).data_to_send.size() << endl;
    }

    // *************************************** Displcement Related ****************************************** //
    for (std::size_t i = 0; i < interfaces_.size(); i++)
    {
        // Import the displacement array from the CoSimulation
        CoSimIO::Info connect_info;
        connect_info.Clear();
        connect_info.Set("identifier", interfaces_.at(i).importDataIdentifier[0]);
        connect_info.Set("connection_name", connection_name);
        connect_info = CoSimIO::ImportData(connect_info, interfaces_.at(i).data_to_recv);
        COSIMIO_CHECK_EQUAL(interfaces_.at(i).data_to_recv.size() , (interfaces_.at(i).numNodes) * dim ); //Check size of Receive data = Number of nodes*dim Is it require??

        Pout << runTime_.timeName() << " : Data has been imported from CoSimulation to OpenFOAM: (interface name = " << interfaces_.at(i).nameOfInterface << ") , Disp values with array size = " << interfaces_.at(i).data_to_recv.size() << endl;

        for(std::size_t p=0; p<interfaces_.at(i).data_to_recv.size(); p++)
        {
            Pout<<"disp intput : " <<interfaces_.at(i).data_to_recv.at(p)<<endl;
        }

        Pout << "Displacement replacement started for the interface : " << interfaces_.at(i).nameOfInterface << endl;

        // Get the displacement on the patch(for every patch in the interface) and assign it those values received from CoSimulation
        for (std::size_t j = 0; j < interfaces_.at(i).patchNames.size(); j++)
        {
            Foam::pointVectorField* point_disp;
            point_disp = const_cast<pointVectorField*>( &mesh_.lookupObject<pointVectorField>("pointDisplacement") );
            label patchIndex = mesh_.boundaryMesh().findPatchID(interfaces_.at(i).patchNames[j]);//Remove hardcoded part for finding patchIndex
            fixedValuePointPatchVectorField& pointDisplacementFluidPatch = refCast<fixedValuePointPatchVectorField>(point_disp->boundaryFieldRef()[patchIndex]);

            int iterator = 0;
            forAll(point_disp->boundaryFieldRef()[patchIndex] ,k)
            {
                //Pout<< "Displ node number = " << k <<endl;
                pointDisplacementFluidPatch[k][0] = interfaces_.at(i).data_to_recv[iterator++];
                pointDisplacementFluidPatch[k][1] = interfaces_.at(i).data_to_recv[iterator++];
                if (dim ==3)
                    pointDisplacementFluidPatch[k][2] = interfaces_.at(i).data_to_recv[iterator++];
            }
        }

        Pout << "Displacement replacement ended for the interface : " << interfaces_.at(i).nameOfInterface << endl;
    }

    return true;
}

bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::end()
{
    // Dicsonect from CoSimIO
    Pout << "CoSimulation Adapter's function object : end()" << endl;

    // DisConnection between OpenFOAM and Kratos-CoSimulation using CoSimIO
    CoSimIO::Info connect_info;
    CoSimIO::Info disconnect_settings;
    disconnect_settings.Set("connection_name", connection_name);
    connect_info = CoSimIO::Disconnect(disconnect_settings);
    COSIMIO_CHECK_EQUAL(connect_info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Disconnected);

    return true;
}

bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::write()
{
    Pout << "CoSimulation Adapter's function object : write()" << endl;

    return true;
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
            Pout << "Solver Type : Compressible " << endl;
        }
        else if (p_.dimensions() == pressureDimensionsIncompressible)
        {
            solverType_ = "incompressible";
            Pout << "Solver Type : InCompressible " << endl;
        }
    }

    if(solverType_  == "unknown")
    {
        Pout << "Solver Type: Neither Compressible nor Incompresible" << endl;
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
    std::vector<int> patchIDs;
    // For every patch that participates in the coupling interface
    for (std::size_t i = 0; i < interfaces_.at(interface_index).patchNames.size(); i++)
    {
        // Get the patchID
        int patchID = mesh_.boundaryMesh().findPatchID((interfaces_.at(interface_index).patchNames).at(i));

        // Throw an error if the patch was not found
        if (patchID == -1){
            Pout << "ERROR: Patch " << (interfaces_.at(interface_index).patchNames).at(i) << " does not exist." << endl;
        }

    // Add the patch in the list
    patchIDs.push_back(patchID);
    }

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
    for(std::size_t j=0; j< patchIDs.size(); j++)
    {
        int patchID = patchIDs.at(j);

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

}//namespace Foam

// ************************************************************************* //
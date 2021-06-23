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

    //Connect to CoSimIO
    //CoSimulationAdapter_.configure();
    //std::cout << "CoSimulation Adapter : Configuration" << std::endl;

    std::cout << "CoSimulation Adapter's function object : read()" << std::endl;

    //Reading configuration parameters from ControlDict related to Adapter
    fvMeshFunctionObject::read(dict);

    my_name = dict.lookupOrDefault<word>("participant", "fluid");
    std::cout<< "Name of the participant is: " << my_name <<std::endl;

    //Check the solver type and determine it if needed
    solverType_ = dict.lookupOrDefault<word>("solvertype", "none");
    if (solverType_.compare("compressible") == 0 || solverType_.compare("incompressible") == 0)
    {
        std::cout << "Known solver type: " << solverType_ << std::endl;
    }
    else if (solverType_.compare("none") == 0)
    {
        std::cout << "Determining the solver type..." << std::endl;
        solverType_ = determineSolverType();
    }
    else
    {
        std::cout << "Unknown solver type. Determining the solver type..." << std::endl;
        solverType_ = determineSolverType();
    }

    //Check the total number of registered objects in the PolyMesh Object Registry related to Given Solver
    /*std::cout<< "Name of all registered objects in Foam::PolyMesh object Registry are:" << std::endl;
    Foam::wordList Objectnames_ = mesh_.names(); //List of all Objects in the polymesh class::mesh_object
    forAll(Objectnames_,i)
    {
        std::cout << Objectnames_[i] << ", ";
    }
    std::cout<<"\n";
    */

    // Every interface is a subdictionary of "interfaces",
    // each with an arbitrary name. Read all of them and create a list of dictionaries.
    const dictionary * interfaceSubdictPtr = dict.subDictPtr("interfaces");

    std::cout << "Interfaces Reading: Start" << std::endl;

    if(!interfaceSubdictPtr)
    {
        std::cout << "No Interfaces found" << std::endl;
    }
    else
    {
        for(const entry& interfaceSubdictEntry : *interfaceSubdictPtr)
        {
            if(interfaceSubdictEntry.isDict())
            {
                dictionary interfaceSubdict = interfaceSubdictEntry.dict();
                struct InterfaceData interfacedata;

                interfacedata.meshName = interfaceSubdict.lookupType<word>("mesh");

                interfacedata.locationsType = interfaceSubdict.lookupOrDefault<word>("locations", "faceCenters"); //Default value is faceCenter

                wordList patches = interfaceSubdict.lookupType<wordList>("patches");
                for(auto patch : patches)
                {
                    interfacedata.patchNames.push_back(patch);
                }

                wordList readData = interfaceSubdict.lookupType<wordList>("readData");
                for(auto rData : readData)
                {
                    interfacedata.readData.push_back(rData);
                }

                wordList writeData = interfaceSubdict.lookupType<wordList>("writeData");
                for(auto wData : writeData)
                {
                    interfacedata.writeData.push_back(wData);
                }

                //Add this interface in the Array of all interfaces
                interfaces_.push_back(interfacedata);

                num_interfaces_++;
            }
        }
    }

    std::cout << "Interfaces Reading: Done" << std::endl;
    std::cout << "Number of interfaces found: " << num_interfaces_<< std::endl;

    //Connection between openFOAM and Kratos-CoSimulation using CoSimIO
    CoSimIO::Info settings;
    settings.Set("my_name", "Openfoam_Adapter");
    settings.Set("connect_to", "Openfoam_Kratos_Wrapper");
    settings.Set("echo_level", 0);
    settings.Set("version", "1.25");

    auto connect_info = CoSimIO::Connect(settings);
    COSIMIO_CHECK_EQUAL(connect_info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Connected);
    connection_name = connect_info.Get<std::string>("connection_name");

    //*********************** Export Fluid-OpenFOAM Mesh to CoSimulation using CoSimIO *********************//
    /* //Create one CoSimIO::ModelPart for each coupling interface
    std::cout << "\n" <<"********************** Exporting InterfaceMesh using ModelPart: Start **********************" << "\n" <<std::endl;

    for(std::size_t j=0; j < num_interfaces_; j++)
    {
        std::string interface_name = "interface" + std::to_string(j+1);

        //model_part_interfaces_.push_back(CoSimIO::ModelPart(interface_name));
        model_part_interfaces_.push_back(CoSimIO::make_unique<CoSimIO::ModelPart>(interface_name));

        std::cout<< "Name of an Interface/model part: " << model_part_interfaces_.at(j)->Name() << std::endl;

        // **************************Create a mesh as a ModelPart********************************/
/*        std::cout << "Accessing Mesh from openFOAM" << std::endl;
        std::vector<int> patchIDs;

        // For every patch that participates in the coupling interface
        for (uint i = 0; i < interfaces_.at(j).patchNames.size(); i++)
        {
            // Get the patchID
            int patchID = mesh_.boundaryMesh().findPatchID((interfaces_.at(j).patchNames).at(i));

            // Throw an error if the patch was not found
            if (patchID == -1)
            {
                std::cout<< "ERROR: Patch " << (interfaces_.at(j).patchNames).at(i) << " does not exist." << std::endl;
            }

        // Add the patch in the list
        patchIDs.push_back(patchID);
        }

        // Count the Nodes for all the patches in that interface
        for (uint i = 0; i < patchIDs.size(); i++)
        {
            interfaces_.at(j).numNodes += mesh_.boundaryMesh()[patchIDs.at(i)].localPoints().size();
        }
        std::cout << "Total Number of Nodes in this interface: " << interfaces_.at(j).numNodes  << std::endl;

        // Count the number of elements/faces for all the patches in that interface
        for (uint i = 0; i < patchIDs.size(); i++)
        {
            interfaces_.at(j).numElements += mesh_.boundary()[patchIDs[i]].size();
        }
        std::cout << "Total Number of Elements/faces in this interface: " << interfaces_.at(j).numElements << std::endl;

        std::cout << "Creating Model Part (Nodes and Elements) for CoSimIO" << std::endl;

        //-For Nodes and Element IDs for CoSimIO
        std::vector<int> NodeIDs;
        NodeIDs.resize( interfaces_.at(j).numNodes );
        int nodeIndex = 1; //As Node indexing starts with 1 in CoSimIO

        std::vector<int> ElemIDs;
        ElemIDs.resize(interfaces_.at(j).numElements);
        int elemIndex = 1; //As element indexing starts with 1 in CoSimIO

        vector pointX(0,0,0); //For accessing the Co-ordinates of Nodes

        for(uint i = 0; i < patchIDs.size(); i++)
        {
            //const word& patchName = mesh_.boundary()[patchIDs[i]].name();
            //std::cout << "patchID: " << patchIDs[i] << " with name "  << patchName << " With size = "<< mesh_.boundary()[patchIDs[i]].size() << std::endl;

            forAll (mesh_.boundary()[patchIDs[i]],facei)
            {
                const label& faceID = mesh_.boundaryMesh()[patchIDs[i]].start() + facei;
                //std::cout << "ID of this face: " << faceID << " Which face: "<< facei << std::endl;
                std::vector<CoSimIO::IdType> connectivity;
                forAll (mesh_.faces()[faceID], nodei)
                {
                    const label& nodeID = mesh_.faces()[faceID][nodei]; //for OpenFOAM
                    pointX = mesh_.points()[nodeID];
                    //std::cout << "Node "<< nodei << " with Id=" << nodeID << ":" << pointX[0] << "," << pointX[1] << "," << pointX[2] << std::endl;
                    //-Make CoSimIO Nodes
                    NodeIDs.push_back(nodeIndex); //Later used to make CoSimIO::Element
                    CoSimIO::Node& node = model_part_interfaces_.at(j)->CreateNewNode( nodeIndex, pointX[0], pointX[1], pointX[2]);
                    connectivity.push_back(nodeIndex); //connectivity to make that element /face
                    nodeIndex++;
                }
                //-Make CoSimIO Element = currently hardcoded to make Quadrilateral3D4
                ElemIDs.push_back(elemIndex); //Maybe useful
                CoSimIO::Element& element = model_part_interfaces_.at(j)->CreateNewElement( elemIndex, CoSimIO::ElementType::Quadrilateral2D4, connectivity );
                elemIndex++;
            }
        }
        std::cout << "Converting InterfaceMesh to a CoSimIO::ModelPart -> End" << std::endl;
        // **************************Done: Create a mesh as a ModelPart********************************/
/*
        //Export InterfaceMesh/ModelPart to CoSimulation using CoSimIO
        std::cout << "Exporting Mesh as a ModelPart: Start for an interface = " << j+1 << std::endl;
        CoSimIO::Info info;
        info.Clear();
        info.Set("identifier", "interface_flap");
        info.Set("connection_name", connection_name);
        auto export_info = CoSimIO::ExportMesh(info, *model_part_interfaces_.at(j));
        std::cout << "Exporting Mesh as a ModelPart: End for an interface = " << j+1 << "\n " <<std::endl;
    }

    std::cout << "*********************** Exporting InterfaceMesh using ModelPart: End ************************" << "\n" <<std::endl;
*/

    // To test the CoSimulation Import Mesh Fuctationality with one interface only

    std::cout << "\n" <<"********************** Exporting InterfaceMesh using ModelPart: Start **********************" << "\n" <<std::endl;
    for(std::size_t j=0; j < num_interfaces_; j++)
    {
        std::string interface_name = "interface" + std::to_string(j+1);

        //No need to this ...AS we will create only single Model Part later in the code
        //model_part_interfaces_.push_back(CoSimIO::ModelPart(interface_name));
        //model_part_interfaces_.push_back(CoSimIO::make_unique<CoSimIO::ModelPart>(interface_name));

        //std::cout<< "Name of an Interface/model part: " << model_part_interfaces_.at(j)->Name() << std::endl;

        // **************************Create a mesh as a ModelPart********************************/
        std::cout << "Accessing Mesh from openFOAM" << std::endl;
        std::vector<int> patchIDs;

        // For every patch that participates in the coupling interface
        for (uint i = 0; i < interfaces_.at(j).patchNames.size(); i++)
        {
            // Get the patchID
            int patchID = mesh_.boundaryMesh().findPatchID((interfaces_.at(j).patchNames).at(i));

            // Throw an error if the patch was not found
            if (patchID == -1)
            {
                std::cout<< "ERROR: Patch " << (interfaces_.at(j).patchNames).at(i) << " does not exist." << std::endl;
            }

        // Add the patch in the list
        patchIDs.push_back(patchID);
        }

        // Count the Nodes for all the patches in that interface
        for (uint i = 0; i < patchIDs.size(); i++)
        {
            interfaces_.at(j).numNodes += mesh_.boundaryMesh()[patchIDs.at(i)].localPoints().size();
        }
        std::cout << "Total Number of Nodes in this interface: " << interfaces_.at(j).numNodes  << std::endl;

        // Count the number of elements/faces for all the patches in that interface
        for (uint i = 0; i < patchIDs.size(); i++)
        {
            interfaces_.at(j).numElements += mesh_.boundary()[patchIDs[i]].size();
            //interfaces_.at(j).numElements += mesh.boundaryMesh()[patchIDs.at(i)].faceCentres().size();
        }
        std::cout << "Total Number of Elements/faces in this interface: " << interfaces_.at(j).numElements << std::endl;

        //Considering 1st interface ONLY
        if(j==0)
        {
            // create CoSimIO::ModelPart- Put in public Member function
            CoSimIO::ModelPart model_part_interface_flap("interface_flap_model_part");

            //-For Nodes and Element IDs for CoSimIO
            std::cout << "Creating Model Part (Nodes and Elements) for CoSimIO : start" << std::endl;
            std::vector<int> NodeIDs;
            NodeIDs.resize( interfaces_.at(j).numNodes );
            int nodeIndex = 2; //As Node indexing starts with 1 in CoSimIO, 1st push is made

            std::vector<int> ElemIDs;
            ElemIDs.resize(interfaces_.at(j).numElements);
            int elemIndex = 1; //As element indexing starts with 1 in CoSimIO

            // For accessing the Co-ordinates of Nodes : Make 1st Node manually and push it
            //vector pointX(5.5000000e+00, 5.9700000e+00, -5.0000000e-06);
            vector pointX(5.5000000e+00, 5.9700000e+00, -5.0000000e-03);
            NodeIDs.push_back(1);
            CoSimIO::Node& node = model_part_interface_flap.CreateNewNode( 1, pointX[0], pointX[1], pointX[2]);
            array_of_nodes.push_back(pointX); //Initial pushback in array of nodes

            for(uint i = 0; i < patchIDs.size(); i++)
            {
                //const word& patchName = mesh_.boundary()[patchIDs[i]].name();
                //std::cout << "patchID: " << patchIDs[i] << " with name "  << patchName << " With size = "<< mesh_.boundary()[patchIDs[i]].size() << std::endl;

                forAll(mesh_.boundary()[patchIDs[i]],facei)
                {
                    const label& faceID = mesh_.boundaryMesh()[patchIDs[i]].start() + facei;
                    //std::cout << "ID of this face: " << faceID << " Which face: "<< facei << "; with number of nodes = " << mesh_.faces()[faceID].size() << std::endl;

                    std::vector<CoSimIO::IdType> connectivity;
                    forAll(mesh_.faces()[faceID], nodei)
                    {
                        const label& nodeID = mesh_.faces()[faceID][nodei]; //for OpenFOAM
                        pointX = mesh_.points()[nodeID];
                        //std::cout << "nodei: " << nodei << " With points: "<< pointX[0] << " | " << pointX[1] << " | " << pointX[2] << std::endl;

                        int result = compare_nodes(pointX); //nodeIndex if node is already present and (-1) if node is not present
                        if(result == (-1)) //new node
                        {
                            //-Make CoSimIO Nodes
                            NodeIDs.push_back(nodeIndex); //Later used to make CoSimIO::Element //NodeINdex = (facei*4)+ (nodei+1)
                            CoSimIO::Node& node = model_part_interface_flap.CreateNewNode( nodeIndex, pointX[0], pointX[1], pointX[2]);
                            array_of_nodes.push_back(pointX);//Push new element in the list to compare
                            //std::cout << "Node with index= " << nodeIndex << " is created "<<std::endl;
                            connectivity.push_back(nodeIndex); //connectivity to make that element /face
                            nodeIndex++;
                        }
                        else //old node index just push in connectivity to ake new element
                        {
                            connectivity.push_back(result); //connectivity to make that element /face
                        }
                    }
                    //-Make CoSimIO Element = currently hardcoded to make Quadrilateral2D4
                    ElemIDs.push_back(elemIndex); //Maybe useful
                    CoSimIO::Element& element = model_part_interface_flap.CreateNewElement( elemIndex, CoSimIO::ElementType::Quadrilateral2D4, connectivity );
                    //std::cout << "Element with index= " << elemIndex << " is created "<<std::endl;//elemIndex =facei+1
                    elemIndex++;
                }
            }
            std::cout << "\n" << "Number of nodes: " << nodeIndex -1 << " and Elements: " << elemIndex -1 << std::endl;
            //std::cout << "\n" << "Size of Array of nodes = " << array_of_nodes.size() <<std::endl;

            std::cout << "Converting InterfaceMesh to a CoSimIO::ModelPart -> End" << std::endl;

            //Printing each node and element sequence wise
            /* std::cout << "Nodes in Model Part representing Mesh" << std::endl;
            for (auto node_it=model_part_interface_flap.NodesBegin(); node_it!=model_part_interface_flap.NodesEnd(); ++node_it)
            {
                CoSimIO::Node& node = **node_it;
                std::cout<< "nodeID= " << node.Id() << " = "<< node.X() << " | " << node.Y() << " | " << node.Z() << std::endl;
            } */

            //std::cout << "Elements and corresponding Nodes in Model Part representing Mesh" << std::endl;
            for (auto elem_it=model_part_interface_flap.ElementsBegin(); elem_it!=model_part_interface_flap.ElementsEnd(); ++elem_it)
            {
                CoSimIO::Element& element = **elem_it;
                //std::cout << "All nodes in element with ID =" << element.Id() << " , with number of nodes =  " << element.NumberOfNodes() << std::endl;
                my_element_class new_elem;
                new_elem.ID = element.Id()-1;

                for (auto node_it=element.NodesBegin(); node_it!=element.NodesEnd(); ++node_it)
                {
                    CoSimIO::Node& node = **node_it;
                    //std::cout<< "nodeID= " << node.Id() << " = "<< node.X() << " | " << node.Y() << " | " << node.Z() << std::endl;

                    my_node_class new_node;
                    new_node.ID = node.Id() -1;
                    new_node.posx = node.X();
                    new_node.posy = node.Y();
                    new_node.posz = node.Z();
                    new_node.fx = 0.0;
                    new_node.fy = 0.0;
                    new_node.fz = 0.0;
                    new_elem.node_list.push_back(new_node);
                }

                AllElements_.push_back(new_elem);
            }

            std::cout << "Nodal Force Calculation" << std::endl;
            for (std::size_t i = 0 ; i < AllElements_.size(); i++)
            {
                std::cout << "All nodes in element with ID =" << AllElements_.at(i).ID << " , with number of nodes =  " << AllElements_.at(i).node_list.size() << std::endl;
                for (std::size_t j = 0 ; j< AllElements_.at(i).node_list.size() ; j++)
                {
                    std::cout<< "nodeID= " << AllElements_.at(i).node_list[j].ID << " = "<< AllElements_.at(i).node_list[j].posx << " | " << AllElements_.at(i).node_list[j].posy << " | " << AllElements_.at(i).node_list[j].posz << std::endl;
                }
            }

            //Export InterfaceMesh/ModelPart to CoSimulation using CoSimIO
            std::cout << "Exporting Mesh as a ModelPart: Start for an interface flap" << std::endl;
            CoSimIO::Info info;
            info.Clear();
            info.Set("identifier", "interface_flap");
            info.Set("connection_name", connection_name);
            auto export_info = CoSimIO::ExportMesh(info, model_part_interface_flap);
            std::cout << "*********************** Exporting InterfaceMesh using ModelPart: End ************************" << "\n" <<std::endl;
        }

    }


    //Making Data Vectors on the interfaces only
    for(std::size_t i=0; i < num_interfaces_; i++)
    {
        //For "Wirte Data" variables which need to send to CoSimulation
        for(std::size_t j=0; j< interfaces_.at(i).writeData.size(); j++)
        {
            std::string dataName = interfaces_.at(i).writeData.at(j);

            if(dataName.find("Force") == 0 || dataName.find("Stress") == 0) //If "force" or "stress" string is found it will return 0
            {
                if(interfaces_.at(i).locationsType == "faceNodes")
                {
                    data_to_send.resize((interfaces_.at(i).numNodes) * dim);
                }
                else if(interfaces_.at(i).locationsType == "faceCenters")
                {
                    data_to_send.resize((interfaces_.at(i).numElements) * dim);
                    Nodal_Force_data_to_send.resize((interfaces_.at(i).numNodes) * dim); //Changed on 17/06/2021. To check How nodal Forces works on FSI problem
                }
            }
            //else if() //if some other variables
        }

        //For "Read Data" variables which need to receive from CoSimulation
        //No need to resize it. Structual Solver will send this data
        for(std::size_t j=0; j< interfaces_.at(i).readData.size(); j++)
        {
            std::string dataName = interfaces_.at(i).readData.at(j);

            if(dataName.find("Displacement") == 0 || dataName.find("DisplacementDelta") == 0) //If "Displacement" or "DisplacementDelta" string is found it will return 0
            {
                if(interfaces_.at(i).locationsType == "faceNodes")
                {
                    //data_to_recv.resize((interfaces_.at(i).numNodes) * dim);
                }
                else if(interfaces_.at(i).locationsType == "faceCenters")
                {
                    //data_to_recv.resize((interfaces_.at(i).numElements) * dim);
                }
            }
            //else if() //if some other variables
        }
    }

    return true;
}

bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::execute()
{
    std::cout << "CoSimulation Adapter's function object : execution()" << std::endl;

    if(mesh_.foundObject<pointVectorField>("pointDisplacement"))
    {
        // *************************************** Displcement Related ****************************************** //
        //Import the displacement array from the CoSimulation
        CoSimIO::Info connect_info;
        connect_info.Clear();
        connect_info.Set("identifier", "disp");
        connect_info.Set("connection_name", connection_name);
        connect_info = CoSimIO::ImportData(connect_info, data_to_recv);
        //Check size of Receive data = Expected Receive data. Get it from top.
        /* for(uint i = 0; i < data_to_recv.size() ; i++ )
        {
            std::cout << "i = " << i << " ; disp data = " << data_to_recv.at(i) << std::endl;
        } */

        std::cout<< runTime_.timeName() << " : Data has been imported from CoSimulation to OpenFOAM: Disp values with array size = " << data_to_recv.size() << std::endl;

        // Get the displacement on the patch and assign it those values received from CoSimulation,
        Foam::pointVectorField* point_disp;
        point_disp = const_cast<pointVectorField*>( &mesh_.lookupObject<pointVectorField>("pointDisplacement") );
        //std::cout<< "Size of the pointDisplacement array is " << point_disp->size() << std::endl;
        label patchIndex = mesh_.boundaryMesh().findPatchID("flap");

        std::cout<< "Patch for displacement " << patchIndex << std::endl;

        fixedValuePointPatchVectorField& pointDisplacementFluidPatch = refCast<fixedValuePointPatchVectorField>(point_disp->boundaryFieldRef()[patchIndex]);

        int iterator = 0 ; //iterator goes from 0 till (size of recv_array)
        forAll(point_disp->boundaryFieldRef()[patchIndex] ,i)
        {
            //std::cout << i << " Before : "<< pointDisplacementFluidPatch[i][0] << ", " << pointDisplacementFluidPatch[i][1] << ", " << pointDisplacementFluidPatch[i][2] << std::endl;
            pointDisplacementFluidPatch[i][0] = data_to_recv[iterator++];
            pointDisplacementFluidPatch[i][1] = data_to_recv[iterator++];
            //if (dim ==3)
                //pointDisplacementFluidPatch[i][2] = data_to_recv[iterator++]; //ignoring z direction displacement = Anyways the value is zero for 2D case.
            iterator++;
            //std::cout << i << " After : "<< pointDisplacementFluidPatch[i][0] << ", " << pointDisplacementFluidPatch[i][1] << ", " << pointDisplacementFluidPatch[i][2] << std::endl;
        }
        //std::cout << "\n" << "Size of the iterator = " << iterator <<std::endl;

        // *************************************** Force/Load Related ****************************************** //
        //At the end of time loop, Calculation of the Forces, which need to send to structural solver in FSI problem
        //Total force = Viscous force + Pressure force
        std::cout<< "Force calculation : start" << std::endl;
        for(std::size_t i=0; i < num_interfaces_; i++)
        {
            //For "Wirte Data" variables which need to send to CoSimulation
            for(std::size_t j=0; j< interfaces_.at(i).writeData.size(); j++)
            {
                std::string dataName = interfaces_.at(i).writeData.at(j);
                std::cout<< "interface = "<< i+1 << " with DataName = " << dataName << std::endl;

                if(dataName.find("Force") == 0 )
                {
                    calculateForces(j);
                }
            }
        }

        std::cout<< "Force calculation : Done" << std::endl;
        //std::cout<< "Data to be send from OpenFOAM: with array size = " << data_to_send.size() << std::endl;

        /* for(uint i = 0; i < data_to_send.size() ; i++ )
        {
            std::cout << "i = " << i << " ; load data = " << data_to_send.at(i) << std::endl;
        } */

        //Force Mapping on Nodes - 17.06.2021
        //************************************ ONLY WORK FOR THIS MESH ***********************************//
        int mapping_required = 0;
        while(mapping_required)
        {
            std::cout << "DO NOT CHANGE THE MESH : Nodal Force Calculation : Start" << std::endl;
            for (std::size_t i = 0 ; i < AllElements_.size(); i++)
            {
                if(i==0)
                {
                    for(std::size_t j = 0; j < AllElements_[i].node_list.size(); j++)
                    {
                        if (j==0 || j==1)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 ; //No Cell On Right side
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 ; //No Cell On Right side
                            //Z component No Change
                        }
                        else if (j==2 || j==3)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 + (data_to_send[3*(i+1) + 0])/4 ;
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 + (data_to_send[3*(i+1) + 1])/4 ;
                            //Z component No Change
                        }
                    }
                }
                //for i=1 till
                else if(i>=1 || i<=79)
                {
                    for(std::size_t j = 0; j < AllElements_[i].node_list.size(); j++)
                    {
                        if (j==0 || j==1)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 + (data_to_send[3*(i-1) + 0])/4 ;
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 + (data_to_send[3*(i-1) + 1])/4 ;
                            //Z component No Change
                        }
                        else if (j==2 || j==3)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 + (data_to_send[3*(i+1) + 0])/4 ;
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 + (data_to_send[3*(i+1) + 1])/4 ;
                            //Z component No Change
                        }
                    }
                }
                else if(i==80)
                {
                    for(std::size_t j = 0; j < AllElements_[i].node_list.size(); j++)
                    {
                        if (j==0 || j==1)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 + (data_to_send[3*(i-1) + 0])/4 ;
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 + (data_to_send[3*(i-1) + 1])/4 ;
                            //Z component No Change
                        }
                        else if (j==2 || j==3)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 + (data_to_send[3*(150) + 0])/4 ; //looking for cell# 150
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 + (data_to_send[3*(150) + 1])/4 ; //looking for cell# 150
                            //Z component No Change
                        }
                    }
                }
                else if(i==81)
                {
                    for(std::size_t j = 0; j < AllElements_[i].node_list.size(); j++)
                    {
                        if (j==0 || j==3)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 + (data_to_send[3*(160) + 0])/4 ; //looking for cell# 160
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 + (data_to_send[3*(160) + 1])/4 ; //looking for cell# 160
                            //Z component No Change
                        }
                        else if (j==1 || j==2)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 + (data_to_send[3*(i+1) + 0])/4 ;
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 + (data_to_send[3*(i+1) + 1])/4 ;
                            //Z component No Change
                        }
                    }
                }
                else if((i>=82 && i<=149) || (i>=152 && i<=159))
                {
                    for(std::size_t j = 0; j < AllElements_[i].node_list.size(); j++)
                    {
                        if (j==0 || j==3)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 + (data_to_send[3*(i-1) + 0])/4 ;
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 + (data_to_send[3*(i-1) + 1])/4 ;
                            //Z component No Change
                        }
                        else if (j==1 || j==2)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 + (data_to_send[3*(i+1) + 0])/4 ;
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 + (data_to_send[3*(i+1) + 1])/4 ;
                            //Z component No Change
                        }
                    }
                }
                else if(i==150)
                {
                    for(std::size_t j = 0; j < AllElements_[i].node_list.size(); j++)
                    {
                        if (j==0 || j==3)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 + (data_to_send[3*(i-1) + 0])/4 ;
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 + (data_to_send[3*(i-1) + 1])/4 ;
                            //Z component No Change
                        }
                        else if (j==1 || j==2)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 + (data_to_send[3*(80) + 0])/4 ; //looking for cell# 80
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 + (data_to_send[3*(80) + 1])/4 ; //looking for cell# 80
                            //Z component No Change
                        }
                    }
                }
                else if(i==151)
                {
                    for(std::size_t j = 0; j < AllElements_[i].node_list.size(); j++)
                    {
                        if (j==0 || j==3)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 ; //No Cell On left side
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 ; //No Cell On left side
                            //Z component No Change
                        }
                        else if (j==1 || j==2)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 + (data_to_send[3*(i+1) + 0])/4 ;
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 + (data_to_send[3*(i+1) + 1])/4 ;
                            //Z component No Change
                        }
                    }

                }
                else if(i==160)
                {
                    for(std::size_t j = 0; j < AllElements_[i].node_list.size(); j++)
                    {
                        if (j==0 || j==3)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 + (data_to_send[3*(i-1) + 0])/4 ;
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 + (data_to_send[3*(i-1) + 1])/4 ;
                            //Z component No Change
                        }
                        else if (j==1 || j==2)
                        {
                            AllElements_[i].node_list[j].fx = (data_to_send[3*i + 0])/4 + (data_to_send[3*(81) + 0])/4 ;
                            AllElements_[i].node_list[j].fy = (data_to_send[3*i + 1])/4 + (data_to_send[3*(81) + 1])/4 ;
                            //Z component No Change
                        }
                    }
                }
            }

            //Write new Modal Data into new Buffer- Nodal_Force_data_to_send
            for (std::size_t i = 0 ; i < AllElements_.size(); i++)
            {
                for(std::size_t j = 0; j < AllElements_[i].node_list.size(); j++)
                {
                    Nodal_Force_data_to_send[3*(AllElements_[i].node_list[j].ID) + 0] = AllElements_[i].node_list[j].fx ;
                    Nodal_Force_data_to_send[3*(AllElements_[i].node_list[j].ID) + 1] = AllElements_[i].node_list[j].fy ;
                    Nodal_Force_data_to_send[3*(AllElements_[i].node_list[j].ID) + 2] = 0.0;//Z Component for each Node
                }
            }

            std::cout << "DO NOT CHANGE THE MESH : Nodal Force Calculation: Done" << std::endl;
            /* for(uint i = 0; i < Nodal_Force_data_to_send.size() ; i++ )
            {
                std::cout << "i = " << i << " ; load data = " << Nodal_Force_data_to_send.at(i) << std::endl;
            } */
        }

        //************************************ ONLY WORK FOR THIS MESH ***********************************//

        //Export this force array to CoSimulation
        //CoSimIO::Info connect_info;
        connect_info.Clear();
        connect_info.Set("identifier", "load_cells");
        connect_info.Set("connection_name", connection_name);
        connect_info = CoSimIO::ExportData(connect_info, data_to_send);//Elemental Force Data
        //connect_info = CoSimIO::ExportData(connect_info, Nodal_Force_data_to_send); //Nodal Force Data

        std::cout<< runTime_.timeName() << " : Data has been exported from OpenFOAM to CoSimulation: Force values with array size = " << data_to_send.size() << std::endl;
        //std::cout<< runTime_.timeName() << " : Data has been exported from OpenFOAM to CoSimulation: Force values with array size = " << Nodal_Force_data_to_send.size() << std::endl;
    }

    return true;
}

bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::end()
{
    //Dicsonect from CoSimIO
    //CoSimulationAdapter_.end();

    std::cout << "\n" << "CoSimulation Adapter's function object : end()" << std::endl;

    //DisConnection between openFOAM and Kratos-CoSimulation using CoSimIO
    CoSimIO::Info connect_info;
    CoSimIO::Info disconnect_settings;
    disconnect_settings.Set("connection_name", connection_name);
    connect_info = CoSimIO::Disconnect(disconnect_settings); // disconnect afterwards
    COSIMIO_CHECK_EQUAL(connect_info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Disconnected);

    return true;
}

bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::write()
{
    //adapter_.write();
    //cout<<"Welcome to the Kratos_OpenFOAM_Adapter\n";

    return true;
}

// *********************************************** Some Auxillary Functions **************************************************//
//Calculate the Solver Type - according to the pressure dimension
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
            std::cout << "Solver Type : Compressible " << std::endl;
        }
        else if (p_.dimensions() == pressureDimensionsIncompressible)
        {
            solverType_ = "incompressible";
            std::cout << "Solver Type : InCompressible " << std::endl;
        }
    }

    if(solverType_  == "unknown")
    {
        std::cout << "Solver Type: Neither Compressible nor Incompresible" << std::endl;
    }

    return solverType_;
}

// To compare the Foam::Vector . Check whether it is really required?
bool is_same_points(Foam::vector& pointX, Foam::vector& pointY)// Working perfectly . DO NOT CHECK AGAIN
{
    if(pointX[0] == pointY[0] && pointX[1] == pointY[1] && pointX[2] == pointY[2])
        return true;
    else
        return false;
}

//To Compare the new node with all previous nodes before creating the new CoSimIO::Node
int Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::compare_nodes(Foam::vector& pointX)
{
    int answer = 0;

    for(Foam::vector& nodei : array_of_nodes)
    {
        if(is_same_points(pointX , nodei))//current node == previous all nodes
        {
            auto itr = std::find( array_of_nodes.begin(), array_of_nodes.end(), nodei ) ;
            answer = std::distance(array_of_nodes.begin(), itr) + 1 ; //nodeindex starts from 1 in CoSimIO
            break; //Once repeated node is found immediate break it. No need to check it later
        }
        else
        {
            answer = -1; // Compare with all the available nodes and then return (-1) if pointX is not found. Then new node is created in the main()
        }
    }

    return answer;
}

// ************************** Force / load calculation ************************************//
//Calculate viscous Force
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

//Finding correct rho
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

//Finding correct mu
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

//Normal vectors multiplied by face area
Foam::vectorField Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::getFaceVectors(const unsigned int patchID) const
{
    // Normal vectors multiplied by face area
    return mesh_.boundary()[patchID].Sf();
}

//Calculate Total Force
bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::calculateForces(std::size_t interface_index)
{
    std::vector<int> patchIDs;
    // For every patch that participates in the coupling interface
    for (uint i = 0; i < interfaces_.at(interface_index).patchNames.size(); i++)
    {
        // Get the patchID
        int patchID = mesh_.boundaryMesh().findPatchID((interfaces_.at(interface_index).patchNames).at(i));

        std::cout<< "Patch for force " << patchID << std::endl;

        // Throw an error if the patch was not found
        if (patchID == -1)
        {
            std::cout<< "ERROR: Patch " << (interfaces_.at(interface_index).patchNames).at(i) << " does not exist." << std::endl;
        }

    // Add the patch in the list
    patchIDs.push_back(patchID);
    }

    //Get different force fields from OpenFOAM, See Force Function Object
    //1. Normal vectors on the boundary, multiplied with the face areas
    /* const surfaceVectorField::Boundary& Sfb
    (
        mesh_.Sf().boundaryField()
    ); */

    //2. Stress tensor boundary field
    tmp<volSymmTensorField> tdevRhoReff(devRhoReff());
    const volSymmTensorField::Boundary& devRhoReffb
    (
        tdevRhoReff().boundaryField()
    );

    //3. Density boundary field
    tmp<volScalarField> trho(rho());
    const volScalarField::Boundary& rhob = trho().boundaryField();

    //4. Pressure boundary field
    tmp<volScalarField> tp = mesh_.lookupObject<volScalarField>("p");
    const volScalarField::Boundary& pb
    (
        tp().boundaryField()
    );



    // For every boundary patch of the interface
    for(uint j=0; j< patchIDs.size(); j++)
    {
        int patchID = patchIDs.at(j);

        //printing pressure values
        /* std::cout << "Pressure values: "<< std::endl;
        forAll(pb[patchID],i)
        {
            std::cout << "i= " << i << "Pressure = " << pb(i)->evaluate() << std::endl;
        } */

        const auto& surface = getFaceVectors(patchID); //newly modified on 2.6.2021

        //Pressure forces
        if(solverType_.compare("incompressible") == 0)
        {
           // Force_->boundaryFieldRef()[patchID] = Sfb[patchID] * pb[patchID] * rhob[patchID];
           Force_->boundaryFieldRef()[patchID] = surface * pb[patchID] * rhob[patchID];
        }
        else if(solverType_.compare("compressible") == 0)
        {
            //Force_->boundaryFieldRef()[patchID] = Sfb[patchID] * pb[patchID];
            Force_->boundaryFieldRef()[patchID] = surface * pb[patchID];
        }
        else
        {
            FatalErrorInFunction << "Forces calculation does only support compressible or incompressible solver type." << exit(FatalError);
        }

        //Viscous forces
        //Force_->boundaryFieldRef()[patchID] += Sfb[patchID] & devRhoReffb[patchID];
        Force_->boundaryFieldRef()[patchID] += surface & devRhoReffb[patchID];

        // Now write this forces into Buffer to send further to the Strctural Solver
        int bufferIndex = 0;
        // For every cell of the patch
        forAll(Force_->boundaryField()[patchID], i)
        {
            // Copy the force into the buffer
            // x-dimension
            data_to_send[bufferIndex++] = Force_->boundaryField()[patchID][i].x();


            // y-dimension
            data_to_send[bufferIndex++] = Force_->boundaryField()[patchID][i].y();

            if(dim == 3)
            {
                // z-dimension
                //data_to_send[bufferIndex++] = Force_->boundaryField()[patchID][i].z();
                data_to_send[bufferIndex++] = 0.0; //We will skip 3rd dimension
            }
        }

    }

    /* int i=0;
    for(auto& value : data_to_send)
    {
        std::cout << "id = " << i << ", Value = " << value << std::endl;
        i++;
    } */

    return true;
}

}//namespace Foam

// ************************************************************************* //
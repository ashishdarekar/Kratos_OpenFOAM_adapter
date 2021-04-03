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
fvMeshFunctionObject(name, runTime, dict)//, CoSimulationAdapter_()
{
    read(dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::~KratosOpenfoamAdapterFunctionObject()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//const pointField& p = points();

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

    //Check the total number of registered objects in the PolyMesh Object Registry related to Given Solver
    std::cout<< "Name of all registered objects in Foam::PolyMesh object Registry are:" << std::endl;
    Foam::wordList Objectnames_ = mesh_.names(); //List of all Objects in the polymesh class::mesh_object
    forAll(Objectnames_,i)
    {
        std::cout << Objectnames_[i] << ", ";
    }
    std::cout<<"\n";

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

    //Connection between OpneFOAM and Kratos-CoSimulation using CoSimIO
    CoSimIO::Info settings;
    settings.Set("my_name", "Openfoam_adapter");
    settings.Set("connect_to", "Kratos_CoSimulation");
    settings.Set("echo_level", 1);
    settings.Set("version", "1.25");

    auto connect_info = CoSimIO::Connect(settings);
    COSIMIO_CHECK_EQUAL(connect_info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Connected);
    connection_name = connect_info.Get<std::string>("connection_name");

    //*********************** Export Fluid-OpenFOAM Mesh to CoSimulation using CoSimIO *********************//
    //Create one CoSimIO::ModelPart for each coupling interface
    std::cout << "Exporting Mesh: Start" << std::endl;

    for(std::size_t j=0; j < num_interfaces_; j++)
    {
        std::string interface_name = "interface" + std::to_string(j+1);

        //model_part_interfaces_.push_back(CoSimIO::ModelPart(interface_name));
        model_part_interfaces_.push_back(CoSimIO::make_unique<CoSimIO::ModelPart>(interface_name));

        std::cout<< "name of the model part: " << model_part_interfaces_.at(j)->Name() << std::endl;

        // **************************Create a mesh as a ModelPart********************************//
        std::cout << "Creating Mesh as a ModelPart: Start" << std::endl;
        std::vector<int> patchIDs;
        uint numNodes = 0;
        uint numElements = 0;

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
            numNodes += mesh_.boundaryMesh()[patchIDs.at(i)].localPoints().size();
        }
        std::cout << "Number of Nodes: " << numNodes << std::endl;

        // Count the number of elements/faces for all the patches in that interface
        for (uint i = 0; i < patchIDs.size(); i++)
        {
            numElements += mesh_.boundary()[patchIDs[i]].size();
        }
        std::cout << "Number of Elements: " << numElements << std::endl;

        //-For Nodes and Element IDs for CoSimIO
        std::vector<int> NodeIDs;
        NodeIDs.resize(numNodes);
        int nodeIndex = 1; //As Node indexing starts with 1 in CoSimIO

        std::vector<int> ElemIDs;
        ElemIDs.resize(numElements);
        int elemIndex = 1; //As element indexing starts with 1 in CoSimIO

        vector pointX(0,0,0); //For accessing the Co-ordinates of Nodes

        for(uint i = 0; i < patchIDs.size(); i++)
        {
            const word& patchName = mesh_.boundary()[patchIDs[i]].name();
            std::cout << "patchID: " << patchIDs[i] << " with name "  << patchName << " With size = "<< mesh_.boundary()[patchIDs[i]].size() << std::endl;

            forAll (mesh_.boundary()[patchIDs[i]],facei)
            {
                const label& faceID = mesh_.boundaryMesh()[patchIDs[i]].start() + facei;
                std::cout << "ID of this face: " << faceID << " Which face: "<< facei << std::endl;
                std::vector<CoSimIO::IdType> connectivity;
                forAll (mesh_.faces()[faceID], nodei)
                {
                    const label& nodeID = mesh_.faces()[faceID][nodei]; //for OpenFOAM
                    pointX = mesh_.points()[nodeID];
                    std::cout << "Node "<< nodei << " with Id=" << nodeID << ":" << pointX[0] << "," << pointX[1] << "," << pointX[2] << std::endl;
                    //-Make CoSimIO Nodes
                    NodeIDs.push_back(nodeIndex); //Later used to make CoSimIO::Element
                    CoSimIO::Node& node = model_part_interfaces_.at(j)->CreateNewNode( nodeIndex, pointX[0], pointX[1], pointX[2]);
                    connectivity.push_back(nodeIndex); //connectivity to make that element /face
                    nodeIndex++;
                }
                //-Make CoSimIO Element = currently hardcoded to make Quadrilateral3D4
                ElemIDs.push_back(elemIndex); //Maybe useful
                CoSimIO::Element& element = model_part_interfaces_.at(j)->CreateNewElement( elemIndex, CoSimIO::ElementType::Quadrilateral3D4, connectivity );
                elemIndex++;
            }
        }
        std::cout << "Creating Mesh as a ModelPart: Done" << std::endl;
        // **************************Done: Create a mesh as a ModelPart********************************//

        //Import InterfaceMesh/ModelPart to CoSimulation using CoSimIO
        CoSimIO::Info info;
        info.Clear();
        info.Set("identifier", "fluid_mesh");
        info.Set("connection_name", connection_name);
        auto export_info = CoSimIO::ExportMesh(info, *model_part_interfaces_.at(j));
        std::cout << "Exporting Mesh as a ModelPart: Done for interface = " << j+1 << std::endl;
    }

    std::cout << "Exporting Mesh: Done" << std::endl;

    return true;
}


bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::execute()
{
    //CoSimulationAdapter_.execute();

    //Currently, exporting Data on 3rd timestep
    if(time_step_ == 3)
    {
        std::cout << "CoSimulation Adapter's function object : execution()" << std::endl;

        //pressure Data values in the simulation
        if(mesh_.foundObject<volScalarField>("p"))
        {
            std::cout<< "Pressure field is found" << std::endl;

            const volScalarField& pressure_ = mesh_.lookupObject<volScalarField>("p");
            std::cout<< "Size of the array is " << pressure_.size() << std::endl;

            std::vector<double> data_to_send;
            data_to_send.resize(pressure_.size());

            forAll(pressure_, i)
            {
                //std::cout << pressure_[i] << std::endl;
                data_to_send[i] = pressure_[i];
            }

            //Export this pressure arrat to CoSimulation
            CoSimIO::Info connect_info;
            connect_info.Clear();
            connect_info.Set("identifier", "pressure_values");
            connect_info.Set("connection_name", connection_name);
            connect_info = CoSimIO::ExportData(connect_info, data_to_send);

            std::cout<< "Data has been exported from OpenFOAM to CoSimulation" <<std::endl;

        }

        /*//Velocity Data values in the simulation
         if(mesh_.foundObject<volVectorField>("U"))
        {
            std::cout<< "Velocity field is found" << std::endl;

            const volVectorField& velocity_ = mesh_.lookupObject<volVectorField>("U");
            std::cout<< "Size of the array is " << velocity_.size() << std::endl;
            forAll(velocity_, i)
            {
                std::cout << "(" << velocity_[i][1] <<  "," << velocity_[i][2] << "," << velocity_[i][3] << ")" << std::endl;
            }
        }*/
    }

    time_step_++;

    return true;
}


bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::end()
{
    //Dicsonect from CoSimIO
    //CoSimulationAdapter_.end();

    std::cout << "CoSimulation Adapter's function object : end()" << std::endl;

    //DisConnection between OpneFOAM and Kratos-CoSimulation using CoSimIO
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

}//namespace Foam

// ************************************************************************* //
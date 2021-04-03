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
        //int numCells = 0;

        //******************Information about all cells/faces in the boundary patches *************************//
        const pointField& p = mesh_.points();
        const labelList& cellPts = mesh_.boundaryMesh()[1].localFaces()[0]; //for 1st cell
        auto numfaces = mesh_.boundaryMesh()[1].localFaces();
        std::cout<< typeid(numfaces).name()<< ", " << mesh_.boundaryMesh()[1].meshPoints().size() << "," << mesh_.boundaryMesh()[1].localFaces().size() << std::endl;
        std::cout<< "Indexes of the nodes for face/cell-0 in this patch are :( " << numfaces[0][0] << "," << numfaces[0][1] <<"," << numfaces[0][2]<< "," << numfaces[0][3] <<" )"  << std::endl;
        std::cout<< "Indexes of the nodes for face/cell-1 in this patch are :( " << numfaces[1][0] << "," << numfaces[1][1] <<"," << numfaces[1][2]<< "," << numfaces[1][3] <<" )"  << std::endl;
        std::cout<< "Indexes of the nodes for face/cell-2 in this patch are :( " << numfaces[2][0] << "," << numfaces[2][1] <<"," << numfaces[2][2]<< "," << numfaces[2][3] <<" )"  << std::endl;
        std::cout<< "Indexes of the nodes for face/cell-3 in this patch are :( " << numfaces[3][0] << "," << numfaces[3][1] <<"," << numfaces[3][2]<< "," << numfaces[3][3] <<" )"  << std::endl;
        std::cout<< "Indexes of the nodes for face/cell-8 in this patch are :( " << numfaces[8][0] << "," << numfaces[8][1] <<"," << numfaces[8][2]<< "," << numfaces[8][3] <<" )"  << std::endl;
        std::cout<< "Indexes of the nodes for face/cell-13 in this patch are :( " << numfaces[13][0] << "," << numfaces[13][1] <<"," << numfaces[13][2]<< "," << numfaces[13][3] <<" )"  << std::endl;
        std::cout<< "Indexes of the nodes for face/cell-23 in this patch are :( " << numfaces[23][0] << "," << numfaces[23][1] <<"," << numfaces[23][2]<< "," << numfaces[23][3] <<" )"  << std::endl;
        std::cout<< "Indexes of the nodes for face/cell-69 in this patch are :( " << numfaces[69][0] << "," << numfaces[69][1] <<"," << numfaces[69][2]<< "," << numfaces[69][3] <<" )"  << std::endl;


        std::cout<< "Coordinates of the node-0 of the 1st face/cell are :( " << p[numfaces[0][0]][0] << "," << p[numfaces[0][0]][1] << "," << p[numfaces[0][0]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-1 ot the 1st face/cell are :( " << p[numfaces[1][0]][0] << "," << p[numfaces[1][0]][1] << "," << p[numfaces[1][0]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-2 ot the 1st face/cell are :( " << p[numfaces[2][0]][0] << "," << p[numfaces[2][0]][1] << "," << p[numfaces[2][0]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-3 ot the 1st face/cell are :( " << p[numfaces[3][0]][0] << "," << p[numfaces[3][0]][1] << "," << p[numfaces[3][0]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-0 of the 1st face/cell are :( " << p[numfaces[4][0]][0] << "," << p[numfaces[4][0]][1] << "," << p[numfaces[4][0]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-1 ot the 1st face/cell are :( " << p[numfaces[5][0]][0] << "," << p[numfaces[5][0]][1] << "," << p[numfaces[5][0]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-2 ot the 1st face/cell are :( " << p[numfaces[6][0]][0] << "," << p[numfaces[6][0]][1] << "," << p[numfaces[6][0]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-3 ot the 1st face/cell are :( " << p[numfaces[7][0]][0] << "," << p[numfaces[7][0]][1] << "," << p[numfaces[7][0]][2] << " )" << std::endl;


        std::cout<< "Coordinates of the node-0 of the 1st face/cell are :( " << p[cellPts[0]][0] << "," << p[cellPts[0]][1] << "," << p[cellPts[0]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-1 ot the 1st face/cell are :( " << p[cellPts[1]][0] << "," << p[cellPts[1]][1] << "," << p[cellPts[1]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-2 ot the 1st face/cell are :( " << p[cellPts[2]][0] << "," << p[cellPts[2]][1] << "," << p[cellPts[2]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-3 ot the 1st face/cell are :( " << p[cellPts[3]][0] << "," << p[cellPts[3]][1] << "," << p[cellPts[3]][2] << " )" << std::endl;


        std::cout<< "Coordinates of the node-0 are :( " << mesh_.points()[0][0] << "," << mesh_.points()[0][1] << "," << mesh_.points()[0][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-1 are :( " << mesh_.points()[1][0] << "," << mesh_.points()[1][1] << "," << mesh_.points()[1][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-2 are :( " << mesh_.points()[2][0] << "," << mesh_.points()[2][1] << "," << mesh_.points()[2][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-3 are :( " << mesh_.points()[3][0] << "," << mesh_.points()[3][1] << "," << mesh_.points()[3][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-4 are :( " << mesh_.points()[4][0] << "," << mesh_.points()[4][1] << "," << mesh_.points()[4][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-5 are :( " << mesh_.points()[5][0] << "," << mesh_.points()[5][1] << "," << mesh_.points()[5][2] << " )" << std::endl;

        std::cout<< "Coordinates of the node-0 are :( " << mesh_.boundaryMesh()[1].points()[0][0] << "," << mesh_.boundaryMesh()[1].points()[0][1] << "," << mesh_.boundaryMesh()[1].points()[0][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-1 are :( " << mesh_.boundaryMesh()[1].points()[1][0] << "," << mesh_.boundaryMesh()[1].points()[1][1] << "," << mesh_.boundaryMesh()[1].points()[1][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-2 are :( " << mesh_.boundaryMesh()[1].points()[2][0] << "," << mesh_.boundaryMesh()[1].points()[2][1] << "," << mesh_.boundaryMesh()[1].points()[2][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-3 are :( " << mesh_.boundaryMesh()[1].points()[3][0] << "," << mesh_.boundaryMesh()[1].points()[3][1] << "," << mesh_.boundaryMesh()[1].points()[3][2] << " )" << std::endl;

        vector pointX(0,0,0);

        forAll(mesh_.boundary(), patchID)
        {
            const word& patchName = mesh_.boundary()[patchID].name();
            std::cout << "patchID: " << patchID << " with name "  << patchName << " With size = "<< mesh_.boundary()[patchID].size() << std::endl;

            forAll (mesh_.boundary()[patchID],facei)
            {
                const label& faceID = mesh_.boundaryMesh()[patchID].start() + facei;
                std::cout << "ID of this face: " << faceID << " Which face: "<< facei << std::endl;
                forAll (mesh_.faces()[faceID], nodei)
                {
                    const label& nodeID = mesh_.faces()[faceID][nodei];
                    pointX = mesh_.points()[nodeID];
                    std::cout << "Node "<< nodei << " with Id=" << nodeID << ":" << pointX[0] << "," << pointX[1] << "," << pointX[2] << std::endl;
                }
            }
        }


        //******************Information about all cells in the domain/Internal Mesh *************************//
        /* const pointField& p = mesh_.points();
        auto numCells = mesh_.cellPoints()[0];
        const labelList& cellPts = mesh_.cellPoints()[0];

        std::cout<< typeid(numCells).name()<< ", " << mesh_.cellPoints().size() << std::endl;
        std::cout<< "Indexes of the nodes for cell-0 are :( " << numCells[0] << "," << numCells[1] << "," << numCells[2] << "," << numCells[3] << "," << numCells[4] << "," << numCells[5] << "," << numCells[6] << "," << numCells[7] << " )" << std::endl;
        std::cout<< "Coordinates of the node-0 of that cell are :( " << p[cellPts[0]][0] << "," << p[cellPts[0]][1] << "," << p[cellPts[0]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-1 ot that cell are :( " << p[cellPts[1]][0] << "," << p[cellPts[1]][1] << "," << p[cellPts[1]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-2 ot that cell are :( " << p[cellPts[2]][0] << "," << p[cellPts[2]][1] << "," << p[cellPts[2]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-3 ot that cell are :( " << p[cellPts[3]][0] << "," << p[cellPts[3]][1] << "," << p[cellPts[3]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-4 ot that cell are :( " << p[cellPts[4]][0] << "," << p[cellPts[4]][1] << "," << p[cellPts[4]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-5 ot that cell are :( " << p[cellPts[5]][0] << "," << p[cellPts[5]][1] << "," << p[cellPts[5]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-4 ot that cell are :( " << p[cellPts[6]][0] << "," << p[cellPts[6]][1] << "," << p[cellPts[6]][2] << " )" << std::endl;
        std::cout<< "Coordinates of the node-5 ot that cell are :( " << p[cellPts[7]][0] << "," << p[cellPts[7]][1] << "," << p[cellPts[7]][2] << " )" << std::endl;
 */


        // For every patch that participates in the coupling
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

        // Count the Nodes for all the patches
        for (uint i = 0; i < patchIDs.size(); i++)
        {
            numNodes += mesh_.boundaryMesh()[patchIDs.at(i)].localPoints().size();
        }
        std::cout << "Number of Nodes: " << numNodes << std::endl;

        // Array of the indices of the Nodes. Each node has one index.
        int * NodeIDs;
        NodeIDs = new int[numNodes];

        // Get the locations of the mesh vertices, for all the patches
        for (uint i = 0; i < patchIDs.size(); i++)
        {
            // Get the face centers of the current patch
            const pointField Nodes = mesh_.boundaryMesh()[patchIDs.at(i)].localPoints();

            //-Make CoSimIO Nodes
            for(uint k = 0; k< numNodes; k++)
            {
                NodeIDs[k] = k+1; //CoSimIO Starts with 1 while OpenFOAM starts with 0
                CoSimIO::Node& node = model_part_interfaces_.at(j)->CreateNewNode( NodeIDs[k], Nodes[k].x(), Nodes[k].y(), Nodes[k].z());
            }
        }

        // Count the Elements for all the patches
        for (uint i = 0; i < patchIDs.size(); i++)
        {
            numElements += mesh_.boundaryMesh()[patchIDs.at(i)].faceCentres().size();
        }
        std::cout << "Number of Elements: " << numElements << std::endl;

        //-Make CoSimIO Elements
        int * elem_IDs;
        elem_IDs = new int[numElements];
        for(uint k = 0; k< numElements; k++)
        {
            elem_IDs[k] = k+1;//CoSimIO Starts with 1 while OpenFOAM starts with 0
            std::vector<CoSimIO::IdType> connectivity {NodeIDs[k], NodeIDs[k+1], NodeIDs[k+2], NodeIDs[k+3] }; //Ids of the nodes
            CoSimIO::Element& element = model_part_interfaces_.at(j)->CreateNewElement( elem_IDs[k], CoSimIO::ElementType::Quadrilateral3D4, connectivity );
        }

        std::cout << "Creating Mesh as a ModelPart: Done" << std::endl;
        // **************************Done: Create a mesh as a ModelPart********************************//

        //Import mesh to Cosimulation using CoSimIO
        CoSimIO::Info info;
        info.Clear();
        info.Set("identifier", "fluid_mesh");
        info.Set("connection_name", connection_name);
        auto export_info = CoSimIO::ExportMesh(info, *model_part_interfaces_.at(j));
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
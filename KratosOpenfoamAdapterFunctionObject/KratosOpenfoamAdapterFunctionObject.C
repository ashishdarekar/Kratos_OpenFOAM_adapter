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
    std::cout << "CoSimulation Adapter's function object : read()" << std::endl;

    // Reading configuration parameters from ControlDict related to function object
    fvMeshFunctionObject::read(dict);

    my_name = dict.lookupOrDefault<word>("participant", "fluid");
    std::cout<< "Name of the participant is: " << my_name <<std::endl;

    dim = dict.lookupOrDefault<double>("dim", 3);
    std::cout<< "Dimension of a problem is: " << dim <<std::endl;

    thick = dict.lookupOrDefault<double>("thick", 1.0);
    std::cout<< "Thickness of a domain is: " << thick <<std::endl;

    // Check the solver type and determine it if needed
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

    // Every interface is a subdictionary of "interfaces"
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

                wordList exportData = interfaceSubdict.lookupType<wordList>("exportData");
                for(auto wData : exportData)
                {
                    interfacedata.exportData.push_back(wData);
                }

                //Add this interface in the Array of all interfaces
                interfaces_.push_back(interfacedata);

                num_interfaces_++;
            }
        }
    }

    std::cout << "Interfaces Reading: Done" << std::endl;
    std::cout << "Number of interfaces found: " << num_interfaces_<< std::endl;

    // Connection between openFOAM and Kratos-CoSimulation using CoSimIO
    CoSimIO::Info settings;
    settings.Set("my_name", "Openfoam_Adapter");
    settings.Set("connect_to", "Openfoam_Kratos_Wrapper");
    settings.Set("echo_level", 0);
    settings.Set("version", "1.25");

    auto connect_info = CoSimIO::Connect(settings);
    COSIMIO_CHECK_EQUAL(connect_info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Connected);
    connection_name = connect_info.Get<std::string>("connection_name");

    // To test the CoSimulation Import Mesh Fuctationality with one interface only
    std::cout << "\n" <<"********************** Exporting InterfaceMesh using ModelPart: Start **********************" << "\n" <<std::endl;
    for(std::size_t j = 0; j < num_interfaces_; j++)
    {
        std::string interface_name = "interface" + std::to_string(j+1);

        // *******************************Create a mesh as a ModelPart************************************ //
        std::cout << "Accessing Mesh from openFOAM" << std::endl;
        std::vector<int> patchIDs;

        // For every patch that participates in the coupling interface. We are keeping one patch for one interface
        for (std::size_t i = 0; i < interfaces_.at(j).patchNames.size(); i++)
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
        for (std::size_t i = 0; i < patchIDs.size(); i++)
        {
            interfaces_.at(j).numNodes += mesh_.boundaryMesh()[patchIDs.at(i)].localPoints().size();
        }
        std::cout << "Total Number of Nodes in this interface: " << interfaces_.at(j).numNodes  << std::endl;

        // Count the number of elements/faces for all the patches in that interface
        for (std::size_t i = 0; i < patchIDs.size(); i++)
        {
            interfaces_.at(j).numElements += mesh_.boundary()[patchIDs[i]].size();
        }
        std::cout << "Total Number of Elements/faces in this interface: " << interfaces_.at(j).numElements << std::endl;

        // Considering 1st interface ONLY
        // Create CoSimIO::ModelPart- Put in public Member function
        CoSimIO::ModelPart model_part_interface_flap("interface_flap_model_part");

        // For Nodes and Element IDs for CoSimIO
        std::cout << "Creating Model Part (Nodes and Elements) for CoSimIO : start" << std::endl;
        std::vector<int> NodeIDs;
        NodeIDs.resize( interfaces_.at(j).numNodes );
        int nodeIndex = 1; //As Node indexing starts with 1 in CoSimIO

        std::vector<int> ElemIDs;
        ElemIDs.resize(interfaces_.at(j).numElements);
        int elemIndex = 1; //As element indexing starts with 1 in CoSimIO

        // Accessing the Co-ordinates of nodes in the Inteface and making CoSimIO nodes and elements
        for(std::size_t i = 0; i < patchIDs.size(); i++)
        {
            forAll(mesh_.boundary()[patchIDs[i]],facei)
            {
                const label& faceID = mesh_.boundaryMesh()[patchIDs[i]].start() + facei;

                std::vector<CoSimIO::IdType> connectivity;
                forAll(mesh_.faces()[faceID], nodei)
                {
                    const label& nodeID = mesh_.faces()[faceID][nodei]; //for OpenFOAM
                    auto pointX = mesh_.points()[nodeID];

                    int result = compare_nodes(pointX); // return nodeIndex if node is already present and (-1) if node is not present
                    if(result == (-1)) // For new node
                    {
                        // Make CoSimIO Nodes
                        NodeIDs.push_back(nodeIndex); // Later used to make CoSimIO::Element
                        CoSimIO::Node& node = model_part_interface_flap.CreateNewNode( nodeIndex, pointX[0], pointX[1], pointX[2]);

                        array_of_nodes.push_back(pointX);// Push new element in the list to compare
                        connectivity.push_back(nodeIndex); // connectivity to make that element
                        nodeIndex++;
                    }
                    else // Old node index just push in connectivity to make new element
                    {
                        connectivity.push_back(result); // connectivity to make that element
                    }
                }

                ElemIDs.push_back(elemIndex); // For future use
                CoSimIO::Element& element = model_part_interface_flap.CreateNewElement( elemIndex, CoSimIO::ElementType::Quadrilateral2D4, connectivity );
                elemIndex++;
            }
        }
        std::cout << "Converting InterfaceMesh to a CoSimIO::ModelPart -> End" << std::endl;

        // Export InterfaceMesh/ModelPart to CoSimulation using CoSimIO
        std::cout << "Exporting Mesh as a ModelPart: Start for an interface flap" << std::endl;
        CoSimIO::Info info;
        info.Clear();
        info.Set("identifier", "interface_flap");
        info.Set("connection_name", connection_name);
        auto export_info = CoSimIO::ExportMesh(info, model_part_interface_flap);
        std::cout << "*********************** Exporting InterfaceMesh using ModelPart: End ************************" << "\n" <<std::endl;

    }


    // Resizing the Data Vectors for Import Export operations
    for(std::size_t i=0; i < num_interfaces_; i++)
    {
        // Import Data from CoSimulation (Displacement/Delta) present only on faceNodes/Nodes. No need to resize it
        for(std::size_t j=0; j< interfaces_.at(i).importData.size(); j++)
        {
            std::string dataName = interfaces_.at(i).importData.at(j);

            if(dataName.find("Displacement") == 0 || dataName.find("DisplacementDelta") == 0) //If "Displacement" or "DisplacementDelta" string is found it will return 0
            {
                //data_to_recv.resize((interfaces_.at(i).numNodes) * dim);
            }
            //else if() //if some other variables
        }

        // Export Data to CoSimulation (force/Stress) present only on faceCenters/Elements
        for(std::size_t j=0; j< interfaces_.at(i).exportData.size(); j++)
        {
            std::string dataName = interfaces_.at(i).exportData.at(j);

            if(dataName.find("Force") == 0 || dataName.find("Stress") == 0) //If "force" or "stress" string is found it will return 0
            {
                data_to_send.resize((interfaces_.at(i).numElements) * dim);
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
        // Import the displacement array from the CoSimulation
        CoSimIO::Info connect_info;
        connect_info.Clear();
        connect_info.Set("identifier", "disp");
        connect_info.Set("connection_name", connection_name);
        connect_info = CoSimIO::ImportData(connect_info, data_to_recv);
        //Check size of Receive data = Expected Receive data. Get it from top.
        //COSIMIO_CHECK_EQUAL(data_to_recv.size(), );

        std::cout<< runTime_.timeName() << " : Data has been imported from CoSimulation to OpenFOAM: Disp values with array size = " << data_to_recv.size() << std::endl;

        // Get the displacement on the patch and assign it those values received from CoSimulation,
        // For every patch that participates in the coupling interface
        std::cout<< "Displacement replacement : start" << std::endl;
        for (std::size_t i = 0; i < interfaces_.size(); i++)
        {
            for (std::size_t j = 0; j < interfaces_.at(i).patchNames.size(); j++)
            {
                Foam::pointVectorField* point_disp;
                point_disp = const_cast<pointVectorField*>( &mesh_.lookupObject<pointVectorField>("pointDisplacement") );
                label patchIndex = mesh_.boundaryMesh().findPatchID(interfaces_.at(i).patchNames[j]);//Remove hardcoded part for finding patchIndex
                std::cout << "Interface number = " << i << " and patch number = " << patchIndex << " and patch name = " << interfaces_.at(i).patchNames[j] << std::endl;
                fixedValuePointPatchVectorField& pointDisplacementFluidPatch = refCast<fixedValuePointPatchVectorField>(point_disp->boundaryFieldRef()[patchIndex]);

                int iterator = 0;
                forAll(point_disp->boundaryFieldRef()[patchIndex] ,i)
                {
                    pointDisplacementFluidPatch[i][0] = data_to_recv[iterator++];
                    pointDisplacementFluidPatch[i][1] = data_to_recv[iterator++];
                    if (dim ==3)
                        pointDisplacementFluidPatch[i][2] = data_to_recv[iterator++];
                }
            }
        }
        std::cout<< "Displacement replacement : End" << std::endl;

        // *************************************** Force/Load Related ****************************************** //
        std::cout<< "Force calculation : start" << std::endl;
        for(std::size_t i=0; i < num_interfaces_; i++)
        {
            // For "Wirte Data" variables which need to send to CoSimulation
            for(std::size_t j=0; j< interfaces_.at(i).exportData.size(); j++)
            {
                std::string dataName = interfaces_.at(i).exportData.at(j);
                std::cout<< "interface number = "<< i << " with Export DataName = " << dataName << std::endl;

                if(dataName.find("Force") == 0 )
                {
                    calculateForces(j);
                }
            }
        }
        std::cout<< "Force calculation : End" << std::endl;

        // Export this force array to CoSimulation //Elemental Force Data
        connect_info.Clear();
        connect_info.Set("identifier", "load_cells");
        connect_info.Set("connection_name", connection_name);
        connect_info = CoSimIO::ExportData(connect_info, data_to_send);

        std::cout<< runTime_.timeName() << " : Data has been exported from OpenFOAM to CoSimulation: Force values with array size = " << data_to_send.size() << std::endl;
    }

    return true;
}

bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::end()
{
    // Dicsonect from CoSimIO
    std::cout << "\n" << "CoSimulation Adapter's function object : end()" << std::endl;

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
    std::cout << "\n" << "CoSimulation Adapter's function object : write()" << std::endl;

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

// To compare the Foam::Vector. Check whether it is really required?
bool is_same_points(Foam::vector& pointX, Foam::vector& pointY)
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

        std::cout<< "PatchNumber for force calculation " << patchID << std::endl;

        // Throw an error if the patch was not found
        if (patchID == -1){
            std::cout<< "ERROR: Patch " << (interfaces_.at(interface_index).patchNames).at(i) << " does not exist." << std::endl;
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
            data_to_send[bufferIndex++] = Force_->boundaryField()[patchID][i].x();

            // y-dimension
            data_to_send[bufferIndex++] = Force_->boundaryField()[patchID][i].y();

            if(dim == 3)
            {
                // z-dimension
                data_to_send[bufferIndex++] = Force_->boundaryField()[patchID][i].z();
            }
        }

    }

    return true;
}

}//namespace Foam

// ************************************************************************* //
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

    /* //Check the total number of registered objects in the PolyMesh Object Registry related to Given Solver
    std::cout<< "Name of all registered objects in Foam::PolyMesh object Registry are:" << std::endl;
    Foam::wordList Objectnames_ = mesh_.names(); //List of all Objects in the polymesh class::mesh_object
    forAll(Objectnames_,i)
    {
        std::cout << Objectnames_[i] << ", ";
    }
    std::cout<<"\n"; */

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
    //Create one CoSimIO::ModelPart for each coupling interface
    std::cout << "\n" <<"********************** Exporting InterfaceMesh using ModelPart: Start **********************" << "\n" <<std::endl;

    for(std::size_t j=0; j < num_interfaces_; j++)
    {
        std::string interface_name = "interface" + std::to_string(j+1);

        //model_part_interfaces_.push_back(CoSimIO::ModelPart(interface_name));
        model_part_interfaces_.push_back(CoSimIO::make_unique<CoSimIO::ModelPart>(interface_name));

        std::cout<< "Name of an Interface/model part: " << model_part_interfaces_.at(j)->Name() << std::endl;

        // **************************Create a mesh as a ModelPart********************************//
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
                CoSimIO::Element& element = model_part_interfaces_.at(j)->CreateNewElement( elemIndex, CoSimIO::ElementType::Quadrilateral3D4, connectivity );
                elemIndex++;
            }
        }
        std::cout << "Converting InterfaceMesh to a CoSimIO::ModelPart -> End" << std::endl;
        // **************************Done: Create a mesh as a ModelPart********************************//

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
                }
            }
            //else if() //if some other variables
        }

        //For "Read Data" variables which need to receive from CoSimulation
        for(std::size_t j=0; j< interfaces_.at(i).readData.size(); j++)
        {
            std::string dataName = interfaces_.at(i).readData.at(j);

            if(dataName.find("Displacement") == 0 || dataName.find("DisplacementDelta") == 0) //If "Displacement" or "DisplacementDelta" string is found it will return 0
            {
                if(interfaces_.at(i).locationsType == "faceNodes")
                {
                    data_to_recv.resize((interfaces_.at(i).numNodes) * dim);
                }
                else if(interfaces_.at(i).locationsType == "faceCenters")
                {
                    data_to_recv.resize((interfaces_.at(i).numElements) * dim);
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
        connect_info.Set("identifier", "disp_values");
        connect_info.Set("connection_name", connection_name);
        connect_info = CoSimIO::ImportData(connect_info, data_to_recv);

        std::cout<< runTime_.timeName() << " : Data has been imported from CoSimulation to OpenFOAM: Disp values with array size = " << data_to_recv.size() << std::endl;

        /* // Get the displacement on the patch and assign it those values received from CoSimulation,
        // currently values are random hence if it may break the solution.
        Foam::pointVectorField* point_disp;
        point_disp = const_cast<pointVectorField*>( &mesh_.lookupObject<pointVectorField>("pointDisplacement") );
        //std::cout<< "Size of the pointDisplacement array is " << point_disp->size() << std::endl;
        label patchIndex = mesh_.boundaryMesh().findPatchID("flap");
        fixedValuePointPatchVectorField& pointDisplacementFluidPatch = refCast<fixedValuePointPatchVectorField>(point_disp->boundaryFieldRef()[patchIndex]);

        int iterator = 0 ; //iterator goes from 0 till (size of recv_array)
        forAll(point_disp->boundaryFieldRef()[patchIndex] ,i)
        {
            pointDisplacementFluidPatch[i][0] = data_to_recv[iterator++];
            pointDisplacementFluidPatch[i][1] = data_to_recv[iterator++];
            if (dim ==3)
                pointDisplacementFluidPatch[i][2] = data_to_recv[iterator++];

            std::cout << i << " : "<< pointDisplacementFluidPatch[i][0] << ", " << pointDisplacementFluidPatch[i][1] << ", " << pointDisplacementFluidPatch[i][2] << std::endl;
        } */

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
        std::cout<< "Data to be send from OpenFOAM: with array size = " << data_to_send.size() << std::endl;

        //Export this force array to CoSimulation
        //CoSimIO::Info connect_info;
        connect_info.Clear();
        connect_info.Set("identifier", "force_values");
        connect_info.Set("connection_name", connection_name);
        connect_info = CoSimIO::ExportData(connect_info, data_to_send);

        std::cout<< runTime_.timeName() << " : Data has been exported from OpenFOAM to CoSimulation: Force values" <<std::endl;
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
            std::cout << "Solver Type : Compressible " << std::endl;
        }
    }

    if(solverType_  == "unknown")
    {
        std::cout << "Solver Type: Neither Compressible nor Incompresible" << std::endl;
    }

    return solverType_;
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

//Calculate Total Force
bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::calculateForces(std::size_t interface_index)
{
    std::vector<int> patchIDs;
    // For every patch that participates in the coupling interface
    for (uint i = 0; i < interfaces_.at(interface_index).patchNames.size(); i++)
    {
        // Get the patchID
        int patchID = mesh_.boundaryMesh().findPatchID((interfaces_.at(interface_index).patchNames).at(i));

        // Throw an error if the patch was not found
        if (patchID == -1)
        {
            std::cout<< "ERROR: Patch " << (interfaces_.at(interface_index).patchNames).at(i) << " does not exist." << std::endl;
        }

    // Add the patch in the list
    patchIDs.push_back(patchID);
    }

    //- Force field
    Foam::volVectorField * Force_; //Access this real force values from OpenFOAM

    //Initialize the Force -> Need to check
    Force_ = new volVectorField
    (
        IOobject
        (
            "Force",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "fdim",
            dimensionSet(1,1,-2,0,0,0,0),
            Foam::vector::zero
        )
    );

    //Get different force fields from OpenFOAM, See Force Function Object
    //1. Normal vectors on the boundary, multiplied with the face areas
    const surfaceVectorField::Boundary& Sfb
    (
        mesh_.Sf().boundaryField()
    );

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

        //Pressure forces
        if(solverType_.compare("incompressible") == 0)
        {
            Force_->boundaryFieldRef()[patchID] = Sfb[patchID] * pb[patchID] * rhob[patchID];
        }
        else if(solverType_.compare("compressible") == 0)
        {
            Force_->boundaryFieldRef()[patchID] = Sfb[patchID] * pb[patchID];
        }
        else
        {
            FatalErrorInFunction << "Forces calculation does only support compressible or incompressible solver type." << exit(FatalError);
        }

        //Viscous forces
        Force_->boundaryFieldRef()[patchID] += Sfb[patchID] & devRhoReffb[patchID];

        // Now write this forces into Buffer to send further to the Strctural Solver
        int bufferIndex = 0;
        // For every cell of the patch
        forAll(Force_->boundaryFieldRef()[patchID], i)
        {
            // Copy the force into the buffer
            // x-dimension
            data_to_send[bufferIndex++] = Force_->boundaryFieldRef()[patchID][i].x();

            // y-dimension
            data_to_send[bufferIndex++] = Force_->boundaryFieldRef()[patchID][i].y();

            if(dim == 3)
            {
                // z-dimension
                data_to_send[bufferIndex++] = Force_->boundaryFieldRef()[patchID][i].z();
            }
        }

    }

/*     int i=0;
    for(auto& value : data_to_send)
    {
        std::cout << "id = " << i << ", Value = " << value << std::endl;
        i++;
    }
 */
    return true;
}

}//namespace Foam

// ************************************************************************* //
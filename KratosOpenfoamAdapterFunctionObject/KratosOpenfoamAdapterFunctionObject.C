/*-----------------------------------------------------------------------*\

Master-Thesis Work
Ashish Darekar

Sourcefile for the KratosOpenfoamAdpterFunctionObject.H

\*-----------------------------------------------------------------------*/

#include "KratosOpenfoamAdapterFunctionObject.H"

// OpenFOAM header files
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstreams.H" //Useful for IO operations

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(KratosOpenfoamAdapterFunctionObject, 0);
    addToRunTimeSelectionTable(functionObject, KratosOpenfoamAdapterFunctionObject, dictionary);
}
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

bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::read(const dictionary& dict)
{

    //Connect to CoSimIO
    //CoSimulationAdapter_.configure();
    //std::cout << "CoSimulation Adapter : Configuration" << std::endl;

    std::cout << "CoSimulation Adapter's function object : read()" << std::endl;

    //Reading configuration parameters from ControlDict related to Adapter
    fvMeshFunctionObject::read(dict);

    participantName_ = dict.lookupOrDefault<word>("participant", "fluid");

    std::cout<< "Name of the participant is: " << participantName_ <<std::endl;

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
            }
        }
    }

    std::cout << "Interfaces Reading: Done" << std::endl;

    //Importing data at the beginning from Kratos Cosimulation using CoSimIO
    CoSimIO::Info settings;
    settings.Set("my_name", "Openfoam_adapter");
    settings.Set("connect_to", "Kratos_CoSimulation");
    settings.Set("echo_level", 1);
    settings.Set("version", "1.25");

    connect_info = CoSimIO::Connect(settings);
    COSIMIO_CHECK_EQUAL(connect_info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Connected);
    connection_name = connect_info.Get<std::string>("connection_name");

    std::vector<double> receive_data;
    connect_info.Clear();
    connect_info.Set("identifier", "vector_of_pi");
    connect_info.Set("connection_name", connection_name);
    connect_info = CoSimIO::ImportData(connect_info, receive_data);

    for(auto& value : receive_data)
    {
        std::cout<< value << std::endl;
        COSIMIO_CHECK_EQUAL(value, 3.14);
    }

    return true;
}


bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::execute()
{
    //CoSimulationAdapter_.execute();

    return true;
}


bool Foam::functionObjects::KratosOpenfoamAdapterFunctionObject::end()
{
    //Dicsonect from CoSimIO
    //CoSimulationAdapter_.end();

    std::cout << "CoSimulation Adapter's function object : end()" << std::endl;

    //Exporting the dummy vector data from OpenFOAM to the CoSimulation using CoSimIO
    std::vector<double> data_to_send(4,3.14);
    connect_info.Clear();
    connect_info.Set("identifier", "vector_of_pi");
    connect_info.Set("connection_name", connection_name);
    connect_info = CoSimIO::ExportData(connect_info, data_to_send);

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

// ************************************************************************* //
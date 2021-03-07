/*-----------------------------------------------------------------------*\

Master-Thesis Work
Ashish Darekar

Sourcefile for the CoSimulationAdapter.H

\*-----------------------------------------------------------------------*/

#include "CoSimulationAdapter.H"

#define COSIMIO_CHECK_EQUAL(a, b)                            \
if (a != b) {                                                \
    std::cout << "in line " << __LINE__ << " : " << a        \
            << " is not equal to " << b << std::endl;        \
}

//using namespace Foam;

/* CoSimulationAdapter::CoSimulationAdapter(const Foam::Time& runTime, const Foam::fvMesh& mesh)
:
runTime_(runTime),
mesh_(mesh)
{
    std::cout << "CoSimulation Adapter is loaded" << std::endl;
    return;
}
 */

/* CoSimulationAdapter::CoSimulationAdapter()
{
    std::cout << "CoSimulation Adapter is loaded" << std::endl;
    return;
}

void CoSimulationAdapter::configure()
{
    std::cout << "CoSimulation Adapter : Configuration" << std::endl;

    CoSimIO::Info settings;
    settings.Set("my_name", "Openfoam_adapter");
    settings.Set("connect_to", "Kratos_CoSimulation");
    settings.Set("echo_level", 1);
    settings.Set("version", "1.25");

    connect_info = CoSimIO::Connect(settings);
    COSIMIO_CHECK_EQUAL(connect_info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Connected);
    connection_name = connect_info.Get<std::string>("connection_name");

    return;
}

void CoSimulationAdapter::execute()
{
    std::cout << "CoSimulation Adapter : Execute" << std::endl;

    std::vector<double> data_to_send(4,3.14);
    connect_info.Clear();
    connect_info.Set("identifier", "vector_of_pi");
    connect_info.Set("connection_name", connection_name);
    connect_info = CoSimIO::ExportData(connect_info, data_to_send);

    return;
}

void CoSimulationAdapter::end()
{
    std::cout << "CoSimulation Adapter : end" << std::endl;

    CoSimIO::Info disconnect_settings;
    disconnect_settings.Set("connection_name", connection_name);
    connect_info = CoSimIO::Disconnect(disconnect_settings);
    COSIMIO_CHECK_EQUAL(connect_info.Get<int>("conection_status"),CoSimIO::ConnectionStatus::Disconnected);

    return;
}

CoSimulationAdapter::~CoSimulationAdapter()
{
    std::cout << "CoSimulation Adapter : Destructor" << std::endl;
    return;
} */
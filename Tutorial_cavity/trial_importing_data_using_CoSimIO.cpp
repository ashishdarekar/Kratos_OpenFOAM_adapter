//Importing Data from Kratos_OpenFoam_adapter

// CoSimulation includes
#include "/home/ashish/Documents/MS/CoSimIO/co_sim_io/co_sim_io.hpp"

#define COSIMIO_CHECK_EQUAL(a, b)                                \
    if (a != b) {                                                \
        std::cout << "in line " << __LINE__ << " : " << a        \
                  << " is not equal to " << b << std::endl;      \
        return 1;                                                \
    }

int main()
{
    CoSimIO::Info settings;
    settings.Set("my_name", "Kratos_CoSimulation");
    settings.Set("connect_to", "Openfoam_adapter");
    settings.Set("echo_level", 1);
    settings.Set("version", "1.25");

    auto info = CoSimIO::Connect(settings);
    COSIMIO_CHECK_EQUAL(info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Connected);
    const std::string connection_name = info.Get<std::string>("connection_name");

    std::vector<double> receive_data;
    info.Clear();
    info.Set("identifier", "vector_of_pi");
    info.Set("connection_name", connection_name);
    info = CoSimIO::ImportData(info, receive_data);

    for(auto& value : receive_data)
    {
        std::cout<< value << std::endl;
        COSIMIO_CHECK_EQUAL(value, 3.14);
    }

    CoSimIO::Info disconnect_settings;
    disconnect_settings.Set("connection_name", connection_name);
    info = CoSimIO::Disconnect(disconnect_settings); // disconnect afterwards
    COSIMIO_CHECK_EQUAL(info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Disconnected);

    return 0;
}
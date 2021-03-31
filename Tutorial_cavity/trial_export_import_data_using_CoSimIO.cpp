//Exporting Data to Kratos_OpenFoam_adapter

// CoSimIO includes
#include "co_sim_io.hpp"

#define COSIMIO_CHECK_EQUAL(a, b)                                \
    if (a != b) {                                                \
        std::cout << "in line " << __LINE__ << " : " << a        \
                  << " is not equal to " << b << std::endl;      \
        return 1;                                                \
    }

int main()
{
    // ************ Connection setting ******************//
    CoSimIO::Info settings;
    settings.Set("my_name", "Kratos_CoSimulation");
    settings.Set("connect_to", "Openfoam_adapter");
    settings.Set("echo_level", 1);
    settings.Set("version", "1.25");

    auto info = CoSimIO::Connect(settings);
    COSIMIO_CHECK_EQUAL(info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Connected);
    const std::string connection_name = info.Get<std::string>("connection_name");

    // *********** Import Mesh **************//
    std::cout << "Importing All interface meshes: Start" << std::endl;
    std::vector<std::unique_ptr<CoSimIO::ModelPart>> model_part_interfaces_;
    int num_interfaces_ = 2 ; //We need to transfer this data as well to CoSimulation

    for(std::size_t j=0; j < num_interfaces_; j++)
    {
        std::string interface_name = "interface" + std::to_string(j+1);
        model_part_interfaces_.push_back(CoSimIO::make_unique<CoSimIO::ModelPart>(interface_name));

        info.Clear();
        info.Set("identifier", "fluid_mesh");
        info.Set("connection_name", connection_name); // connection_name is obtained from calling "Connect"

        auto import_info = CoSimIO::ImportMesh(info, *model_part_interfaces_.at(j));

        std::cout << "Name of the imported mesh is " << model_part_interfaces_.at(j)->Name() << std::endl;
    }
    std::cout << "Importing All interface meshes: Done" << std::endl;

    // ************ Importing Data ******************//
    std::cout << "CoSimIO is tyring to Import the data" << std::endl;
    std::vector<double> receive_data;
    info.Clear();
    info.Set("identifier", "pressure_values");
    info.Set("connection_name", connection_name);
    info = CoSimIO::ImportData(info, receive_data);

    std::cout<< "Data Received from OpneFOAM:" << std::endl;

    for(auto& value : receive_data){
        //std::cout<< value << std::endl;
        //COSIMIO_CHECK_EQUAL(value, 3.14);
    }

    // ************ DisConnection setting ******************//
    CoSimIO::Info disconnect_settings;
    disconnect_settings.Set("connection_name", connection_name);
    info = CoSimIO::Disconnect(disconnect_settings); // disconnect afterwards
    COSIMIO_CHECK_EQUAL(info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Disconnected);

    return 0;
}


//-Dummy file To Commnicate between Openfoam_Adapter and Openfoam_Kratos_Wrapper

//-CoSimIO includes
#include "co_sim_io.hpp"

#define COSIMIO_CHECK_EQUAL(a, b)                                \
    if (a != b) {                                                \
        std::cout << "in line " << __LINE__ << " : " << a        \
                  << " is not equal to " << b << std::endl;      \
        return 1;                                                \
    }

int makedummydisplacemetvalues(std::vector<double>& send_data, int size)
{
    int i = 0;
    while(i < size)
    {
        send_data[i++] = i * 1.01;
        send_data[i++] = i * 1.08;
        send_data[i++] = i * 1.05;
    }
    return 0;
}

int main()
{
    // ********************************* Connection setting **********************************************//
    CoSimIO::Info settings;
    settings.Set("my_name", "Openfoam_Kratos_Wrapper");
    settings.Set("connect_to", "Openfoam_Adapter");
    settings.Set("echo_level", 1);
    settings.Set("version", "1.25");

    auto info = CoSimIO::Connect(settings);
    COSIMIO_CHECK_EQUAL(info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Connected);
    const std::string connection_name = info.Get<std::string>("connection_name");

    // ***************************************** Import Mesh *********************************************//
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

    // ***************************************** Export and Import Data *********************************************//
    int num_of_nodes = 182;
    int dim = 3;

    double time_step = 0.001;
    while(time_step < 0.011)
    {
        //Exporting Displacement values from Structural Solver to OpenFOAM
        std::cout << "CoSimIO is trying to Export the data: Displacement Values" << std::endl;
        std::vector<double> send_data;

        // Resize it according to the number of the nodes in the Cooupling Interface "Flap" and then Fill the dummy values in it
        send_data.resize(num_of_nodes * dim);
        makedummydisplacemetvalues(send_data, num_of_nodes * dim);

        info.Clear();
        info.Set("identifier", "disp_values");
        info.Set("connection_name", connection_name);
        info = CoSimIO::ExportData(info, send_data);
        std::cout << time_step << " : Data sent to OpenFOAM: displacement Values with array size = " << send_data.size() << "\n" << std::endl;

        //In between we should have one iteration of Coupling Solver/ CoSimulaton in our case, Need to write the code for it

        //Importing Force Values from OpenFOAM
        std::cout << "CoSimIO is trying to Import the data: Force Values" << std::endl;
        std::vector<double> receive_data;
        info.Clear();
        info.Set("identifier", "force_values");
        info.Set("connection_name", connection_name);
        info = CoSimIO::ImportData(info, receive_data);
        std::cout << time_step << " : Data Received from OpenFOAM: Force Values with array size = " << receive_data.size() << "\n" << "\n" << std::endl;

        time_step+=0.001;

        /* int i = 0;
        for(auto& value : receive_data)
        {
            std::cout << "id = " << i << ", Value = " << value << std::endl;
            i++;
        } */
    }

    // ************************************** DisConnection setting ******************************************//
    CoSimIO::Info disconnect_settings;
    disconnect_settings.Set("connection_name", connection_name);
    info = CoSimIO::Disconnect(disconnect_settings); // disconnect afterwards
    COSIMIO_CHECK_EQUAL(info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Disconnected);

    return 0;
}




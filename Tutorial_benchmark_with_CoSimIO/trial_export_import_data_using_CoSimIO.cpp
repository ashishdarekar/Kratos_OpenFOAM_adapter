//-Dummy file To Commnicate between Openfoam_Adapter and Openfoam_Kratos_Wrapper

//-CoSimIO includes
#include "co_sim_io.hpp"
#include <math.h>

#define COSIMIO_CHECK_EQUAL(a, b)                                \
    if (a != b) {                                                \
        std::cout << "in line " << __LINE__ << " : " << a        \
                  << " is not equal to " << b << std::endl;      \
        return 1;                                                \
    }


//Hard Coded for Number of nodes = 182, Just to test the working of Pimple Solver
int makedummydisplacemetvalues(std::vector<double>& send_data, double time_step, int size)
{
    int num_nodes = size / 3;
    int i=0;
    int itr =0;
    double a = 1.0;

    while(i < num_nodes)
    {
        if( i < 20 )//i > 161
        {
            itr = i*3;
            send_data[itr] = 0; //X
            send_data[itr+1] = (pow( a, time_step))/10; //Y
            send_data[itr+2] = 0; //Z

            itr = (i+162)*3;
            send_data[itr] = 0; //X
            send_data[itr+1] = (pow( a, time_step))/10; //Y
            send_data[itr+2] = 0; //Z

            a+=0.04;
        }
        else if( i < 88) //i>93 covered
        {
            itr = i*3;
            send_data[itr] = 0; //X
            send_data[itr+1] = (pow( a, time_step))/10; //Y
            send_data[itr+2] = 0; //Z

            itr = (i+74)*3;
            send_data[itr] = 0; //X
            send_data[itr+1] = (pow( a, time_step))/10; //Y
            send_data[itr+2] = 0; //Z

            a+=0.04;
        }
        else if( i == 88 || i == 89)
        {
            itr = i*3;
            send_data[itr] = 0; //X
            send_data[itr+1] = (pow( a, time_step))/10; //Y
            send_data[itr+2] = 0; //Z

            itr = (i+4)*3;
            send_data[itr] = 0; //X
            send_data[itr+1] = (pow( a, time_step))/10; //Y
            send_data[itr+2] = 0; //Z

            a+=0.04;
        }
        else if(i==90 || i == 91)
        {
            itr = i*3;
            send_data[itr] = 0; //X
            send_data[itr+1] = (pow( a, time_step))/10; //Y
            send_data[itr+2] = 0; //Z
            a+=0.04;
        }
        i++;

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
/*     std::cout << "Importing All interface meshes: Start" << std::endl;
    std::vector<std::unique_ptr<CoSimIO::ModelPart>> model_part_interfaces_;
    int num_interfaces_ = 2 ; //We need to transfer this data as well to CoSimulation

    for(std::size_t j=0; j < num_interfaces_; j++)
    {
        std::string interface_name = "interface" + std::to_string(j+1);
        model_part_interfaces_.push_back(CoSimIO::make_unique<CoSimIO::ModelPart>(interface_name));

        info.Clear();
        info.Set("identifier", "interface_flap");
        info.Set("connection_name", connection_name); // connection_name is obtained from calling "Connect"

        auto import_info = CoSimIO::ImportMesh(info, *model_part_interfaces_.at(j));

        std::cout << "Name of the imported mesh is " << model_part_interfaces_.at(j)->Name() << std::endl;
    }
    std::cout << "Importing All interface meshes: Done" << std::endl; */


    // ***************************************** Iport Mesh for 1 interface case *************************************//

    std::cout << "Importing one interface meshe: Start" << std::endl;
    int num_interfaces_ = 1 ; //We need to transfer this data as well to CoSimulation

    // create CoSimIO::ModelPart
    CoSimIO::ModelPart model_part_interface_flap("interface_flap_model_part");

    for(std::size_t j=0; j < num_interfaces_; j++)
    {
        info.Clear();
        info.Set("identifier", "interface_flap");
        info.Set("connection_name", connection_name); // connection_name is obtained from calling "Connect"

        auto import_info = CoSimIO::ImportMesh(info, model_part_interface_flap);
    }
    std::cout << "Nodes of the part imported: " << model_part_interface_flap.NumberOfNodes() <<std::endl;
    std::cout << "Elements of the part imported: " << model_part_interface_flap.NumberOfElements() <<std::endl;

    std::cout << "Importing All interface meshes: Done" << std::endl;


    // ***************************************** Export and Import Data *********************************************//
    int num_of_nodes = model_part_interface_flap.NumberOfNodes();
    int dim = 3;

    double time_step = 0.0;
    while(time_step < 25)
    {
        //Exporting Displacement values from Structural Solver to OpenFOAM
        std::cout << "CoSimIO is trying to Export the data: Displacement Values" << std::endl;
        std::vector<double> send_data;

        // Resize it according to the number of the nodes in the Cooupling Interface "Flap" and then Fill the dummy values in it
        send_data.resize(num_of_nodes * dim);
        makedummydisplacemetvalues(send_data, time_step, num_of_nodes * dim);

        info.Clear();
        info.Set("identifier", "disp");
        info.Set("connection_name", connection_name);
        info = CoSimIO::ExportData(info, send_data);
        std::cout << time_step << " : Data sent to OpenFOAM: displacement Values with array size = " << send_data.size() << "\n" << std::endl;

        //In between we should have one iteration of Coupling Solver/ CoSimulaton in our case, Need to write the code for it

        //Importing Force Values from OpenFOAM
        std::cout << "CoSimIO is trying to Import the data: Force Values" << std::endl;
        std::vector<double> receive_data;
        info.Clear();
        info.Set("identifier", "load");
        info.Set("connection_name", connection_name);
        info = CoSimIO::ImportData(info, receive_data);
        std::cout << time_step << " : Data Received from OpenFOAM: Force Values with array size = " << receive_data.size() << "\n" << "\n" << std::endl;

        time_step+=0.005;

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




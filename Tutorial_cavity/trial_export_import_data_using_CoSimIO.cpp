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

    // ****************************** Code to test the working of CosimIO::ImportMesh ***********************************************//
    //Checking the nodal coordinates of the received meshes
    for(int i=1; i<5; i++)
    {
        std::cout << "Coordinates of node with Id "<< i << " in interface1 are: (" << model_part_interfaces_.at(0)->GetNode(i).X() << "," << model_part_interfaces_.at(0)->GetNode(i).Y()
        << "," << model_part_interfaces_.at(0)->GetNode(i).Z() << ")." << std::endl;
    }

    //Checking the Elements of the received meshes
    for(int i=1; i<20; i+=2)
    {
        CoSimIO::Element& element = model_part_interfaces_.at(0)->GetElement(i);
        std::cout << "Element info in the interface1 is: ( Id= " << element.Id() << ", Number of Nodes = " << element.NumberOfNodes() << ")" << std::endl;

        // access Id of element:
        CoSimIO::IdType element_id = element.Id();
        // the type can be accessed:
        CoSimIO::ElementType element_type = element.Type(); // e.g. CoSimIO::ElementType::Point2D or CoSimIO::ElementType::Line2D2
        std::cout << "Element Type:( " << typeid(element_type).name() << " )" << std::endl;

        // iterate the nodes of the element:
        for (auto node_it=element.NodesBegin(); node_it!=element.NodesEnd(); ++node_it)
        {
            CoSimIO::Node& node = **node_it;
            std::cout << "Id:( " << node.X() << "," << node.Y() << "," << node.Z() << " )" << std::endl;
        }
    }
    // ****************************** Code to test the working of CosimIO::ImportMesh ***********************************************//

    // ************ Importing Data ******************//
    std::cout << "CoSimIO is trying to Import the data: Pressure Values" << std::endl;
    std::vector<double> receive_data;
    info.Clear();
    info.Set("identifier", "pressure_values");
    info.Set("connection_name", connection_name);
    info = CoSimIO::ImportData(info, receive_data);

    std::cout<< "Data Received from OpenFOAM: Pressure values with array size =" << receive_data.size() << std::endl;

    int i = 0;

    for(auto& value : receive_data){
        std::cout << "id = " << i << ", Value = " << value << std::endl;
        i++;
        //COSIMIO_CHECK_EQUAL(value, 3.14);
    }

    // ************ DisConnection setting ******************//
    CoSimIO::Info disconnect_settings;
    disconnect_settings.Set("connection_name", connection_name);
    info = CoSimIO::Disconnect(disconnect_settings); // disconnect afterwards
    COSIMIO_CHECK_EQUAL(info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Disconnected);

    return 0;
}


// OpenFOAM header files
//#include "Time.H"
//#include "fvMesh.H"

//CoSimIO header file
#include "co_sim_io.hpp"

//namespace Kratos
//{

#define COSIMIO_CHECK_EQUAL(a, b)                            \
if (a != b) {                                                \
    std::cout << "in line " << __LINE__ << " : " << a        \
            << " is not equal to " << b << std::endl;        \
}


class CoSimulationAdapter
{

private:
        //- Variables required  for CoSimIO

        CoSimIO::Info connect_info;
        std::string connection_name;

/*     //- OpenFOAM runTime object
    const Foam::Time& runTime_;

    //- OpenFOAM fvMesh object
    const Foam::fvMesh& mesh_;

    // Vectors of pointers to the fields(These are Registry of Objects)

        //- volScalarField fields
        std::vector<Foam::volScalarField*> volScalarFields_;

        //- volVectorField fields
        std::vector<Foam::volVectorField*> volVectorFields_;

        //- surfaceScalarField fields
        std::vector<Foam::surfaceScalarField*> surfaceScalarFields_;

        //- surfaceVectorField fields
        std::vector<Foam::surfaceVectorField*> surfaceVectorFields_;

        //- pointScalarField fields
        std::vector<Foam::pointScalarField*> pointScalarFields_;

        //- pointVectorField fields
        std::vector<Foam::pointVectorField*> pointVectorFields_;

    // Methods communicating with CoSimulation

        //- Initialize CoSimulation and exchange the first data
        void initialize();

        //- Finalize and destroy CoSimulation
        void finalize();

        //- Advance CoSimulaton
        void advance();
*/


public:

    // Methods called by the functionObject

        //- Constructor
        //CoSimulationAdapter(const Foam::Time& runTime, const Foam::fvMesh& mesh);

        //Default Constructor
        CoSimulationAdapter(){
            std::cout << "CoSimulation Adapter : Constructor" << std::endl;
            return;
        }

        //- Setup the adapter's configuration
        void configure(){
            std::cout << "CoSimulation Adapter : Configuration" << std::endl;

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

            CoSimIO::Info disconnect_settings;
            disconnect_settings.Set("connection_name", connection_name);
            connect_info = CoSimIO::Disconnect(disconnect_settings); // disconnect afterwards
            COSIMIO_CHECK_EQUAL(connect_info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Disconnected);

            return;
        }

        //- Called by the functionObject's execute()
        void execute()
        {
            std::cout << "CoSimulation Adapter : Execute" << std::endl;

            //Send the dummy vector data at the end of each iteration

            /* std::vector<double> data_to_send(4,3.14);
            connect_info.Clear();
            connect_info.Set("identifier", "vector_of_pi");
            connect_info.Set("connection_name", connection_name);
            connect_info = CoSimIO::ExportData(connect_info, data_to_send); */

            return;
        }

        //- Called by the functionObject's end()
        void end()
        {
            std::cout << "CoSimulation Adapter : end" << std::endl;

            CoSimIO::Info settings;
            settings.Set("my_name", "Openfoam_adapter");
            settings.Set("connect_to", "Kratos_CoSimulation");
            settings.Set("echo_level", 1);
            settings.Set("version", "1.25");

            connect_info = CoSimIO::Connect(settings);
            COSIMIO_CHECK_EQUAL(connect_info.Get<int>("connection_status"), CoSimIO::ConnectionStatus::Connected);
            connection_name = connect_info.Get<std::string>("connection_name");

            //Send the dummy vector data at the end of the simulation
            std::vector<double> data_to_send(4,3.14);
            connect_info.Clear();
            connect_info.Set("identifier", "vector_of_pi");
            connect_info.Set("connection_name", connection_name);
            connect_info = CoSimIO::ExportData(connect_info, data_to_send);

            CoSimIO::Info disconnect_settings;
            disconnect_settings.Set("connection_name", connection_name);
            connect_info = CoSimIO::Disconnect(disconnect_settings);
            COSIMIO_CHECK_EQUAL(connect_info.Get<int>("conection_status"),CoSimIO::ConnectionStatus::Disconnected);

            return;
        }

        //- Destructor
        ~CoSimulationAdapter()
        {
            std::cout << "CoSimulation Adapter : Destructor" << std::endl;
            return;
        }

};

//}


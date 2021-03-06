// OpenFOAM header files
//#include "Time.H"
//#include "fvMesh.H"

//CoSimIO header file
#include "/home/ashish/Documents/MS/CoSimIO/co_sim_io/co_sim_io.hpp"

//namespace Kratos
//{

class CoSimulationAdapter
{

private:

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
        CoSimulationAdapter();

        //- Setup the adapter's configuration
        void configure();

        //- Called by the functionObject's execute()
        void execute();

        //- Called by the functionObject's end()
        void end();

        //- Destructor
        ~CoSimulationAdapter();

};

//}


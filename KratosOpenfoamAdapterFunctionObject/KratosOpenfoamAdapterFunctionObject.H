/*-----------------------------------------------------------------------*\

Master-Thesis Work
Ashish Darekar

Class
    Foam::functionObjects::KratosOpenfoamAdapterFunctionObject

Description
    This function object tries to access the real time field data
    from OpenFOAM solver and make it available to the COsimIO
    for further coupling in Kratos-CoSimualtion application

    NOTE:This function object is created using the template given in the $FOAM_ETC/codeTemplates/ directory.

    Example of function object specification:
    \verbatim
    CoSimulation_Adapter
    {
        type KratosOpenfoamAdapterFunctionObject;
        libs ("libKratosOpenfoamAdapterFunctionObject.so");
    }
    \endverbatim

Usage
    \table
        Property            | Description                                                                       | Required | Default value
        type                | type name: KratosOpenfoamAdapterFunctionObject                                    | yes      |
        my_name             | Name of the participant                                                           | no       | Fluid
        interfaces          | To save the coupling interfaces(Each interface is a subdirectory)                 | yes      |
            mesh            | Name of the coupling mesh                                                         | yes      |
            patches         | WordList of the openFOAM patches which contain the coupling faces                 | yes      |
            locations       | Name of the location in mesh to get the required data (faceCenter/facenode)       | yes      | faceCenter
            importData        | WordList of the variables to Read by OpenFOAM                                     | yes      |
            exportData       | WordList of the variables to write by OpenFOAM                                    | yes      |
        parameters          | List of Parameters required for FSI: rho, mu                                      | no       |
    \endtable

SourceFiles
    KratosOpenfoamAdapterFunctionObject.C

\*-----------------------------------------------------------------------*/

#ifndef KRATOSOPENFOAMADAPTERFUNCTIONOBJECT_H
#define KRATOSOPENFOAMADAPTERFUNCTIONOBJECT_H

// OpenFOAM header files
#include "fvMeshFunctionObject.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H" //To use volumeVectorfield

#include "addToRunTimeSelectionTable.H"
#include "IOstreams.H" //Useful for IO operations

//Displacement Calculation
#include "pointMesh.H" //To use pointVectorfield
#include "pointPatchField.H" //To use pointVectorfield
#include "fvCFD.H"
#include "fixedValuePointPatchFields.H"

//Force Calculation
#include "pointFields.H"
#include "vectorField.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

//For MPI
/* #include "UPstream.H"
#include "PstreamReduceOps.H"
#include "OSspecific.H"
#include "PstreamGlobals.H"
#include "SubList.H"
#include "allReduce.H"
 */

// Main adapter header files
// #include "CoSimulationAdapter.H"

//CoSimIO header file
#include "co_sim_io_mpi.hpp"

#define COSIMIO_CHECK_EQUAL(a, b)                                 \
if (a != b) {                                                     \
std::cout << "in line " << __LINE__ << " : " << a                 \
    << " is not equal to " << b << std::endl;                     \
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class KRATOSOPENFOAMADAPTERFUNCTIONOBJECT Declaration
\*---------------------------------------------------------------------------*/

class KratosOpenfoamAdapterFunctionObject
:
    public fvMeshFunctionObject
{

private:

    const Foam::Time& runTime_; //Copy of RunTime Time object
    const dictionary& dict_; //Copy of Control Dict

    //- Member variables required for CoSimIO
    std::string connection_name;

    // CoSimIO info object- one is sufficient
    CoSimIO::Info info;

    //- Member for OpenFOAM
    std::size_t time_step_ = 0 ;
    std::size_t num_interfaces_ = 0;

    //- Solver type compressible/incompressible
    std::string solverType_ = "unknown";

    int dim = 3; //Dimension of the problem, need to check how we can get this? We will provide this in config file.
    double thick = 1.0; //Thickness of the problem , Added to make the force value independent of the Thickess for 2D

    // Private Member Functions

        //- Disallow default bitwise copy construct
        KratosOpenfoamAdapterFunctionObject(const KratosOpenfoamAdapterFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const KratosOpenfoamAdapterFunctionObject&) = delete;

        //- Private function for CoSimIO
        std::vector<double> data_to_send;
        std::vector<double> data_to_recv;

        //- For making geometry of a interface Mesh(Model Part)
        std::vector<vector> array_of_nodes;

        // Parallel Related
        int TotalNumOfProcesse = Pstream::nProcs();
        int MyRank = Pstream::myProcNo();

protected:
    //- Protected Data to set the parameters required for coupling

        //- Parameters for Adapter's configuration
        Foam::word my_name;

        //- For Force Field
        //- Initialize the Force -> Need to check
        Foam::volVectorField*  Force_ = new volVectorField
        (
            IOobject
            (
                "Force",
                mesh_.time().timeName(),
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

        //- Create different interface for different FSI coupling interface (or OpenFOAM patches)
        //- This strucutre is made to save all the data related to an interface
        struct InterfaceData
        {
            std::string nameOfInterface;
            std::vector<std::string> patchNames;
            std::vector<std::string> importData;
            std::vector<std::string> importDataIdentifier;
            std::vector<std::string> exportData;
            std::vector<std::string> exportDataIdentifier;
            std::vector<double> data_to_send;
            std::vector<double> data_to_recv;
            uint numNodes = 0;
            uint numElements = 0;
        };

        //- Array of strctures to store data for all coupling interfaces
        std::vector<struct InterfaceData> interfaces_;

        //- Array of CoSimIO::Model parts for all coupling interfaces
        std::vector<std::unique_ptr<CoSimIO::ModelPart>> model_part_interfaces_;

        //- Auxiliar functions
        bool calculateForces(std::size_t interface_index);
        std::string determineSolverType();

        bool is_same_points(Foam::vector& pointX, Foam::vector& pointY);
        int compare_nodes(Foam::vector& pointX);

        //- Stress tensor (see the OpenFOAM "Forces" function object)
        Foam::tmp<Foam::volSymmTensorField> devRhoReff() const;

        Foam::tmp<Foam::volScalarField> rho() const;

        Foam::tmp<Foam::volScalarField> mu() const;

        Foam::vectorField getFaceVectors(const unsigned int patchID) const;

public:

    //- Runtime type information
    TypeName("KratosOpenfoamAdapterFunctionObject");

    //- Constructors

        //- Construct from Time and dictionary
        KratosOpenfoamAdapterFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

    //- Destructor
    virtual ~KratosOpenfoamAdapterFunctionObject();

    // Member Functions

        //- Read the data/after each change and connection with CoSimIO
        virtual bool read(const dictionary&);

        //- Execute, Exchaging required data on the coupling interfaces
        virtual bool execute();

        //- Execute at the final time-loop/ Disconection from CoSimIO
        virtual bool end();

        //- Write data after each loop??/ Currently not doing anything
        virtual bool write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

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

// Main adapter header files
// #include "CoSimulationAdapter.H"

//CoSimIO header file
#include "co_sim_io.hpp"

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
        //std::vector<std::unique_ptr<vector>> array_of_nodes;


        //- Implicit Coupling related vector fields and thier copies
        //- Implicit related volScalarField fields
        std::vector<Foam::volScalarField*> volScalarFields_;

        //- Implicit related volScalarField fields (copies)
        std::vector<Foam::volScalarField*> volScalarFieldCopies_;

        //- Implicit related volVectorField fields
        std::vector<Foam::volVectorField*> volVectorFields_;

        //- Implicit related volVectorField fields (copies)
        std::vector<Foam::volVectorField*> volVectorFieldCopies_;

        //- Implicit related surfaceScalarField fields
        std::vector<Foam::surfaceScalarField*> surfaceScalarFields_;

        //- Implicit related surfaceScalarField fields (copies)
        std::vector<Foam::surfaceScalarField*> surfaceScalarFieldCopies_;

        //- Implicit related surfaceVectorField fields
        std::vector<Foam::surfaceVectorField*> surfaceVectorFields_;

        //- Implicit related surfaceVectorField fields (copies)
        std::vector<Foam::surfaceVectorField*> surfaceVectorFieldCopies_;

        //- Implicit related pointScalarField fields
        std::vector<Foam::pointScalarField*> pointScalarFields_;

        //- Implicit related pointScalarField fields (copies)
        std::vector<Foam::pointScalarField*> pointScalarFieldCopies_;

        //- Implicit related pointVectorField fields
        std::vector<Foam::pointVectorField*> pointVectorFields_;

        //- Implicit related pointVectorField fields (copies)
        std::vector<Foam::pointVectorField*> pointVectorFieldCopies_;

        //- Implicit related volTensorField fields
        std::vector<Foam::volTensorField*> volTensorFields_;

        //- Implicit related volTensorField fields (copies)
        std::vector<Foam::volTensorField*> volTensorFieldCopies_;

        //- Implicit related surfaceTensorField fields
        std::vector<Foam::surfaceTensorField*> surfaceTensorFields_;

        //- Implicit related surfaceTensorField fields (copies)
        std::vector<Foam::surfaceTensorField*> surfaceTensorFieldCopies_;

        //- Implicit related pointTensorField fields
        std::vector<Foam::pointTensorField*> pointTensorFields_;

        //- Implicit related pointTensorField fields (copies)
        std::vector<Foam::pointTensorField*> pointTensorFieldCopies_;

        //- Implicit related volSymmTensorField fields
        std::vector<Foam::volSymmTensorField*> volSymmTensorFields_;

        //- Implicit related volSymmTensorField fields (copies)
        std::vector<Foam::volSymmTensorField*> volSymmTensorFieldCopies_;

        //- Checkpointed time (value)
        Foam::scalar couplingIterationTimeValue_;

        //- Checkpointed time (index)
        Foam::label couplingIterationTimeIndex_;

        //- Checkpointed mesh points
        Foam::pointField meshPoints_;
        Foam::pointField oldMeshPoints_;

        bool isSimulationConverged = 0;
        int nMaxImplicit;
        int nIterationImplicit = 0;


protected:
    //- Protected Data to set the parameters required for coupling

        //- Parameters for Adapter's configuration
        Foam::word my_name;
        Foam::word coupling_scheme;

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
            std::vector<std::string> patchNames;
            std::vector<std::string> importData;
            std::vector<std::string> exportData;
            uint numNodes = 0;
            uint numElements = 0;
        };

        //- Array of strctures to save all interfaces
        std::vector<struct InterfaceData> interfaces_;

        //- Array of CoSimIO::Model parts for all coupling interfaces
        //std::vector<CoSimIO::ModelPart> model_part_interfaces_;
        std::vector<std::unique_ptr<CoSimIO::ModelPart>> model_part_interfaces_;

        //- create CoSimIO::ModelPart -> Not working
        //CoSimIO::ModelPart model_part_interface_flap("interface_flap_model_part");

        //- Auxiliar functions
        bool calculateForces(std::size_t interface_index);
        std::string determineSolverType();

        int compare_nodes(Foam::vector& pointX);

        //- Stress tensor (see the OpenFOAM "Forces" function object)
        Foam::tmp<Foam::volSymmTensorField> devRhoReff() const;

        Foam::tmp<Foam::volScalarField> rho() const;

        Foam::tmp<Foam::volScalarField> mu() const;

        Foam::vectorField getFaceVectors(const unsigned int patchID) const;

        //- For Implicit Coupling
        void setupImplicitCoupling();

        // Make a copies of Registered Object fields, depending on the type

        //- Copy of volScalarField
        void makeCopyOfObject(volScalarField& field);

        //- Copy of volVectorField
        void makeCopyOfObject(volVectorField& field);

        //- Copy of surfaceScalarField
        void makeCopyOfObject(surfaceScalarField& field);

        //- Copy of surfaceVectorField
        void makeCopyOfObject(surfaceVectorField& field);

        //- Copy of pointScalarField
        void makeCopyOfObject(pointScalarField& field);

        //- Copy of pointVectorField
        void makeCopyOfObject(pointVectorField& field);

        //- Copy of volTensorField
        void makeCopyOfObject(volTensorField& field);

        //- Copy of surfaceTensorField
        void makeCopyOfObject(surfaceTensorField& field);

        //- Copy of pointTensorField
        void makeCopyOfObject(pointTensorField& field);

        //- Copy of volSymmTensorField
        void makeCopyOfObject(volSymmTensorField& field);

        //- Writing all the fields for that time step required in an Implicit coupling
        void writeDataFieldsForImplicitCoupling();

        //- Reading all the fields for that time step required in an Implicit coupling
        void readDataFieldsForImplicitCoupling();

        void storeCheckpointTime();

        void reloadCheckpointTime();

        void storeMeshPoints();

        void reloadMeshPoints();

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

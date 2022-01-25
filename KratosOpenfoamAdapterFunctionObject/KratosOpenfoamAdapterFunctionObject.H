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
        libs ("<Path_to_KratosOpenfoamAdapterFunctionObject>/libKratosOpenfoamAdapterFunctionObject.so");
    }
    \endverbatim

Usage
    \table
        Property                    | Description                                                                       | Required | Default value
        type                        | type name: KratosOpenfoamAdapterFunctionObject                                    | yes      |
        libs                        | Excutable: ("libKratosOpenfoamAdapterFunctionObjectFunctionObjects.so")           | yes      |
        dim                         | Dimension of a problem                                                            | yes      |      3
        interfaces                  | To save the coupling interfaces(Each interface is a subdirectory)                 | yes      |
            name                    | Name of the coupling interface                                                    | yes      |
            patches                 | WordList of the OpenFOAM patches which contain the coupling faces                 | yes      |
            importDataIdentifier    | Import DataItendifiers required for communication with Kratos                     | yes      |
            exportDataIdentifier    | Export DataItendifiers required for communication with Kratos                     | yes      |
        parameters                  | List of Parameters required for FSI: rho, nu                                      | no       |
        debugLevel                  | Debug level either "info" or "debug" to see printed information                   | no       |     info
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

    //- Copy of OpenFOAM's RunTime Time object
    const Foam::Time& runTime_;

    //- Copy of OpenFOAM's Control Dict
    const dictionary& dict_; //Copy of Control Dict

    //- Member variables required for CoSimIO
    std::string connection_name;

    //- CoSimIO info object
    CoSimIO::Info info;

    //- For comparing with the nodes, while making mesh of a interface(Model Part)
    std::vector<vector> array_of_nodes;

    //- Total number of processes in MPI Parallel Simulations
    int TotalNumOfProcesses = Pstream::nProcs();

    //- Rank of the process
    int MyRank = Pstream::myProcNo();

    //- Node class : for an Adapter to save Node data of an interface mesh
    class Node
    {
        private:
            int localNodeIndex_;
            int nodeIndexforCoSim_;
            vector position_;
            std::vector<bool>commonWithWhichRanks; //With how many Ranks this node is shared
            std::vector<double>force_value; //To Store Nodal force Values

        public:
            Node(vector position, int localNodeIndex, int nodeIndexforCoSim, int commonRankFirstTime) : commonWithWhichRanks(Pstream::nProcs(),0){
                localNodeIndex_ = localNodeIndex;
                nodeIndexforCoSim_ = nodeIndexforCoSim;
                position_ = position;
                setCommonWithRank(commonRankFirstTime);
                force_value.resize(3);
                force_value[0] = 0.0;
                force_value[1] = 0.0;
                force_value[2] = 0.0;
            };

            vector getNodePosition(){
                return position_;
            };

            uint getNodeIndexForCoSim(){
                return nodeIndexforCoSim_;
            };

            void setNodeIndexForCoSim(uint index){
                nodeIndexforCoSim_ = index;
            };

            uint getLocalNodeIndex(){
                return localNodeIndex_;
            };

            void setLocalNodeIndex(uint index){
                localNodeIndex_ = index;
            };

            void setCommonWithRank(int rank)
            {
                if(rank != Pstream::myProcNo()) //As if rank ==MyRank, it means  node is local
                {
                    commonWithWhichRanks.at(rank) = 1 ;
                }
            };

            std::vector<bool> getCommonWithRank(){
                return commonWithWhichRanks;
            };

            int getHighestCommRank()
            {
                int highest_Rank = Pstream::nProcs(); //for local number returning highest number
                for(std::size_t i = 0 ;  i < commonWithWhichRanks.size() ; i++)
                {
                    if(commonWithWhichRanks[i] == 1)
                    {
                        highest_Rank = i;
                    }
                }

                if(highest_Rank > Pstream::myProcNo())
                {
                    return highest_Rank;
                }
                else
                {
                    return Pstream::nProcs(); //Returning highest number
                }
            };

            void setLoadValues(std::vector<double>& value){
                force_value[0] = value[0];
                force_value[1] = value[1];
                force_value[2] = value[2];
            };

            std::vector<double>& getLoadValues(){
                return force_value;
            };

            void setLoadValuesToZeros(){
                force_value[0] = 0.0;
                force_value[1] = 0.0;
                force_value[2] = 0.0;
            };
    };

    //- Element class : for an Adapter to save Element data of an interface mesh
    class Element
    {
        private:
            uint localElementIndex_;
            std::vector<uint> elementalNodeIndexes_;
            std::vector<Node> elementalNodes_; //Not useful unless it contains the actual nodes make it Pointer
            double forceValue;

        public:
            Element(uint localElementIndex){
                localElementIndex_ = localElementIndex;
            };

            uint getLocalElementIndex(){
                return localElementIndex_;
            };

            void addNodeIndexInList(uint nodeIndex){
                elementalNodeIndexes_.push_back(nodeIndex);
            };

            std::vector<uint>& getElementalNodeIndexes(){
                return elementalNodeIndexes_;
            };

            void addNodesInList(Node nodei){
                elementalNodes_.push_back(nodei);
            };

            std::vector<Node> getElementalNodes(){
                return elementalNodes_;
            };

            int getNumberOfNodes(){
                return elementalNodeIndexes_.size();
            };
    };

    //- Disallow default bitwise copy construct
    KratosOpenfoamAdapterFunctionObject(const KratosOpenfoamAdapterFunctionObject&);

    //- Disallow default bitwise assignment
    void operator=(const KratosOpenfoamAdapterFunctionObject&) = delete;

protected:
    //- Dimension of the problem(default = 3)
    int dim ;

    //- Thickness of the problem(default = 1.0)
    double thick ;

    //- Debug Level for printing the Adapter related information, 2 Levels availale: "debug" and "info"(default)
    std::string debugLevel ;

    //- Solver type compressible/incompressible
    std::string solverType_;

    //- Number of Coupling interfaces
    std::size_t num_interfaces_ = 0;

    //- Create different interface for different FSI coupling interface (or OpenFOAM patches)
    //- This strucutre is made to store all the data related to an interface
    struct InterfaceData
    {
        std::string nameOfInterface;
        std::vector<std::string> patchNames;
        std::vector<int> patchIDs;
        std::vector<std::string> importDataIdentifier;
        std::vector<std::string> exportDataIdentifier;
        std::vector<double> data_to_send;
        std::vector<double> data_to_recv;
        uint numNodes = 0;
        uint numElements = 0;
        std::vector<Node> Interface_nodes;                  // Contain Node objects
        std::vector<Element> Interface_elements;            // Contain elemental Objects
        std::vector<int> neighbour_ids_comm_num_of_nodes;   // Max size of a array =  TotalNumOfProcesses * 2
        int globalNodeIndexBegin = 0;
        int globalElementIndexBegin = 0;
    };

    //- Array of interface-strctures to store data for all coupling interfaces
    std::vector<struct InterfaceData> interfaces_;

    //- Array of CoSimIO::Model parts(to store mesh info) for all coupling interfaces
    std::vector<std::unique_ptr<CoSimIO::ModelPart>> model_part_interfaces_;

    //- For Force Field calculations
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

    // ********************************* Auxiliar functions for configuration ****************************************//
    //- To read all user sepcified input from controlDict
    bool readConfig(const dictionary& );

    //- To print the console ouput from this adapter
    void debugInfo(const std::string, const std::string );

    //- Determine the type of the solver (using pressure dimension)
    std::string determineSolverType();

    //- To compare two points
    bool is_same_points(Foam::vector& , Foam::vector& );

    //- To Compare the new node with all previous nodes before creating the new CoSimIO::Node
    int compare_nodes(Foam::vector& , std::size_t );


    // *********************************** Auxiliar functions for load calculation ************************************//
    //- Calculate Total Force on each element or cell in the interface mesh
    bool calculateForces(std::size_t);

    //- Stress tensor (see the OpenFOAM "Forces" function object)
    Foam::tmp<Foam::volSymmTensorField> devRhoReff() const;

    //- To get or read the value of density
    Foam::tmp<Foam::volScalarField> rho() const;

    //- To get or read the value of dynamic viscosity
    Foam::tmp<Foam::volScalarField> mu() const;

    //- To find the cell or face normal vectors
    Foam::vectorField getFaceVectors(const unsigned int ) const;

    //- Conversion of Elemental Force values to the nodal force vales
    void conversionElementalToNodalValues(std::size_t);


    // *********************************** Auxiliar functions for interaction with CoSimulation *************************//
    //- Connection with CoSimulation
    void connectKratos();

    //- DisConnect from CoSimulation
    void disconnectKratos();

    //- Export Mesh to CoSimulation
    void exportMeshToCosim();

    //- Resize the data vectors for communication with CoSimulation
    void resizeDataVectors();

    //- Export vector data to CoSimulation
    void exportDataToKratos();

    //- Import vector data from CoSimulation
    void importDataFromKratos();

public:

    //- Runtime type information
    TypeName("KratosOpenfoamAdapterFunctionObject");

    //- Construct from Time and dictionary
    KratosOpenfoamAdapterFunctionObject
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    //- Destructor
    virtual ~KratosOpenfoamAdapterFunctionObject();

    // ********************************* FunctionObject's Member Functions *********************************** //
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

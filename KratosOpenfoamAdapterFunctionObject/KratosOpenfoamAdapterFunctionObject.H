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
        Property            | Description                                       | Required | Default value
        type                | type name: KratosOpenfoamAdapterFunctionObject    | yes      |
        participant name    | Name of the participant                           | no       | Fluid
        readData            | Name of the parameters to Read by OpenFOAM        | yes      |
        writeData           | Name of the parameters to write by OpenFOAM       | yes      |
    \endtable

SourceFiles
    KratosOpenfoamAdapterFunctionObject.C

\*-----------------------------------------------------------------------*/
/*
#include <iostream>
using namespace std;
#include <iomanip> */

#ifndef KRATOSOPENFOAMADAPTERFUNCTIONOBJECT_H
#define KRATOSOPENFOAMADAPTERFUNCTIONOBJECT_H

// OpenFOAM header files
#include "fvMeshFunctionObject.H"

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

    //- Adapter's main class instance
    //- CoSimulationAdapter CoSimulationAdapter_;
    //- Member variables required for CoSimIO
    CoSimIO::Info connect_info;
    std::string connection_name;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        KratosOpenfoamAdapterFunctionObject(const KratosOpenfoamAdapterFunctionObject&);

        //- Disallow default bitwise assignment
        void operator=(const KratosOpenfoamAdapterFunctionObject&) = delete;

        //- Private function for CoSimIO

protected:
    //- Protected Data to set the parameters required for coupling

        //Parameters for Adapter's configuration
        Foam::word participantName_;

        //Vector of string to save the name of the data to be read by Fluid Solver-OpenFOAM
        Foam::wordList readData_;

        //Vector of string to save the name of the data to be write by Fluid Solver-OpenFOAM
        Foam::wordList writeData_;

public:
    //- Runtime type information
    TypeName("KratosOpenfoamAdapterFunctionObject");

    //Constructors

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

        //- Read the data/after each change
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write data after each loop??
        virtual bool write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

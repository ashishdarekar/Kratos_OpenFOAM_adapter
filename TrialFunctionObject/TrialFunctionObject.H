/*-----------------------------------------------------------------------*\

Master Thesis
Ashish Darekar

Class
    Foam::functionObjects::TrialFunctionObject

Description
    This function object prints "Hello World" on the terminal

    Example of function object specification:
    \verbatim
    Trial_
    {
        type           TrialFunctionObject;
        libs ("libfieldFunctionObject.so");
        ...
        wordData       someWord;
        scalarData     1.0;
        labelData      1;
    }
    \endverbatim

Usage
    \table
        Property     | Description               | Required | Default value
        type         | type name: FUNCTIONOBJECT | yes      |
        wordData     | some word option...       | no       | defaultWord
        scalarData   | some scalar value...      | yes      |
        labelData    | some label value...       | yes      |
    \endtable

SourceFiles
    TrialFunctionObject.C

\*---------------------------------------------------------------------------*/

#include <iostream>
using namespace std;
#include <iomanip>

#ifndef TRIALFUNCTIONOBJECT_H
#define TRIALFUNCTIONOBJECT_H

//OpenFOAM header files
#include "fvMeshFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class TRIALFUNCTIONOBJECT Declaration
\*---------------------------------------------------------------------------*/

class TrialFunctionObject
:
    public fvMeshFunctionObject
{
    // Private Data

    // Private Member Functions

    //- Disallow default bitwise copy construct

    //- Disallow default bitwise assignment
        void operator=(const TrialFunctionObject&) = delete;

public:

    //- Runtime type information
    TypeName("TrialFunctionObject");


    // Constructors

        //- Construct from Time and dictionary
        TrialFunctionObject
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

    //- Destructor
    virtual ~TrialFunctionObject();


    // Member Functions

        //- Read the FUNCTIONOBJECT data
        //virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        //virtual bool end();

        //- Write the TrialFunctionObject
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

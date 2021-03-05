/*-----------------------------------------------------------------------*\

Master Thesis
Ashish Darekar

\*---------------------------------------------------------------------------*/

#include "TrialFunctionObject.H"

// OpenFOAM header files
#include "Time.H"
#include "fvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(TrialFunctionObject, 0);
    addToRunTimeSelectionTable(functionObject, TrialFunctionObject, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::TrialFunctionObject::TrialFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::TrialFunctionObject::~TrialFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::functionObjects::TrialFunctionObject::execute()
{
    return true;
}

bool Foam::functionObjects::TrialFunctionObject::write()
{
    cout<<"Hello World - Ashish here\n";

    cout<<"Welcome to the Adapter\n";

    std::cout << std::setw(6) << "C++ version used: " <<  __cplusplus << std::endl;

    return true;
}



// ************************************************************************* //

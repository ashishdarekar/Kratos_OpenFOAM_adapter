/*-----------------------------------------------------------------------*\

Master-Thesis Work
Ashish Darekar

Sourcefile for the CoSimulationAdapter.H

\*-----------------------------------------------------------------------*/

#include "CoSimulationAdapter.H"

//using namespace Foam;

/* CoSimulationAdapter::CoSimulationAdapter(const Foam::Time& runTime, const Foam::fvMesh& mesh)
:
runTime_(runTime),
mesh_(mesh)
{
    std::cout << "CoSimulation Adapter is loaded" << std::endl;
    return;
}
 */

CoSimulationAdapter::CoSimulationAdapter()
{
    std::cout << "CoSimulation Adapter is loaded" << std::endl;
    return;
}

void CoSimulationAdapter::configure()
{
    std::cout << "CoSimulation Adapter : Configuration" << std::endl;
    return;
}

void CoSimulationAdapter::execute()
{
    std::cout << "CoSimulation Adapter : Execute" << std::endl;
    return;
}

void CoSimulationAdapter::end()
{
    std::cout << "CoSimulation Adapter : Execute" << std::endl;
    return;
}

CoSimulationAdapter::~CoSimulationAdapter()
{
    std::cout << "CoSimulation Adapter : Destructor" << std::endl;
    return;
}
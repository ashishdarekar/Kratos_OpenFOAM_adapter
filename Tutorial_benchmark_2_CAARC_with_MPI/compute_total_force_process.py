from KratosMultiphysics import *

## This proces sets the value of a vector variable component-by-component.
## In this case, the fixicity is given by the user and some of the components may not be fixed.

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeTotalForceProcess(Model, settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class ComputeTotalForceProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.output_file = open(settings["output_filename"].GetString(),'w')
        self.output_file.write("#Drag for group " + settings["model_part_name"].GetString() + "\n")
        self.output_file.write("#time RX RY RZ MX MY MZ\n")
        self.output_file.flush()

    def ExecuteInitialize(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        rx = 0.0
        ry = 0.0
        rz = 0.0

        mx = 0.0
        my = 0.0
        mz = 0.0

        for node in self.model_part.Nodes:
            reaction = node.GetSolutionStepValue(REACTION,0)
            rx += reaction[0]
            ry += reaction[1]
            rz += reaction[2]

            x = node.X - 0.0
            y = node.Y - 0.0
            z = node.Z - 0.0
            mx -= y * reaction[2] - z * reaction[1]
            my -= z * reaction[0] - x * reaction[2]
            mz -= x * reaction[1] - y * reaction[0]

        time = self.model_part.ProcessInfo[TIME]
        self.output_file.write(str(time) + " " + str(rx) + " " + str(ry) + " " + str(rz) + " " + str(mx) + " " + str(my) + " " + str(mz) + "\n")
        self.output_file.flush()

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass





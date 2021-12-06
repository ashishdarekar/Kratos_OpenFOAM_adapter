from KratosMultiphysics import *


# import to have as globals - object creation for each implemented inlet functions
from inlet_velocity_2D_function_implementations import *

## This proces sets the value of a vector variable component-by-component.
## In this case, the fixicity is given by the user and some of the components may not be fixed.

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeCustomInletVelocityProcess(Model, settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class ImposeCustomInletVelocityProcess(Process):
    def __init__(self, Model, settings ):
        Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.inlet_type = settings["inlet_type"].GetString()

        # check if directory exists if not create it
        file_path = settings["output_inlet_coordinates_filename"].GetString()
        directory = os.path.dirname(file_path)
        if not os.path.exists(directory):
            os.makedirs(directory)
        self.output_file = open(settings["output_inlet_coordinates_filename"].GetString(),'w')
        self.output_file.write("#Inlet nodes for group " + settings["model_part_name"].GetString() + " applying " + self.inlet_type + "\n")
        self.output_file.write("# Node-ID  CoordX  CoordY  CoordZ\n")

        for node in self.model_part.Nodes:
            self.output_file.write(str(node.Id) + " " + str(node.X) + " " + str(node.Y) + " " + str(node.Z) + "\n")
            self.output_file.flush()

        self.output_file.close()

        # check if directory exists if not create it
        file_path = settings["output_inlet_velocity_filename"].GetString()
        directory = os.path.dirname(file_path)
        if not os.path.exists(directory):
            os.makedirs(directory)
        self.output_file = open(settings["output_inlet_velocity_filename"].GetString(),'w')
        self.output_file.write("#Inlet velocity for group " + settings["model_part_name"].GetString() + " applying " + self.inlet_type + "\n")
        self.output_file.write("# time  nodal results - Vx - 1...n \n")
        self.output_file.flush()

        self.inlet_settings = settings

    def ExecuteInitialize(self):

        ##using default parameters, customize to use other
        if self.inlet_type == "Constant2D_Inlet":
            self.inlet_velocity = ConstantInletVelocity(self.model_part.Nodes,
                                                        self.inlet_settings["inlet_parameters"]["mean_velocity"].GetDouble(),
                                                        self.inlet_settings["ramp_up"].GetBool(),
                                                        self.inlet_settings["ramp_up_time"].GetDouble())

            print("Constant2D_Inlet is initialized")

        elif self.inlet_type == "ConstantParabola2D_Inlet":
            self.inlet_velocity = ParabolicInletVelocity(self.model_part.Nodes,
                                                         self.inlet_settings["inlet_parameters"]["mean_velocity_at_middle"].GetDouble(),
                                                         self.inlet_settings["ramp_up"].GetBool(),
                                                         self.inlet_settings["ramp_up_time"].GetDouble())
            print("ConstParabola2D Inlet is initialized")

        elif self.inlet_type == "OscillatingParabola2D_Inlet":
            self.inlet_velocity = OscillatingParabolicInletVelocity(self.model_part.Nodes,
                                                                    self.inlet_settings["inlet_parameters"]["mean_velocity_at_middle"].GetDouble(),
                                                                    self.inlet_settings["inlet_parameters"]["velocity_deviation"].GetDouble(),
                                                                    self.inlet_settings["inlet_parameters"]["deviation_frequency"].GetDouble(),
                                                                    self.inlet_settings["ramp_up"].GetBool(),
                                                                    self.inlet_settings["ramp_up_time"].GetDouble())
            print("OscillatingParabola2D is initialized")

        elif self.inlet_type == "PowerLawInletVelocity":
            self.inlet_velocity = PowerLawInletVelocity(self.model_part.Nodes,
                                                        self.inlet_settings["inlet_parameters"]["mean_velocity"].GetDouble(),
                                                        self.inlet_settings["inlet_parameters"]["reference_height_z"].GetDouble(),
                                                        self.inlet_settings["inlet_parameters"]["alpha"].GetDouble(),
                                                        self.inlet_settings["ramp_up"].GetBool(),
                                                        self.inlet_settings["ramp_up_time"].GetDouble())
            print("PowerLaw2D is initialized")

        elif self.inlet_type == "LogarithmicLaw2D_Inlet":
            self.inlet_velocity = LogLawInletVelocity(self.model_part.Nodes,
                                                      self.inlet_settings["inlet_parameters"]["mean_velocity"].GetDouble(),
                                                      self.inlet_settings["inlet_parameters"]["reference_height_z"].GetDouble(),
                                                      self.inlet_settings["ramp_up"].GetBool(),
                                                      self.inlet_settings["ramp_up_time"].GetDouble())
            print("LogarithmicLaw2D is initialized")

        elif self.inlet_type == "ConvolutedSineForPowerLaw2D_Inlet":
            self.inlet_velocity = ConvolutedSinOnPowerLawInletVelocity(self.model_part.Nodes,
                                                                       self.inlet_settings["inlet_parameters"]["mean_velocity"].GetDouble(),
                                                                       self.inlet_settings["inlet_parameters"]["reference_height_z"].GetDouble(),
                                                                       self.inlet_settings["inlet_parameters"]["alpha"].GetDouble(),
                                                                       self.inlet_settings["inlet_parameters"]["velocity_of_sine"].GetDouble(),
                                                                       self.inlet_settings["inlet_parameters"]["period_of_sine"].GetDouble(),
                                                                       self.inlet_settings["ramp_up"].GetBool(),
                                                                       self.inlet_settings["ramp_up_time"].GetDouble())
            print("ConvolutedSineForPowerLaw2D is initialized")

        elif self.inlet_type == "ConvolutedSineForLogarithmicLaw2D_Inlet":
            self.inlet_velocity = ConvolutedSinOnLogLawInletVelocity(self.model_part.Nodes,
                                                                     self.inlet_settings["inlet_parameters"]["mean_velocity"].GetDouble(),
                                                                     self.inlet_settings["inlet_parameters"]["reference_height_z"].GetDouble(),
                                                                     self.inlet_settings["inlet_parameters"]["velocity_of_sine"].GetDouble(),
                                                                     self.inlet_settings["inlet_parameters"]["period_of_sine"].GetDouble(),
                                                                     self.inlet_settings["ramp_up"].GetBool(),
                                                                     self.inlet_settings["ramp_up_time"].GetDouble())
            print("ConvolutedSineForLogarithmicLaw2D is initialized")

        elif self.inlet_type == "ConvolutedParabola2D_Inlet":
            self.inlet_velocity = ConvolutedSinOnParabolicInletVelocity(self.model_part.Nodes,
                                                                        self.inlet_settings["inlet_parameters"]["mean_velocity_at_middle"].GetDouble(),
                                                                        self.inlet_settings["inlet_parameters"]["velocity_deviation"].GetDouble(),
                                                                        self.inlet_settings["inlet_parameters"]["velocity_of_sine"].GetDouble(),
                                                                        self.inlet_settings["inlet_parameters"]["period_of_sine"].GetDouble(),
                                                                        self.inlet_settings["ramp_up"].GetBool(),
                                                                        self.inlet_settings["ramp_up_time"].GetDouble())
            print("ConvolutedParabole2D is initialized")

        elif self.inlet_type == "ConvolutedMultipliedParabola2D_Inlet":
            self.inlet_velocity = ConvolutedSinMultipliedOnParabolicInletVelocity(self.model_part.Nodes,
                                                                                  self.inlet_settings["inlet_parameters"]["mean_velocity_at_middle"].GetDouble(),
                                                                                  self.inlet_settings["inlet_parameters"]["velocity_deviation"].GetDouble(),
                                                                                  self.inlet_settings["inlet_parameters"]["velocity_of_sine"].GetDouble(),
                                                                                  self.inlet_settings["inlet_parameters"]["period_of_sine"].GetDouble(),
                                                                                  self.inlet_settings["ramp_up"].GetBool(),
                                                                                  self.inlet_settings["ramp_up_time"].GetDouble())
            print("ConvolutedMultipliedParabola2D is initialized")

    def ExecuteInitializeSolutionStep(self):

        ## applying custom velocity
        self.inlet_velocity.ApplyInletVelocity(self.model_part.ProcessInfo[TIME])


    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):

        ## printing out the applied velocity done here
        output_result_for_time_step = str(self.model_part.ProcessInfo[TIME]) + " "
        for node in self.model_part.Nodes:
            output_result_for_time_step += str(node.GetSolutionStepValue(VELOCITY_X, 0)) + " " # VELOCITY_X is velocity in stream direction
        output_result_for_time_step += "\n"

        self.output_file.write(output_result_for_time_step)
        self.output_file.flush()

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass

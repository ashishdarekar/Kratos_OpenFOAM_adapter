from KratosMultiphysics import *
from point_output_process import PointOutputProcess

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MultiplePointsOutputProcess(Model, settings["Parameters"])

class MultiplePointsOutputProcess(Process):

    defaults = Parameters('''{
        "locations"         : [],
        "model_part_name"  : "",
        "output_filename" : "",
        "output_variables" : []
    }''')

    def __init__(self,Model,settings):

        self.model_part = Model[settings["model_part_name"].GetString()]

        positions = []

        for i in range(0,settings["locations"].size() ):
            position=[]
            for j in range (0,settings["locations"][i].size()):
                position.append(settings["locations"][i][j].GetDouble())
            positions.append(position)

        # re-using the point_output_process
        print("Multiple points output - considered points")
        self.point_output_objects = []
        for idx,position in enumerate(positions):
            print(idx," ",position)

            variable_names = "["
            for i in range( 0,settings["output_variables"].size() ):
                variable_names += r'"' + settings["output_variables"][i].GetString() + r'"'
                if i < settings["output_variables"].size() -1:
                    variable_names += ","
            variable_names += "]"

            output_filename = settings["output_filename"].GetString() + "_" + str(idx) + settings["output_file_format"].GetString()

            # formatting for a mock JSON input string
            input_data_string = "{" + "\n"
            input_data_string += r'"' + "position" + r'"' + " : " + str(position) + "," + "\n"
            input_data_string += r'"' + "model_part_name" + r'"' + " : " + r'"' + settings["model_part_name"].GetString() + r'"' + "," + "\n"
            input_data_string += r'"' + "output_filename" + r'"' + " : " + r'"' + output_filename + r'"' +  "," + "\n"
            input_data_string += r'"' + "output_variables" + r'"' + " : " + variable_names + "\n"
            input_data_string += "}"

            # parsing the string to the needed JSON format
            point_parameters = Parameters(input_data_string)
            self.point_output_objects.append(PointOutputProcess(Model,point_parameters))

        #print(list(enumerate(positions)))
        #wait = input("check and press key 1...")

    def ExecuteInitialize(self):

        for idx,point_output_object in enumerate(self.point_output_objects):
            point_output_object.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):

        for idx,point_output_object in enumerate(self.point_output_objects):
            point_output_object.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):

        for idx,point_output_object in enumerate(self.point_output_objects):
            point_output_object.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):

        for idx,point_output_object in enumerate(self.point_output_objects):
            point_output_object.ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):

        for idx,point_output_object in enumerate(self.point_output_objects):
            point_output_object.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):

        for idx,point_output_object in enumerate(self.point_output_objects):
            point_output_object.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):

        for idx,point_output_object in enumerate(self.point_output_objects):
            point_output_object.ExecuteFinalize()


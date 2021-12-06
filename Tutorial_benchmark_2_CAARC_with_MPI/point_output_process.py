from KratosMultiphysics import *
import os


def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")
    return PointOutputProcess(Model, settings["Parameters"])


class PointOutputProcess(Process):

    defaults = Parameters('''{
        "position"         : [],
        "model_part_name"  : "",
        "output_filename" : "",
        "output_variables" : []
    }''')

    def __init__(self, Model, settings):

        # self.model_part = model_part
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.settings = settings
        self.positions = []
        self.elements = []
        self.area_coordinates = []
        self.output_files = []
        self.output_variables = []

    def ExecuteInitialize(self):

        settings = self.settings
        settings.ValidateAndAssignDefaults(self.defaults)

        # create directory if non existent
        file_path = settings["output_filename"].GetString()
        directory = os.path.dirname(file_path)

        if not os.path.exists(directory):
            os.makedirs(directory)

        output_filename = settings["output_filename"].GetString()

        # positiion parameters from line output process
        position_data = settings["position"]
        pos = []

        for i in range(0, position_data.size()):
            pos.append(position_data[i].GetDouble())

        position = Vector(pos)  # conversion from list to ublas vector

        # Identify the position of the point within the mesh
        area_coordinates = Vector()
        point_locator = PointLocation(self.model_part)
        try:
            elem_id = point_locator.Find(position, area_coordinates)
        except:
            e = sys.exc_info()[0]
            write_to_page( "<p>Error: %s</p>" % e )
                                     # If successful, fills in area coordinates
                                     # in second argument

        # Check if a point was found, and initalize output
        # NOTE: If the search was not successful, we fail silently and do nothing.
        # This is BY DESIGN, as we are supposed to work on MPI too, and the point
        # in question might lie on a different partition.

        if point_locator.found():

            self.positions.append(position)
            self.elements.append(self.model_part.Elements[
                                 elem_id])  # Note: this is actually a find
            self.area_coordinates.append(area_coordinates)

            variable_data = settings["output_variables"]
            if not variable_data.IsArray():
                raise Exception(
                    "{0} Error: Variable list is unreadable".format(self.__class__.__name__))

            variable_names = [variable_data[i].GetString()
                              for i in range(0, variable_data.size())]
            output_variables = [KratosGlobals.GetVariable(var)
                                for var in variable_names]

            output_file = self.__initialize_output_file(
                output_filename, position, variable_names, output_variables)

            self.output_files.append(output_file)
            self.output_variables.append(output_variables)

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        '''Print results to file.'''

        time = self.model_part.ProcessInfo[TIME]

        for var_list, elem, coord, f in zip(self.output_variables, self.elements, self.area_coordinates, self.output_files):
            out = str(time)
            for var in var_list:
                value = self.__interpolate(var, elem, coord)

                if self.__is_array_var(var):
                    out += " " + " ".join(str(v) for v in value)
                else:
                    out += " " + str(value)

            out += "\n"
            f.write(out)
            f.flush()

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        '''Close output files.'''
        for f in self.output_files:
            f.close()

    def __initialize_output_file(self, output_filename, position, variable_names, output_variables):
        try:
            output_file = open(output_filename, "w")
        except IOError:
            print("Wrong file or file path")
        # Note that zip stops once the shortest list is finished, which is
        # exactly what we want here
        out  = "# Results for position " + \
            " ".join(a + str(b)
                     for a, b in zip([" x: ", "; y: ", "; z: "], position)) + "\n"
        out += "# time"
        for var, varname in zip(output_variables, variable_names):
            # if this is a Variable< array_1d< double,3 > >
            if self.__is_array_var(var):
                out += " {0}_X {0}_Y {0}_Z".format(varname)
            else:
                out += " " + varname

        out += "\n"
        output_file.write(out)

        return output_file

    def __interpolate(self, variable, element, coordinates):

        nodes = element.GetNodes()
        # Initializing 'value' like this, i don't need to know its type
        value = nodes[0].GetSolutionStepValue(variable) * coordinates[0]
        for n, c in zip(nodes[1:], coordinates[1:]):
            value = value + c * n.GetSolutionStepValue(variable)

        return value

    def __is_array_var(self, var):
        return type(var) == type(VELOCITY)

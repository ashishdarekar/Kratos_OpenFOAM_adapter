import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SM

from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
import math


def Factory(params, Model):
    if(type(params) != KratosMultiphysics.Parameters):
        raise Exception(
            'expected input shall be a Parameters object, encapsulating a json string')
    return OutputSumForcesProcess(Model, params['Parameters'])


class OutputSumForcesProcess(KratosMultiphysics.Process):
    '''
    TODO
    '''
    def __init__(self, Model, params):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"       : "",
                "interval"              : [0.0, 1e30],                
                "print_to_screen"       : false,
                "print_format"          : ".8f",
                "write_output_file"     : true,
                "output_file_settings"  : {}
            }
            """)

        # Detect 'End' as a tag and replace it by a large number
        if(params.Has('interval')):
            if(params['interval'][1].IsString()):
                if(params['interval'][1].GetString() == 'End'):
                    params['interval'][1].SetDouble(1e30)
                else:
                    raise Exception('The second value of interval can be \'End\' or a number, interval currently:' +
                                    params['interval'].PrettyPrintJsonString())

        params.ValidateAndAssignDefaults(default_settings)

        self.model_part_name = params['model_part_name'].GetString()
        self.model_part = Model[self.model_part_name]
        self.interval = params["interval"].GetVector()
        self.print_to_screen = params['print_to_screen'].GetBool()
        self.write_output_file = params['write_output_file'].GetBool()
        self.format = params["print_format"].GetString()

        if (self.model_part.GetCommunicator().MyPID() == 0):
            if (self.write_output_file):

                output_file_name = params["model_part_name"].GetString(
                ) + "_force.dat"

                file_handler_params = KratosMultiphysics.Parameters(
                    params["output_file_settings"])

                if file_handler_params.Has("file_name"):
                    warn_msg = 'Unexpected user-specified entry found in "output_file_settings": {"file_name": '
                    warn_msg += '"' + \
                        file_handler_params["file_name"].GetString() + '"}\n'
                    warn_msg += 'Using this specififed file name instead of the default "' + \
                        output_file_name + '"'
                    KratosMultiphysics.Logger.PrintWarning(
                        "OutputSumForcesProcess", warn_msg)
                else:
                    file_handler_params.AddEmptyValue("file_name")
                    file_handler_params["file_name"].SetString(
                        output_file_name)

                file_header = self._GetFileHeader()
                self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part,
                                                                   file_handler_params, file_header).file

    def ExecuteFinalizeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if((current_time >= self.interval[0]) and (current_time < self.interval[1])):
            results = self._EvaluateGlobalForces()

            if (self.model_part.GetCommunicator().MyPID() == 0):
               
                output = []
                output.extend(results)               
                output_vals = [format(val, self.format) for val in output]
                # not formatting time in order to not lead to problems with time recognition
                # in the file writer when restarting               
                output_vals.insert(0, str(current_time))                

                res_labels = ['time: ','sum_forces_x: ', 'sum_forces_y: ', 'sum_forces_z: ']

                if (self.print_to_screen):
                    result_msg = 'Force evaluation for model part ' + \
                        self.model_part_name + '\n'
                    result_msg += ', '.join([a+b for a,
                                             b in zip(res_labels, output_vals)])
                    self._PrintToScreen(result_msg)

                if (self.write_output_file):
                    self.output_file.write(' '.join(output_vals) + '\n')

    def _EvaluateGlobalForces(self):
        forces = [0.0,0.0,0.0]        

        for node in self.model_part.Nodes:
            forces[0] += node.GetSolutionStepValue(SM.POINT_LOAD_X,0)
            forces[1] += node.GetSolutionStepValue(SM.POINT_LOAD_Y,0)
            forces[2] += node.GetSolutionStepValue(SM.POINT_LOAD_Z,0)
                
        forces = self.model_part.GetCommunicator().GetDataCommunicator().SumDoubles(forces,0)

        return forces
                      

    def _GetFileHeader(self):
        header = '# Forces for model part ' + self.model_part_name + '\n'
        header += '# Time sum_forces_x sum_forces_y sum_forces_z \n'
        return header

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo(
            'OutputSumForcesProcess', 'Force results: ')
        KratosMultiphysics.Logger.PrintInfo(
            'OutputSumForcesProcess', 'Current time: ' + result_msg)

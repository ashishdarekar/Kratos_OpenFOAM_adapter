#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS19-20 
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek
        
Author: mate.pentek@tum.de 

Description: Kratos global forces in body- and flow-attached coordinates

Created on:  06.01.2020
Last update: 06.01.2020
'''
#===============================================================================

import KratosMultiphysics
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
import sys

import math

def ccw_rotate_comp_around_z(vect, angle):
    '''
    Counter-clockwise roration around z-axis in a right-handed coordinate system
    '''
    rot_matr_z = [[math.cos(angle), math.sin(angle), 0],
                  [-math.sin(angle), math.cos(angle), 0],
                  [0, 0, 1]]

    vect_rot = [0.0, 0.0, 0.0]

    for i in range(3):
        vect_rot[i] = sum([a*b for a, b in zip(rot_matr_z[i], vect)])

    return vect_rot


def Factory(params, Model):
    if(type(params) != KratosMultiphysics.Parameters):
        raise Exception(
            'expected input shall be a Parameters object, encapsulating a json string')
    return ComputeBaseReactionProcess(Model, params['Parameters'])


class ComputeBaseReactionProcess(KratosMultiphysics.Process):
    '''
    Computes the base reactopm
    '''
    def __init__(self, Model, params):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name_list"  : [""],
                "interval"              : [0.0, 1e30],
                "reference_point"       : [0.0,0.0,0.0],
                "z_rotation_angle"      : 0.0,
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

        self.model_part_name_list = params['model_part_name_list'].GetStringArray()

        self.model_part_name = self.model_part_name_list[0]

        self.model_part_list = []
        for model_part_name in self.model_part_name_list:
            self.model_part_list.append( Model[model_part_name] )

        self.model_part = self.model_part_list[0]


        self.interval = params["interval"].GetVector()
        self.print_to_screen = params['print_to_screen'].GetBool()
        self.write_output_file = params['write_output_file'].GetBool()
        self.format = params["print_format"].GetString()

        # added reference point for moment calculation
        self.reference = params['reference_point'].GetVector()
        if self.reference.Size() != 3:
            raise Exception(
                'The reference point position has to be provided with 3 coordinates!')

        # user inpput expected in degrees, here changing to radians
        self.theta = math.radians(params['z_rotation_angle'].GetDouble())

        if (self.model_part.GetCommunicator().MyPID() == 0):
            if (self.write_output_file):

                output_file_name = self.model_part_name + "_base_reaction.dat"

                file_handler_params = KratosMultiphysics.Parameters(
                    params["output_file_settings"])

                if file_handler_params.Has("file_name"):
                    warn_msg = 'Unexpected user-specified entry found in "output_file_settings": {"file_name": '
                    warn_msg += '"' + \
                        file_handler_params["file_name"].GetString() + '"}\n'
                    warn_msg += 'Using this specififed file name instead of the default "' + \
                        output_file_name + '"'
                    KratosMultiphysics.Logger.PrintWarning(
                        "ComputeBaseReactionProcess", warn_msg)
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
            ff, mf, fb, mb = self._EvaluateGlobalForces()

            if (self.model_part.GetCommunicator().MyPID() == 0):
                output = []
                output.extend(ff)
                output.extend(mf)
                output.extend(fb)
                output.extend(mb)

                output_vals = [format(val, self.format) for val in output]
                # not formatting time in order to not lead to problems with time recognition
                # in the file writer when restarting
                output_vals.insert(0, str(current_time))

                res_labels = ['time: ',
                              'fx: ', 'fy: ', 'fz: ', 'mx: ', 'my: ', 'mz: ',
                              'fx\': ', 'fy\': ', 'fz\': ', 'mx\': ', 'my\': ', 'mz\': ']

                if (self.print_to_screen):
                    result_msg = 'Global force evaluation for model part ' + \
                        self.model_part_name + '\n'
                    result_msg += ', '.join([a+b for a,
                                             b in zip(res_labels, output_vals)])
                    self._PrintToScreen(result_msg)
                    sys.stdout.flush()

                if (self.write_output_file):
                    self.output_file.write(' '.join(output_vals) + '\n')

    def _EvaluateGlobalForces(self):
        # flow-attached forces: in x-y-z coordinate system
        ff = [0.0, 0.0, 0.0]
        mf = [0.0, 0.0, 0.0]

        for model_part_i in self.model_part_list:

            for node in model_part_i.GetCommunicator().LocalMesh().Nodes:
                # sign is flipped to go from reaction to action -> force
                reaction = node.GetSolutionStepValue(KratosMultiphysics.REACTION, 0)
                #moment_reaction = node.GetSolutionStepValue(KratosMultiphysics.REACTION_MOMENT, 0)

                # summing up nodal contributions to get resultant for model_part
                ff[0] += (-1) * reaction[0]
                ff[1] += (-1) * reaction[1]
                ff[2] += (-1) * reaction[2]

                x = node.X - self.reference[0]
                y = node.Y - self.reference[1]
                z = node.Z - self.reference[2]
                mf[0] += y * (-1) * reaction[2] - z * (-1) * reaction[1] #+ (-1) * moment_reaction[0]
                mf[1] += z * (-1) * reaction[0] - x * (-1) * reaction[2] #+ (-1) * moment_reaction[1]
                mf[2] += x * (-1) * reaction[1] - y * (-1) * reaction[0] #+ (-1) * moment_reaction[2]

        # body-attached forces -> here only a rotation around z-axis
        # in x'-y'-z' coordinate system
        # of the summed-up forces and moment

        fb = ccw_rotate_comp_around_z(ff, self.theta)
        mb = ccw_rotate_comp_around_z(mf, self.theta)
             
        ff = self.model_part.GetCommunicator().GetDataCommunicator().SumDoubles(ff,0)
        mf = self.model_part.GetCommunicator().GetDataCommunicator().SumDoubles(mf,0)

        fb = self.model_part.GetCommunicator().GetDataCommunicator().SumDoubles(fb,0)
        mb = self.model_part.GetCommunicator().GetDataCommunicator().SumDoubles(mb,0)
        
        return ff, mf, fb, mb

    def _GetFileHeader(self):
        header = '# Base reaction for model part ' + self.model_part_name + '\n'
        header += '# Time Fx Fy Fz Mx My Mz Fx\' Fy\' Fz\' Mx\' My\' Mz\'\n'
        return header

    def _PrintToScreen(self, result_msg):
        KratosMultiphysics.Logger.PrintInfo(
            'ComputeBaseReactionProcess', 'Base reaction results:')
        KratosMultiphysics.Logger.PrintInfo(
            'ComputeBaseReactionProcess', 'Current time: ' + result_msg)

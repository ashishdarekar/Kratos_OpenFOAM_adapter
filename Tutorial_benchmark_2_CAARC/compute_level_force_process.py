#===============================================================================
'''
Project:Lecture - Structural Wind Engineering WS19-20
        Chair of Structural Analysis @ TUM - R. Wuchner, M. Pentek

Author: mate.pentek@tum.de

Description: Kratos level forces in body- and flow-attached coordinates

Created on:  06.01.2020
Last update: 06.01.2020
'''
#===============================================================================

import KratosMultiphysics
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility
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
    return ComputeLevelForceProcess(Model, params["Parameters"])


class ComputeLevelForceProcess(KratosMultiphysics.Process):
    '''
    Computes the flow- and body-attached forces
    for a model part with body-fitted mesh
    split into a number of intervals
    thus the naming LevelForces

    Takes as input a CCW positive rotation (in degrees) around
    axis z for the body-attached forces
    '''
    def __init__(self, Model, params):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"       : "",
                "interval"              : [0.0, 1e30],
                "z_rotation_angle"      : 0.0,
                "start_point"           : [],
                "end_point"             : [],
                "intervals"             : 1,
                "include_endpoints"     : true,
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

        # modal part params
        self.model_part_name = params['model_part_name'].GetString()
        self.model_part = Model[self.model_part_name]
        self.interval = params["interval"].GetVector()
        self.print_to_screen = params['print_to_screen'].GetBool()
        self.write_output_file = params['write_output_file'].GetBool()
        self.format = params["print_format"].GetString()

        # user inpput expected in degrees, here changing to radians
        self.theta = math.radians(params['z_rotation_angle'].GetDouble())

        # creating parametrized and final start-center-end point for intervals
        start_point_position = params["start_point"].GetVector()
        if start_point_position.Size() != 3:
            raise Exception(
                'The start point position has to be provided with 3 coordinates!')
        end_point_position = params["end_point"].GetVector()
        if end_point_position.Size() != 3:
            raise Exception(
                'The end point position has to be provided with 3 coordinates!')
        number_of_sampling_intervals = params["intervals"].GetInt()

        if number_of_sampling_intervals <= 0:
            raise Exception(
                'The number of sampling points has to be larger than 0!')
        else:
            include_endpoints = params["include_endpoints"].GetBool()

            # setup the parametric space for the internal points on the line
            lower_bound = 0.0
            upper_bound = 1.0

            my_param = [1.0] * (number_of_sampling_intervals)
            if include_endpoints:

                my_param[0] = 0.5
                my_param[-1] = 0.5

                parametrized_internal_points = [lower_bound + x*(upper_bound-lower_bound)/(
                    number_of_sampling_intervals-1) for x in range(number_of_sampling_intervals + 1)]

            else:

                parametrized_internal_points = [lower_bound + x*(upper_bound-lower_bound)/(
                    number_of_sampling_intervals-1) for x in range(number_of_sampling_intervals)]

            parametrized_internal_points = [0.0] * (len(my_param)+1)

            for i in range(len(my_param)):
                for j in range(i+1):
                    parametrized_internal_points[i+1] += my_param[j]

            parametrized_internal_points = [
                x / parametrized_internal_points[-1] for x in parametrized_internal_points]

            # determining the positions of the output points
            direction_vector = [
                x - y for x, y in zip(end_point_position, start_point_position)]

            current_positions = []
            for k in range(len(parametrized_internal_points)):
                current_positions.append(
                    [x + parametrized_internal_points[k]*y for x, y in zip(start_point_position, direction_vector)])

            self.levels = {}
            for idx in range(len(current_positions)-1):
                self.levels[idx] = {}
                self.levels[idx]['start'] = current_positions[idx]
                self.levels[idx]['end'] = current_positions[idx+1]

                # first interval
                if (include_endpoints and idx == 0):
                    self.levels[idx]['center'] = [
                        x1 for x1 in self.levels[idx]['start']]
                elif (include_endpoints and idx == (len(current_positions)-1)-1):
                    self.levels[idx]['center'] = [
                        x2 for x2 in self.levels[idx]['end']]
                else:
                    self.levels[idx]['center'] = [
                        (x1+x2)/2 for x1, x2 in zip(self.levels[idx]['start'], self.levels[idx]['end'])]

                self.levels[idx]['output_file'] = None

                for idx in range(len(self.levels)):
                    self.levels[idx]['node_ids'] = []
                    for node in self.model_part.GetCommunicator().LocalMesh().Nodes:
                        # PMT here only comparing Z coordinate
                        if idx == 0:
                            if self.levels[idx]['start'][2] <= node.Z0 and node.Z0 <= self.levels[idx]['end'][2]:
                                self.levels[idx]['node_ids'].append(node.Id)
                        else:
                            if self.levels[idx]['start'][2] < node.Z0 and node.Z0 <= self.levels[idx]['end'][2]:
                                self.levels[idx]['node_ids'].append(node.Id)

        # initialize output files for each level
        if (self.model_part.GetCommunicator().MyPID() == 0):
            if (self.write_output_file):

                # create/check/assign file name prefix
                output_file_name_prefix = params["model_part_name"].GetString() + "_level_force_"

                file_handler_params = KratosMultiphysics.Parameters(
                    params["output_file_settings"])

                if file_handler_params.Has("file_name"):
                    warn_msg = 'Unexpected user-specified entry found in "output_file_settings": {"file_name": '
                    warn_msg += '"' + \
                        file_handler_params["file_name"].GetString(
                        ) + '"}\n'
                    warn_msg += 'Using this specififed file name instead of the default "' + \
                        output_file_name_prefix + '"'
                    KratosMultiphysics.Logger.PrintWarning(
                        "ComputeLevelForceProcess", warn_msg)

                    output_file_name_prefix = file_handler_params["file_name"].GetString() + "_level_force_"
                else:
                    file_handler_params.AddEmptyValue("file_name")

                for idx in range(len(self.levels)):
                    # file for each level
                    file_handler_params["file_name"].SetString(
                            output_file_name_prefix + str(idx) + '.dat')
                    file_header = self._GetFileHeader(idx)

                    self.levels[idx]['output_file'] = TimeBasedAsciiFileWriterUtility(self.model_part,
                                                                                      file_handler_params, file_header).file

    def ExecuteFinalizeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if((current_time >= self.interval[0]) and (current_time < self.interval[1])):
            for idx in range(len(self.levels)):
                ff, mf, fb, mb = self._EvaluateLevelForces(idx)

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
                        result_msg = 'Level force evaluation for model part ' + \
                            self.model_part_name + '\n' \
                            + ' and level ' + str(idx) + '\n'
                        result_msg += ', '.join([a+b for a,
                                                 b in zip(res_labels, output_vals)])
                        self._PrintToScreen(result_msg, idx)

                    if (self.write_output_file):
                        self.levels[idx]['output_file'].write(
                            ' '.join(output_vals) + '\n')

                        # NOTE: forcing flush
                        # check in TimeBasedAsciiFileWriterUtility why this is not handled properly
                        self.levels[idx]['output_file'].flush()

    def ExecuteFinalize(self):
        '''Close output files.'''
        if (self.model_part.GetCommunicator().MyPID() == 0):
            for idx in range(len(self.levels)):
                self.levels[idx]['output_file'].close()

    def _EvaluateLevelForces(self, idx):
        # flow-attached forces: in x-y-z coordinate system
        ff = [0.0, 0.0, 0.0]
        mf = [0.0, 0.0, 0.0]

        for node_id in self.levels[idx]['node_ids']:
            node = self.model_part.Nodes[node_id]

            # sign is flipped to go from reaction to action -> force
            nodal_force = (-1) * node.GetSolutionStepValue(KratosMultiphysics.REACTION, 0)

            # summing up nodal contributions to get resultant for model_part
            ff[0] += nodal_force[0]
            ff[1] += nodal_force[1]
            ff[2] += nodal_force[2]

            x = node.X - self.levels[idx]['center'][0]
            y = node.Y - self.levels[idx]['center'][1]
            z = node.Z - self.levels[idx]['center'][2]
            mf[0] += y * nodal_force[2] - z * nodal_force[1]
            mf[1] += z * nodal_force[0] - x * nodal_force[2]
            mf[2] += x * nodal_force[1] - y * nodal_force[0]

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

    def _GetFileHeader(self, idx):
        header = "# Level force for level " + str(idx) + "\n"
        header += "# start: " + ', '.join(str(coord)
                                          for coord in self.levels[idx]['start']) + "\n"
        header += "# center: " + ', '.join(str(coord)
                                           for coord in self.levels[idx]['center']) + "\n"
        header += "# end: " + ', '.join(str(coord)
                                        for coord in self.levels[idx]['end']) + "\n"
        header += "# as part of " + \
            self.model_part_name + "\n"
        header += '# Time Fx Fy Fz Mx My Mz Fx\' Fy\' Fz\' Mx\' My\' Mz\'\n'

        return header

    def _PrintToScreen(self, result_msg, idx):
        KratosMultiphysics.Logger.PrintInfo(
            'ComputeLevelForceProcess', 'Level ' + str(idx) + ' - flow- and body-attached:')
        KratosMultiphysics.Logger.PrintInfo(
            'ComputeLevelForceProcess', 'Current time: ' + result_msg)

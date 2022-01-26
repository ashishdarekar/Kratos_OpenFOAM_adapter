#=========================================================================
'''
Project:Lecture - Structural Wind Engineering WS15-16
        Chair of Structural Analysis @ TUM - A. Michalski, R. Wuchner, M. Pentek


Author: nina.korshunova@tum.de, mate.pentek@tum.de

Description: A simplified script for inlet profiles



Created on:  05.12.2015
Last update: 10.12.2015
'''
#=========================================================================


from math import pow, cos, sin, pi, log, e
import sys
import os
# sys.path.append(ProjectParameters.kratos_path)
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *


def RampUpFunction(amplitude, time, timeRampUp):
    if time < timeRampUp:
        newAmplitude =  0.5 * amplitude * \
            (1 + (cos(pi * time / timeRampUp - pi)))
    else:
        newAmplitude = amplitude
    return newAmplitude


class ConstantInletVelocity:

    def __init__(self, inletNodes, inletVelocity, Ramp, timeOfRampUp):
        self.inletNodes = inletNodes
        self.inletVelocity = inletVelocity
        self.Ramp = Ramp
        self.timeRampUp = timeOfRampUp

    def ApplyInletVelocity(self, t):
        if self.Ramp == True:
            Velocity = RampUpFunction(self.inletVelocity, t, self.timeRampUp)
        else:
            Velocity = self.inletVelocity
        for node in self.inletNodes:
            U = Velocity
            node.SetSolutionStepValue(VELOCITY_X, 0, U)
            node.SetSolutionStepValue(VELOCITY_Y, 0, 0.0)


class ParabolicInletVelocity:

    def __init__(self, inletNodes, inletVelocity, RampUp, timeOfRampUp):
        self.inletNodes = inletNodes
        self.inletVelocity = inletVelocity
        self.Ramp = RampUp
        self.timeRampUp = timeOfRampUp
        yMin = 0
        yMax = yMin
        for node in self.inletNodes:
            if(node.Y > yMax):
                yMax = node.Y
            elif(node.Y < yMin):
                yMin = node.Y
        self.yMax = yMax
        self.yMin = yMin

    def ApplyInletVelocity(self, t):
        if self.Ramp == True:
            Velocity = RampUpFunction(self.inletVelocity, t, self.timeRampUp)
        else:
            Velocity = self.inletVelocity
        for node in self.inletNodes:
            Y = node.Y
            U = (Y - self.yMax) * (Y - self.yMin) * \
                Velocity / (self.yMin * self.yMax)
            node.SetSolutionStepValue(VELOCITY_X, 0, U)
            node.SetSolutionStepValue(VELOCITY_Y, 0, 0.0)


class OscillatingParabolicInletVelocity:

    def __init__(self, inletNodes, vAverage, vDeviation, freq, RampUp, timeOfRampUp):
        self.inletNodes = inletNodes
        self.vMax = vAverage + vDeviation
        self.vMin = vAverage - vDeviation
        self.vAvg = vAverage
        self.Frequency = freq
        yMin = 0
        yMax = yMin
        for node in self.inletNodes:
            if(node.Y > yMax):
                yMax = node.Y
            elif(node.Y < yMin):
                yMin = node.Y
        self.yMax = yMax
        self.yMin = yMin
        self.vDeviation = vDeviation
        self.Ramp = RampUp
        self.timeRampUp = timeOfRampUp

    def ApplyInletVelocity(self, t):
        if self.Ramp == True:
            vAvg = RampUpFunction(self.vAvg, t, self.timeRampUp)
            vMax = self.vDeviation + vAvg
            vMin = vAvg - self.vDeviation
        else:
            vMax = self.vMax
            vMin = self.vMin
            vAvg = self.vAvg
        # Unit Frequency components
        f = cos(t * self.Frequency * 2 * pi)
            # Magnitude
        f = f * (vMax - vAvg) + vAvg
        for node in self.inletNodes:
            Y = node.Y
            # Unit Parabolic Base
            NormalParabola = (Y - self.yMax) * (
                Y - self.yMin) / (self.yMin * self.yMax)
            # Compile
            U = NormalParabola * f
            node.SetSolutionStepValue(VELOCITY_X, 0, U)
            node.SetSolutionStepValue(VELOCITY_Y, 0, 0.0)


class PowerLawInletVelocity:

    def __init__(self, inletNodes, vAverage, zRef, alpha, RampUp, timeOfRampUp):
        self.inletNodes = inletNodes
        self.vAverage = vAverage
        self.zRef = zRef
        self.alpha = alpha
        zMin = 0
        zMax = zMin
        for node in self.inletNodes:
            if(node.Z > zMax):
                zMax = node.Z
            elif(node.Z < zMin):
                zMin = node.Z
        self.zMax = zMax
        self.zMin = zMin
        self.Ramp = RampUp
        self.timeRampUp = timeOfRampUp

    def ApplyInletVelocity(self, t):
        if self.Ramp == True:
            Velocity = RampUpFunction(self.vAverage, t, self.timeRampUp)
        else:
            Velocity = self.vAverage
        for node in self.inletNodes:
            Z = node.Z
            U = Velocity * (pow(((Z - self.zMin) / self.zRef), self.alpha))
            node.SetSolutionStepValue(VELOCITY_X, 0, U)
            node.SetSolutionStepValue(VELOCITY_Y, 0, 0.0)
            node.SetSolutionStepValue(VELOCITY_Z, 0, 0.0)


class LogLawInletVelocity:

    def __init__(self, inletNodes, vAverage, zRef, RampUp, timeOfRampUp):
        self.inletNodes = inletNodes
        self.vAverage = vAverage
        self.zRef = zRef
        zMin = 0
        zMax = zMin
        for node in self.inletNodes:
            if(node.Z > zMax):
                zMax = node.Z
            elif(node.Z < zMin):
                zMin = node.Z
        self.zMax = zMax
        self.zMin = zMin
        self.Ramp = RampUp
        self.timeRampUp = timeOfRampUp

    def ApplyInletVelocity(self, t):
        if self.Ramp == True:
            Velocity = RampUpFunction(self.vAverage, t, self.timeRampUp)
        else:
            Velocity = self.vAverage
        for node in self.inletNodes:
            Z = node.Z
            temp = (Z - self.zMin) / self.zRef
            if temp == 0:
                U = Velocity / 0.41
            else:
                U = Velocity / 0.41 * (log(temp))
            node.SetSolutionStepValue(VELOCITY_X, 0, U)
            node.SetSolutionStepValue(VELOCITY_Y, 0, 0.0)
            node.SetSolutionStepValue(VELOCITY_Z, 0, 0.0)


class ConvolutedSinOnPowerLawInletVelocity:

    def __init__(self, inletNodes, vAverage, zRef, alpha, v_sin, X, RampUp, TimeRampUp):
        self.inletNodes = inletNodes
        self.vAverage = vAverage
        self.zRef = zRef
        self.alpha = alpha
        self.v_sin = v_sin
        self.Period = X

        yMin = 0
        yMax = yMin
        for node in self.inletNodes:
            if(node.Y > yMax):
                yMax = node.Y
            elif(node.Y < yMin):
                yMin = node.Y
        self.yMax = yMax
        self.yMin = yMin

        self.inf1 = yMin
        self.dir1 = 1
        self.inf2 = yMax - self.Period
        self.dir2 = 0

        self.tprev = 0
        self.Ramp = RampUp
        self.timeRampUp = TimeRampUp

    def ApplyInletVelocity(self, t):
        if self.Ramp == True:
            Velocity = RampUpFunction(self.vAverage, t, self.timeRampUp)
        else:
            Velocity = self.vAverage
        dt = t - self.tprev
        # Track influence region
        # (Region1)
        if self.dir1 == 1:
            self.inf1 = self.inf1 + dt * self.v_sin
        elif self.dir1 == 0:
            self.inf1 = self.inf1 - dt * self.v_sin
        if self.inf1 < self.yMin:
            self.dir1 = 1
            self.inf1 = self.inf1 + 2 * (self.yMin - self.inf1)
        elif self.inf1 + self.Period > self.yMax:
            self.dir1 = 0
            self.inf1 = self.inf1 - 2 * ((self.inf1 + self.Period) - self.yMax)
        # (Region2)
        if self.dir2 == 1:
            self.inf2 = self.inf2 + dt * self.v_sin
        elif self.dir2 == 0:
            self.inf2 = self.inf2 - dt * self.v_sin
        if self.inf2 < self.yMin:
            self.dir2 = 1
            self.inf2 = self.inf2 + 2 * (self.yMin - self.inf2)
        elif self.inf2 + self.Period > self.yMax:
            self.dir2 = 0
            self.inf2 = self.inf2 - 2 * ((self.inf2 + self.Period) - self.yMax)
        for node in self.inletNodes:
            Y = node.Y
            temp = ((Y - self.yMin) / self.zRef)
            temp2 = Velocity * pow(temp, self.alpha)
            # Add sine if inside area of influence
            if Y >= self.inf1 and Y <= self.inf1 + self.Period:
                S1 = sin((Y - self.inf1) / self.Period * 2 * pi)
            else:
                S1 = 0
            if Y >= self.inf2 and Y <= self.inf2 + self.Period:
                S2 = -sin((Y - self.inf2) / self.Period * 2 * pi)
            else:
                S2 = 0
            S = (S1 + S2) / 2
            U = temp2 + S
            node.SetSolutionStepValue(VELOCITY_X, 0, U)
            node.SetSolutionStepValue(VELOCITY_Y, 0, 0.0)
            self.tprev = t


class ConvolutedSinOnLogLawInletVelocity:

    def __init__(self, inletNodes, vAverage, zRef, v_sin, X, RampUp, TimeRampUp):
        self.inletNodes = inletNodes
        self.vAverage = vAverage
        self.zRef = zRef
        self.v_sin = v_sin
        self.Period = X

        yMin = 0
        yMax = yMin
        for node in self.inletNodes:
            if(node.Y > yMax):
                yMax = node.Y
            elif(node.Y < yMin):
                yMin = node.Y
        self.yMax = yMax
        self.yMin = yMin

        self.inf1 = yMin
        self.dir1 = 1
        self.inf2 = yMax - self.Period
        self.dir2 = 0

        self.tprev = 0
        self.Ramp = RampUp
        self.timeRampUp = TimeRampUp

    def ApplyInletVelocity(self, t):
        if self.Ramp == True:
            Velocity = RampUpFunction(self.vAverage, t, self.timeRampUp)
        else:
            Velocity = self.vAverage
        dt = t - self.tprev
        # Track influence region
        # (Region1)
        if self.dir1 == 1:
            self.inf1 = self.inf1 + dt * self.v_sin
        elif self.dir1 == 0:
            self.inf1 = self.inf1 - dt * self.v_sin
        if self.inf1 < self.yMin:
            self.dir1 = 1
            self.inf1 = self.inf1 + 2 * (self.yMin - self.inf1)
        elif self.inf1 + self.Period > self.yMax:
            self.dir1 = 0
            self.inf1 = self.inf1 - 2 * ((self.inf1 + self.Period) - self.yMax)
        # (Region2)
        if self.dir2 == 1:
            self.inf2 = self.inf2 + dt * self.v_sin
        elif self.dir2 == 0:
            self.inf2 = self.inf2 - dt * self.v_sin
        if self.inf2 < self.yMin:
            self.dir2 = 1
            self.inf2 = self.inf2 + 2 * (self.yMin - self.inf2)
        elif self.inf2 + self.Period > self.yMax:
            self.dir2 = 0
            self.inf2 = self.inf2 - 2 * ((self.inf2 + self.Period) - self.yMax)
        for node in self.inletNodes:
            Y = node.Y
            temp = (Y - self.yMin) / self.zRef
            if temp == 0:
                temp2 = Velocity / 0.41
            else:
                temp2 = Velocity / 0.41 * (log(temp))
            # Add sine if inside area of influence
            if Y >= self.inf1 and Y <= self.inf1 + self.Period:
                S1 = sin((Y - self.inf1) / self.Period * 2 * pi)
            else:
                S1 = 0
            if Y >= self.inf2 and Y <= self.inf2 + self.Period:
                S2 = -sin((Y - self.inf2) / self.Period * 2 * pi)
            else:
                S2 = 0
            S = (S1 + S2) / 2
            U = temp2 + S
            node.SetSolutionStepValue(VELOCITY_X, 0, U)
            node.SetSolutionStepValue(VELOCITY_Y, 0, 0.0)
            self.tprev = t


class ConvolutedSinOnParabolicInletVelocity:

    def __init__(self, inletNodes, vAverage, vDeviation, v_sin, X, RampUp, timeOfRampUp):
        self.inletNodes = inletNodes
        self.vMax = vAverage + vDeviation
        self.vMin = vAverage - vDeviation
        self.vAvg = vAverage
        self.v_sin = v_sin
        self.Period = X
        self.vDeviation = vDeviation

        yMin = 0
        yMax = yMin
        for node in self.inletNodes:
            if(node.Y > yMax):
                yMax = node.Y
            elif(node.Y < yMin):
                yMin = node.Y

        self.yMax = yMax
        self.yMin = yMin

        self.inf1 = yMin
        self.dir1 = 1
        self.inf2 = yMax - self.Period
        self.dir2 = 0
        self.tprev = 0

        self.Ramp = RampUp
        self.timeRampUp = timeOfRampUp

    def ApplyInletVelocity(self, t):
        if self.Ramp == True:
            vAvg = RampUpFunction(self.vAvg, t, self.timeRampUp)
            vMax = self.vDeviation + vAvg
            vMin = vAvg - self.vDeviation
        else:
            vMax = self.vMax
            vMin = self.vMin
            vAvg = self.vAvg
        dt = t - self.tprev

        # Track influence region
        # (Region1)
        if self.dir1 == 1:
            self.inf1 = self.inf1 + dt * self.v_sin
        elif self.dir1 == 0:
            self.inf1 = self.inf1 - dt * self.v_sin
        if self.inf1 < self.yMin:
            self.dir1 = 1
            self.inf1 = self.inf1 + 2 * (self.yMin - self.inf1)
        elif self.inf1 + self.Period > self.yMax:
            self.dir1 = 0
            self.inf1 = self.inf1 - 2 * ((self.inf1 + self.Period) - self.yMax)
        # (Region2)
        if self.dir2 == 1:
            self.inf2 = self.inf2 + dt * self.v_sin
        elif self.dir2 == 0:
            self.inf2 = self.inf2 - dt * self.v_sin
        if self.inf2 < self.yMin:
            self.dir2 = 1
            self.inf2 = self.inf2 + 2 * (self.yMin - self.inf2)
        elif self.inf2 + self.Period > self.yMax:
            self.dir2 = 0
            self.inf2 = self.inf2 - 2 * ((self.inf2 + self.Period) - self.yMax)
        for node in self.inletNodes:
            Y = node.Y
            # Unit Parabolic Base
            NormalParabola = (Y - self.yMax) * (
                Y - self.yMin) / (self.yMin * self.yMax) * vAvg
            # Add sine if inside area of influence
            if Y >= self.inf1 and Y <= self.inf1 + self.Period:
                S1 = sin((Y - self.inf1) / self.Period * 2 * pi)
            else:
                S1 = 0
            if Y >= self.inf2 and Y <= self.inf2 + self.Period:
                S2 = -sin((Y - self.inf2) / self.Period * 2 * pi)
            else:
                S2 = 0
            S = (S1 + S2) / 2 * (vMax - vAvg)
            U = NormalParabola + S
            node.SetSolutionStepValue(VELOCITY_X, 0, U)
            node.SetSolutionStepValue(VELOCITY_Y, 0, 0.0)
            self.tprev = t


class ConvolutedSinMultipliedOnParabolicInletVelocity:

    def __init__(self, inletNodes, vAverage, vDeviation, v_sin, X, RampUp, timeOfRampUp):
        self.inletNodes = inletNodes
        self.vMax = vAverage + vDeviation
        self.vMin = vAverage - vDeviation
        self.vAvg = vAverage
        self.v_sin = v_sin
        self.Period = X
        self.vDeviation = vDeviation

        yMin = 0
        yMax = yMin
        for node in self.inletNodes:
            if(node.Y > yMax):
                yMax = node.Y
            elif(node.Y < yMin):
                yMin = node.Y

        self.yMax = yMax
        self.yMin = yMin

        self.inf1 = yMin
        self.dir1 = 1
        self.inf2 = yMax - self.Period
        self.dir2 = 0
        self.tprev = 0

        self.Ramp = RampUp
        self.timeRampUp = timeOfRampUp

    def ApplyInletVelocity(self, t):
        if self.Ramp == True:
            vAvg = RampUpFunction(self.vAvg, t, self.timeRampUp)
            vMax = self.vDeviation + vAvg
            vMin = vAvg - self.vDeviation
        else:
            vMax = self.vMax
            vMin = self.vMin
            vAvg = self.vAvg
        dt = t - self.tprev

        # Track influence region
        # (Region1)
        if self.dir1 == 1:
            self.inf1 = self.inf1 + dt * self.v_sin
        elif self.dir1 == 0:
            self.inf1 = self.inf1 - dt * self.v_sin
        if self.inf1 < self.yMin:
            self.dir1 = 1
            self.inf1 = self.inf1 + 2 * (self.yMin - self.inf1)
        elif self.inf1 + self.Period > self.yMax:
            self.dir1 = 0
            self.inf1 = self.inf1 - 2 * ((self.inf1 + self.Period) - self.yMax)
        # (Region2)
        if self.dir2 == 1:
            self.inf2 = self.inf2 + dt * self.v_sin
        elif self.dir2 == 0:
            self.inf2 = self.inf2 - dt * self.v_sin
        if self.inf2 < self.yMin:
            self.dir2 = 1
            self.inf2 = self.inf2 + 2 * (self.yMin - self.inf2)
        elif self.inf2 + self.Period > self.yMax:
            self.dir2 = 0
            self.inf2 = self.inf2 - 2 * ((self.inf2 + self.Period) - self.yMax)
        for node in self.inletNodes:
            Y = node.Y
            # Unit Parabolic Base
            NormalParabola = (Y - self.yMax) * (
                Y - self.yMin) / (self.yMin * self.yMax) * vAvg
            # Add sine if inside area of influence
            if Y >= self.inf1 and Y <= self.inf1 + self.Period:
                S1 = sin((Y - self.inf1) / self.Period * 2 * pi)
            else:
                S1 = 0
            if Y >= self.inf2 and Y <= self.inf2 + self.Period:
                S2 = -sin((Y - self.inf2) / self.Period * 2 * pi)
            else:
                S2 = 0
            S = (S1 * S2) * (vMax - vAvg)
            U = NormalParabola + S
            node.SetSolutionStepValue(VELOCITY_X, 0, U)
            node.SetSolutionStepValue(VELOCITY_Y, 0, 0.0)
            self.tprev = t

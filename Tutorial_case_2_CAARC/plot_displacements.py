#!/usr/bin/python

import os
import sys
from math import *
from time import time
import matplotlib.pyplot as plt
import random
import numpy as np

disp_file = "output_csd/point_4.dat"

if not os.path.isfile(disp_file):
	print ("File not found at "+ disp_file)
	sys.exit()

file = open(disp_file,"r")
dis = file.readline()
time_list = []
displacement = []
while dis:
	splited_data = dis.split()
	if splited_data[0] == '#':
		pass
	else:
		time_list_data = splited_data[0]
		time_list.append(float(time_list_data))
		#displacement_data = splited_data[1]
		displacement_data = sqrt(pow(float(splited_data[1]),2) + pow(float(splited_data[2]),2) + pow(float(splited_data[3]),2))
		if(displacement_data != "-1.79769313486e+307"):
			displacement.append(float(displacement_data))
			old_displacement_data = displacement_data
		else:
			displacement.append(float(old_displacement_data))

	dis = file.readline()

#Plotting
plt.plot(time_list, displacement,color='b', linewidth=0.5)
#plt.xlim(0,25)

plt.xlabel('Time')
plt.ylabel('Displacement value')
plt.title("Displacement of a point in CAARC")
plt.suptitle("Kratos-OpenFOAM FSI simulation")
plt.show()



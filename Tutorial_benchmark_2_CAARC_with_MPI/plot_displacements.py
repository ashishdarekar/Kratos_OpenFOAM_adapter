#!/usr/bin/python

import os
import sys
from math import *
from time import time
import matplotlib.pyplot as plt
import random
import numpy as np

disp_file = "output_csd/point_1.dat"

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
		#print(len(splited_data))
		time_list_data = splited_data[0]
		time_list.append(float(time_list_data))
		#print(time_list)
		displacement_data = splited_data[3]
		if(displacement_data != "-1.79769313486e+307"):
			displacement.append(float(displacement_data))
			old_displacement_data = displacement_data
		else:
			displacement.append(float(old_displacement_data))
		#print(displacement)

	dis = file.readline()

#Plotting
plt.plot(time_list, displacement,color='b', linewidth=0.5)
#plt.ylim(-1.053e-2,-1.060e-2) # point 1
#plt.ylim(-1.78e-2,-1.66e-2) # point 2
#plt.ylim(-2.16e-2,-2.06e-2) # point 3
plt.xlim(0,25)

plt.xlabel('Time')
plt.ylabel('Displacement value')
plt.title("Displacement of a point in CAARC")
plt.suptitle("Kratos-OpenFOAM FSI simulation")
#x = [0 , 0]
#y = [3000, 3500]
#plt.plot(y,x,'--', color = 'k')
plt.show()



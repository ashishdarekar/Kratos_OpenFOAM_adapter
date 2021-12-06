#!/usr/bin/python

import os
import sys
from math import *
import matplotlib.pyplot as plt

disp_file = "postProcessing/probes/0/cellDisplacement"

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
		#print(splited_data)
		time_list_data = splited_data[0]
		time_list.append(float(time_list_data))

		displacement_data = splited_data[1]
		for word in displacement_data.split():
			print(word)
			if word.isdigit():
				#displacement.append(float(word))
				print(word)

		print(displacement)
		if(displacement_data != "-1.79769313486e+307"):
			#displacement.append(float(displacement_data))
			old_displacement_data = displacement_data
		else:
			displacement.append(float(old_displacement_data))

	dis = file.readline()

#print(time_list)
#print(displacement)


#Plotting
""" xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim() """

plt.plot(time_list, displacement,color='b', linewidth=0.5)
plt.xlabel('Time')
plt.ylabel('Displacement value')
plt.title(" Displacement of interface flap's tip ")
x = [0 , 0]
y = [0, 50]
plt.plot(y,x,'--', color = 'k')
plt.show()


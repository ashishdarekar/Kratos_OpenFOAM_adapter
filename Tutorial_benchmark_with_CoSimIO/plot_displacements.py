#!/usr/bin/python

import os
import sys
from math import *

disp_file = "postProcessing/probes/0/cellDisplacement"

if not os.path.isfile(disp_file):
	print "Forces file not found at "+ disp_file
	print "Be sure that the case has been run and you have the right directory!"
	print "Exiting."
	sys.exit()

file = open(disp_file,"r")
file_write = open("disp.txt","w")
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
		time_list.append(time_list_data)
		#print(time_list)
		displacement_data = splited_data[2]
		displacement.append(displacement_data)
		#print(displacement)
	dis = file.readline()

for i in range(len(time_list)):
	file_write.write(time_list[i])
	file_write.write(" ")
	file_write.write(displacement[i])
	file_write.write("\n")

os.system("./plot_disp_prep.sh")
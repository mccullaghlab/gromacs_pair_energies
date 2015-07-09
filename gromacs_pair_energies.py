#!/Users/martinmccullagh/anaconda/bin/python
##!/Library/Frameworks/Python.framework/Versions/Current/bin/python

#USAGE : gromacs_pair_energies.py [config file]
import scipy
import sys
import os
import numpy
import math
import MDAnalysis
from MDAnalysis.analysis.align import *

# constants


# Subroutines

# read the configuration file and populate the global variables
def ParseConfigFile(cfg_file):
	global top_file, traj_file, parm_file, gro_file, head_1, head_2, tail_1, tail_2, out_file
	f = open(cfg_file)
	for line in f:
		# first remove comments
		if '#' in line:
			line, comment = line.split('#',1)
		if '=' in line:
			option, value = line.split('=',1)
			option = option.strip()
			value = value.strip()
			# check value
			if option.lower()=='topfile':
				top_file = value
			elif option.lower()=='outfile':
				out_file = value
			elif option.lower()=='trajfile':
				traj_file = value
			elif option.lower()=='grofile':
				gro_file = value
			elif option.lower()=='parmfile':
				parm_file = value
			elif option.lower()=='head_1':
				head_1 = value
			elif option.lower()=='head_2':
				head_2 = value
			elif option.lower()=='tail_1':
				tail_1 = value
			elif option.lower()=='tail_2':
				tail_2 = value
			else :
				print "Option:", option, " is not recognized"

def computePbcDist(r1,r2,box):
	dist = 0

	for j in range(0,3):
		temp = r1[j]-r2[j]
		if temp < -box[j]/2.0:
			temp += box[j]
		elif temp > box[j]/2.0:
			temp -= box[j]
		dist += temp*temp

	dist = math.sqrt(dist)
	return dist;

def computePbcDist2(r1,r2,box):
	dist2 = 0

	for j in range(0,3):
		temp = r1[j]-r2[j]
		if temp < -box[j]/2.0:
			temp += box[j]
		elif temp > box[j]/2.0:
			temp -= box[j]
		dist2 += temp*temp

	return dist2;

# subroutine to parse top file
def parseTopFile(top_file):
	global sel_1, sel_2, sel_3, sel_4

	total_num_atoms = sel_1.numberOfAtoms() + sel_2.numberOfAtoms() + sel_3.numberOfAtoms() + sel_4.numberOfAtoms()

	atom_types = numpy.chararray(total_num_atoms,itemsize=8)
	charges = numpy.empty(total_num_atoms,dtype=float)

	with open(top_file, 'r') as txt:
	#line = txt.next().split()
		for lines in txt:
			line = lines.split()
			if len(line) >= 5:
				counter = 0
				for a in sel_1:
					if a.name == line[4]:
						atom_types[counter] = line[1]
						charges[counter] = float(line[6])*7.62421
						a.charge = charges[counter]
						a.type = atom_types[counter]
						a.mass = float(line[7])
					counter += 1
				for a in sel_2:
					if a.name == line[4]:
						atom_types[counter] = line[1]
						charges[counter] = float(line[6])*7.62421
						a.charge = charges[counter]
						a.type = atom_types[counter]
						a.mass = float(line[7])
					counter += 1
				for a in sel_3:
					if a.name == line[4]:
						atom_types[counter] = line[1]
						charges[counter] = float(line[6])*7.62421
						a.charge = charges[counter]
						a.type = atom_types[counter]
						a.mass = float(line[7])
					counter += 1
				for a in sel_4:
					if a.name == line[4]:
						atom_types[counter] = line[1]
						charges[counter] = float(line[6])*7.62421
						a.charge = charges[counter]
						a.type = atom_types[counter]
						a.mass = float(line[7])
					counter += 1

	txt.close
	return charges, atom_types

def parseParmFile(parm_file, atom_types):

	lj_params = numpy.empty((len(atom_types),2),dtype=float)
	with open(parm_file, 'r') as txt:
		for lines in txt:
			line = lines.split()
			if len(line) >= 1:
				for i in range(len(atom_types)):
					if line[0] == atom_types[i]:
						lj_params[i][0] = float(line[6])
						lj_params[i][1] = float(line[7])*10.0

	txt.close
	return lj_params


def computePairEnergies(charges, lj_params, box):
	global sel_1, sel_2, sel_3, sel_4

	pair_energies = numpy.zeros(6,dtype=float)

	counter1=0
	for a1 in sel_1:
		counter2 = sel_1.numberOfAtoms()
		for a2 in sel_2:
			dist2 = computePbcDist2(a1.pos, a2.pos, box)
			dist6 = dist2*dist2*dist2
			lj_eps = 4*math.sqrt(lj_params[counter1][0]*lj_params[counter2][0])
			lj_sigma2 = lj_params[counter1][1]*lj_params[counter2][1]
			lj_sigma6 = lj_sigma2*lj_sigma2*lj_sigma2
			pair_energies[0]+= lj_eps*(lj_sigma6*lj_sigma6/(dist6*dist6)-lj_sigma6/dist6)
			pair_energies[0]+= charges[counter1]*charges[counter2]/math.sqrt(dist2)
			counter2 += 1
		for a2 in sel_3:
			dist2 = computePbcDist2(a1.pos, a2.pos, box)
			dist6 = dist2*dist2*dist2
			lj_eps = 4*math.sqrt(lj_params[counter1][0]*lj_params[counter2][0])
			lj_sigma2 = lj_params[counter1][1]*lj_params[counter2][1]
			lj_sigma6 = lj_sigma2*lj_sigma2*lj_sigma2
			pair_energies[1]+= lj_eps*(lj_sigma6*lj_sigma6/(dist6*dist6)-lj_sigma6/dist6)
			pair_energies[1]+= charges[counter1]*charges[counter2]/math.sqrt(dist2)
			counter2 += 1
		for a2 in sel_4:
			dist2 = computePbcDist2(a1.pos, a2.pos, box)
			dist6 = dist2*dist2*dist2
			lj_eps = 4*math.sqrt(lj_params[counter1][0]*lj_params[counter2][0])
			lj_sigma2 = lj_params[counter1][1]*lj_params[counter2][1]
			lj_sigma6 = lj_sigma2*lj_sigma2*lj_sigma2
			pair_energies[2]+= lj_eps*(lj_sigma6*lj_sigma6/(dist6*dist6)-lj_sigma6/dist6)
			pair_energies[2]+= charges[counter1]*charges[counter2]/math.sqrt(dist2)
			counter2 += 1
		counter1 += 1
	for a1 in sel_2:
		counter2 = sel_1.numberOfAtoms() + sel_2.numberOfAtoms()
		for a2 in sel_3:
			dist2 = computePbcDist2(a1.pos, a2.pos, box)
			dist6 = dist2*dist2*dist2
			lj_eps = 4*math.sqrt(lj_params[counter1][0]*lj_params[counter2][0])
			lj_sigma2 = lj_params[counter1][1]*lj_params[counter2][1]
			lj_sigma6 = lj_sigma2*lj_sigma2*lj_sigma2
			pair_energies[3]+= lj_eps*(lj_sigma6*lj_sigma6/(dist6*dist6)-lj_sigma6/dist6)
			pair_energies[3]+= charges[counter1]*charges[counter2]/math.sqrt(dist2)
			counter2 += 1
		for a2 in sel_4:
			dist2 = computePbcDist2(a1.pos, a2.pos, box)
			dist6 = dist2*dist2*dist2
			lj_eps = 4*math.sqrt(lj_params[counter1][0]*lj_params[counter2][0])
			lj_sigma2 = lj_params[counter1][1]*lj_params[counter2][1]
			lj_sigma6 = lj_sigma2*lj_sigma2*lj_sigma2
			pair_energies[4]+= lj_eps*(lj_sigma6*lj_sigma6/(dist6*dist6)-lj_sigma6/dist6)
			pair_energies[4]+= charges[counter1]*charges[counter2]/math.sqrt(dist2)
			counter2 += 1
		counter1 += 1
	for a1 in sel_3:
		counter2 = sel_1.numberOfAtoms() + sel_2.numberOfAtoms() + sel_3.numberOfAtoms()
		for a2 in sel_4:
			dist2 = computePbcDist2(a1.pos, a2.pos, box)
			dist6 = dist2*dist2*dist2
			lj_eps = 4*math.sqrt(lj_params[counter1][0]*lj_params[counter2][0])
			lj_sigma2 = lj_params[counter1][1]*lj_params[counter2][1]
			lj_sigma6 = lj_sigma2*lj_sigma2*lj_sigma2
			pair_energies[5]+= lj_eps*(lj_sigma6*lj_sigma6/(dist6*dist6)-lj_sigma6/dist6)
			pair_energies[5]+= charges[counter1]*charges[counter2]/math.sqrt(dist2)
			counter2 += 1
		counter1 += 1

	return pair_energies

def centerOfMassWithWrap(sel,box):

	com = sel.masses()[0]*sel.coordinates()[0]
	total_mass = sel.masses()[0]
	for i in range(1,sel.numberOfAtoms()):
		for k in range(3):
			temp = sel.coordinates()[i][k]-sel.coordinates()[0][k]
			if temp < (-box[k]/2.0):
				temp += box[k]
			elif temp > (box[k]/2.0):
				temp -= box[k]
			com[k]+= sel.masses()[i]*temp
		total_mass += sel.masses()[i]

	com /= total_mass

	return com

# Main Program

# read in command line argument
cfg_file = sys.argv[1]

# read cfg file
ParseConfigFile(cfg_file)

print "Topology file:", top_file
print "Parameter file:", parm_file
print "Gro file:", gro_file
print "Trajectory file:", traj_file
print "Output file:", out_file

# initiate coordinate universe
coord = MDAnalysis.Universe(gro_file, traj_file)

# create atom selections
sel_1 = coord.selectAtoms(head_1)
sel_2 = coord.selectAtoms(head_2)
sel_3 = coord.selectAtoms(tail_1)
sel_4 = coord.selectAtoms(tail_2)
print "Number of atoms in head_1 selection:", sel_1.numberOfAtoms()
print "Number of atoms in head_2 selection:", sel_2.numberOfAtoms()
print "Number of atoms in tail_1 selection:", sel_3.numberOfAtoms()
print "Number of atoms in tail_2 selection:", sel_4.numberOfAtoms()
# obtain charges and atom types from topology file
charges,atom_types = parseTopFile(top_file)

#obtain lj parameters from parameter file for atom types
lj_params = parseParmFile(parm_file,atom_types)

# print some debug stuff
#for i in range(len(atom_types)):
#	print atom_types[i], charges[i], lj_params[i][0], lj_params[i][1]

#open output file
out = open(out_file,"w")
out.write("CV Distance    h1-h2      h1-t1      h1-t2      h2-t1     h2-t2      t1-t2\n")   

# Loop through trajectory
for ts in coord.trajectory:

	# get box
	box = coord.dimensions[0:3]

	# wrap and compute COM of sel_1 and sel_2
	sel_1_com = centerOfMassWithWrap(sel_1,box)
	sel_2_com = centerOfMassWithWrap(sel_2,box)

	# compute CV distance
	cv_dist = computePbcDist(sel_1_com,sel_2_com,box)

	# compute pair energies
	pair_energies = computePairEnergies(charges, lj_params, box)

	# print
	out.write("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n" % (cv_dist, pair_energies[0], pair_energies[1], pair_energies[2], pair_energies[3], pair_energies[4], pair_energies[5]))

out.close

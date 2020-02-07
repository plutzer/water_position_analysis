import Bio.PDB
from Bio.PDB.PDBParser import PDBParser
import numpy as np
from math import sqrt
from matplotlib import pyplot as plt

def get_structure(filename,structure_id = 'default'):
	parser = PDBParser(PERMISSIVE = 1)
	return parser.get_structure(structure_id,filename)

def center_of_mass(structure):
	weighted_coords = np.empty((0,3),'float32')
	total_mass = 0
	for model in structure:
		for chain in model:
			for residue in chain:
				for atom in residue:
					weighted_coords = np.append(weighted_coords,[atom.mass*atom.get_coord()],axis = 0)
					total_mass = total_mass + atom.mass
	center_of_mass = np.sum(weighted_coords,axis = 0)/total_mass
	return center_of_mass

def gyration_radius(structure):
	com = center_of_mass(structure)
	total_mass = 0
	weighted_elements_total = 0
	for model in structure:
		for chain in model:
			for residue in chain:
				for atom in residue:
					total_mass = total_mass + atom.mass
					atom_coords = atom.get_coord()
					dist_from_center = sqrt((com[0]-atom_coords[0])**2 + (com[1]-atom_coords[1])**2 + (com[2]-atom_coords[2])**2)
					weighted_elements_total = weighted_elements_total + atom.mass*(dist_from_center**2)
	return sqrt(weighted_elements_total/total_mass)

def water_distances(structure,coords):
	distances = []
	for chain in structure[0]:	## Uses the first model in the file as default to avoid using water molecules twice.
		for residue in chain:
			if residue.id[0] == 'W':
				for atom in residue:
					atc = atom.get_coord()
					distance = sqrt((coords[0]-atc[0])**2 + (coords[1]-atc[1])**2 + (coords[2]-atc[2])**2)
					distances.append(distance)
	return distances

def water_distance_plot(distances):
	distances.sort()
	plt.scatter(range(len(distances)),distances)
	plt.show()
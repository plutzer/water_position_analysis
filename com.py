import Bio.PDB
from Bio.PDB.PDBParser import PDBParser
import numpy as np
from math import sqrt
from matplotlib import pyplot as plt


def get_com(filename,structure_id = 'default'):
	parser = PDBParser(PERMISSIVE = 1)
	structure = parser.get_structure(structure_id,filename)
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

def water_distance_plot(filename):
	parser = PDBParser(PERMISSIVE = 1)
	structure = parser.get_structure('default',filename) ## Will need to eventually change this so we aren't parsing twice.
	com = get_com(filename)
	distances = []
	watercount = 0
	watercount_2 = 0
	for chain in structure[0]:	## Uses the first model in the file as default to avoid using water molecules twice.
		for residue in chain:
			if residue.id[0] == 'W':
				watercount_2 = watercount_2 + 1
				for atom in residue:
					watercount = watercount + 1
					atc = atom.get_coord()
					water_distance = sqrt((com[0]-atc[0])**2 + (com[1]-atc[1])**2 + (com[2]-atc[2])**2)
					distances.append(water_distance)
	distances.sort()
	plt.scatter(range(len(distances)),distances)
	plt.show()
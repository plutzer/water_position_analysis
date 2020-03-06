from water_functions import *
import sys
import seaborn as sns
import os

def run_grandpa(pdb_name):
	structure_file = download_structure(pdb_name)
	structure = get_structure(structure_file)
	com = center_of_mass(structure)
	rad = gyration_radius(structure)
	distances = water_distances(structure,com)
	distances.sort()
	gyr_norm_dist = distances/rad

	#Produce the plots
	sns.set(style="white", palette="muted", color_codes=True)
	f,axes = plt.subplots(2, 2, figsize=(10, 8), sharex=False) #Change this to resize to grandpas screen
	#Top left - normalized water distances
	axes[0,0].scatter(range(len(gyr_norm_dist)),gyr_norm_dist)
	axes[0,0].title.set_text('Water Distances by Water Index')
	axes[0,0].set_ylabel('Distance from COM (Gyration Radii)')
	axes[0,0].set_xlabel('Water Index Number')
	#Top right - Kernel Density of normalized distances
	sns.distplot(gyr_norm_dist, hist=False, rug=True, color="r", ax=axes[0, 1])
	axes[0,1].title.set_text('Distribution of Water Distances')
	axes[0,1].set_xlabel('Distance from COM (Gyration Radii)')
	axes[0,1].set_ylabel('Density')
	plt.savefig(os.getcwd() + '/Waterplots/' + pdb_name)
	plt.show()

pdb_name = sys.argv[1]
#print(pdb_name)
run_grandpa(pdb_name)
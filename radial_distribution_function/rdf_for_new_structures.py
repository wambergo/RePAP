import numpy as np
import matplotlib.pyplot as plt
import ase
import asap3
from asap3.analysis.rdf import RadialDistributionFunction
from ase.visualize import view
from ase import Atoms


# Import our generated structure
generated_structure = np.load("structures_array_m4_s025_unsort.npy") # dim(n_snaps, n_atoms, 3)
n_generations= len(generated_structure) # obtain number of generations

# renormalize
n_grid_points = 50
grid = np.linspace(0, 10, n_grid_points + 1)
d = grid[:-1]+(grid[1]-grid[0])/2 # taking mid value of linspace

for i in range(n_generations):
    for j in range(len(generated_structure[0])):
        for k in range(3):
            generated_structure[i][j][k] = d[int(generated_structure[i][j][k])]

# Random placement test
random_structure = 10 * np.random.rand(n_generations, len(generated_structure[0]), 3)


# initialize trajectory object
new_structure = []

for i in range(n_generations):
    new_structure.append(Atoms(symbols='Si30', pbc=True, cell=[10.0, 10.0, 10.0])) # initialize objects
    new_structure[i].set_positions(generated_structure[i])

# Determine the RDF of new structure
# first setup RDF out to a distance of X Rng, with Y bins
rng= 10.0
bins = 100

RDFobj = None
for atoms in new_structure:
    if RDFobj is None:
        RDFobj = RadialDistributionFunction(atoms, rng, bins)
    else:
        RDFobj.atoms = atoms  # Fool RDFobj to use the new atoms
    RDFobj.update()           # Collect data
    
rdf_new = RDFobj.get_rdf()

# Get the RDF and plot it.
rdf = RDFobj.get_rdf()
x = np.arange(bins) * rng / bins
plt.plot(x, rdf)
plt.show()

# view structs
view(new_structure[0])
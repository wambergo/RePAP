import numpy as np
import matplotlib.pyplot as plt
import ase
import asap3
from asap3.analysis.rdf import RadialDistributionFunction
from ase.visualize import view
from ase import Atoms


# Import our generated structure
generated_structure = np.load("structures_array.npy")
n_generations= len(generated_structure) # obtain number of generations

# renormalize
generated_structure /= 5

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
    
rdf = RDFobj.get_rdf()

# Get the RDF and plot it.
rdf = RDFobj.get_rdf()
x = np.arange(bins) * rng / bins
plt.plot(x, rdf)
plt.show()
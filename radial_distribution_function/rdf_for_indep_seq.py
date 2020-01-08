import numpy as np
import matplotlib.pyplot as plt
import ase
import asap3
from asap3.analysis.rdf import RadialDistributionFunction
from ase.visualize import view
from ase import Atom

from ase.io import Trajectory
from ase.visualize import view # run view(atom_objects[0])


import os
trajs = [f for f in os.listdir('.') if f.endswith('.traj')]

def load_trajectory(filename):
    """load trajectory
    :param filename:
    :return list of ASE atom objects
    """
    atom_objects = Trajectory(filename)

    return atom_objects

# Load atoms objects from ASE trajectory file
    

n_snaps = len(trajs) # number of used snapshots 

traj = []
trunc_traj = []

for i in range(n_snaps):
    atoms = load_trajectory(trajs[i])
    trunc_traj.append(atoms[0])
    atoms.close()


# Setup RDF out to a distance of X Rng, with Y bins
rng= 6.0
bins = 100

# Calculating the RDF of a previously stored trajectory
# https://wiki.fysik.dtu.dk/asap/Radial%20Distribution%20Functions
RDFobj = None
for atoms in trunc_traj:
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


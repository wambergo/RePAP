import numpy as np
import matplotlib.pyplot as plt
import ase
import asap3
from asap3.analysis.rdf import RadialDistributionFunction
from ase.visualize import view
from ase import Atom

from ase.io import Trajectory

def load_trajectory(filename):
    """load trajectory

    :param filename:
    :return list of ASE atom objects
    """
    atom_objects = Trajectory(filename)
    return atom_objects

# Load atoms objects from ASE trajectory file
traj = "new.traj"
atoms = load_trajectory(traj)
snaps = len(atoms) # 175983

# truncation of trajectory
n_snaps = 8000 # number of used snapshots 
trunc_traj = []
for i in range(n_snaps):
    trunc_traj.append(atoms[i])

# Setup RDF out to a distance of X Rng, with Y bins
rng= 10.0
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

view(trunc_traj[0])
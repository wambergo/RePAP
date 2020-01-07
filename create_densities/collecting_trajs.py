import ase
import numpy as np
import torch
import argparse

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
    #atom_objects.close()
    return atom_objects

def atoms_to_array(atom_objects, n_train):
    """Take an ASE atom object and convert it into a matrix of size
    (n_atoms, 3) such that each row is as atom and the columns are the xyz
    coordinates
    :param atoms: Atoms objects defining atoms position in space
    :type atoms: ase.Atoms
    :return numpy array of dimension (n_atoms, 3)
    """
    pos = []

    for i in range(n_train):
        atoms = atom_objects[i]
        xyz = atoms.get_positions() # get xyz coordinates
        pos.append(xyz)

    return pos

# Create matrix of size (n_atoms, 3) with the coordinates of each atom
n_train = 1
pos = []

for i in range(len(trajs)):
    filename = trajs[i]
    atom_objects = load_trajectory(filename)
    
    pos.append(atoms_to_array(atom_objects, n_train))

# array formating 
ny_pos=[]

for i in range(len(pos)):
    ny_pos.append(pos[i][0])


np.save('ny_indep_trajs.npy', np.asarray(ny_pos))

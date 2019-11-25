import ase
import numpy as np
import torch
import argparse

from ase.io import Trajectory

def load_trajectory(filename):
    """load trajectory

    :param filename:
    :return list of ASE atom objects
    """
    atom_objects = Trajectory(filename)
    return atom_objects

def atoms_to_array(atom_objects, n_train):
    """atoms_to_array
    Take an ASE atom object and convert it into a matrix of size
    (n_train, n_atoms, 3) such that each row is as atom and the columns
    are the xyz coordinates

    :param atoms: Atoms objects defining atoms position in space
    :type atoms: ase.Atoms

    :return numpy array of dimension (n_train, n_atoms, 3)
    """
    pos = []

    for i in range(n_train):
        atoms = atom_objects[i]
        xyz = atoms.get_positions() # get xyz coordinates
        pos.append(xyz)

    return pos


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('trajectories', help='Location of trajectory file')
    parser.add_argument('--n_train', help='Number of training examples', type=int, default=1)
    parser.add_argument('save_name', help="Name of the saved position array")
    args = parser.parse_args()

    # Load atoms objects from ASE trajectory file
    filename = args.trajectories
    atom_objects = load_trajectory(filename)

    # Create matrix of size (n_atoms, 3) with the coordinates of each atom
    pos = atoms_to_array(atom_objects, args.n_train)

    # Save position array
    np.save(args.save_name, pos)

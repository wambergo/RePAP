import ase
import numpy as np
import torch

def atoms_to_array(atoms: ase.Atoms):
    """atoms_to_array
    Take an ASE atoms object and convert it into a matrix of size
    (n_atoms, 3) such that each row is an atom and the columns are the xyz coordinates

    :param atoms: Atoms object defining atoms positions in space
    :type atoms: ase.Atoms

    :return: numpy array of dimension (n_atoms, 3)
    """

    # Sort atoms according to distance from origin

    # Create matrix of size (n_atoms, 3) with the coordinate of each atom

    return putting_sequence

def gaussian_expansion(input_array, n_grid_points, d_max):
    """gaussian_expansion
    Expand the last dimension of input_array such that an
    extra dimension of size n_grid_points is added.
    The extra dimension is a smoothed one-hot encoding of the scalar values of the last dimension.

    :param input_array:
    :return: expanded_array
    """

    return expanded_array

def load_trajectory(filename):
    """load_trajectory

    :param filename:
    :return: list of ASE atom objects
    """
    pass


def main():
    filename = "ase.traj"
    # Load atoms objects from ASE trajectory file
    atom_objects = load_trajectory(filename)

    # Convert to matrices suitable for training
    dataset = [atoms_to_array(at) for at in atom_objects]

if __name__ == "__main__":
    main()

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

def mean_sorting(pos, n_train):
    """mean_sorting
    Take position array from the atoms_to_array function, and sort atoms
    by their mean position to the origin

    :param pos: (n_train, n_atoms, 3) array obtained from atoms_to_array
    :param n_train: number of snapshots used for sorting
    :return numpy array of dimension (n_train, n_atoms, 3)
    """

    mean_xyz = np.mean(pos, axis=0)
    mean_dist = np.linalg.norm(mean_xyz, axis=1)

    pos_mean_sort = []
    putting_sequence = []

    for i in range(n_train):
        pos_snap_i = pos[i]
        pos_mean_sort = pos_snap_i[mean_dist.argsort()]
        putting_sequence.append(pos_mean_sort)

    return putting_sequence

def gaussian_expansion(input_array, n_grid_points, sigma, n_train, n_type):
    """gaussian expansion
    Expand the last dimension of input array such that an extra dimension
    of size n_grid is added.
    The extra dimension is a smoothed one-hot encoding of the scalar values
    of the last dimension.

    :param input_array:
    :n_grid_points: Number of grid points, i.e. resolution of the gaussain expansion
    :sigma: Kernel parameter for density smoothing
    :n_train: size of the expanded array
    :n_type: number of atom types, default = 1
    :return expanded_array
    """
    x_grid = np.linspace(0,10,n_grid_points + 1)

    pos_matrix = input_array[:n_train] # the sequence with n training examples
    d = x_grid[:-1]+(x_grid[1]-x_grid[0])/2 # taking mid value of linspace
    d = np.tile(d,(n_type,1)) # concatenate n_type column
    d = np.transpose(d) # make it dim(n_grid_points,n_type)

    # obtain n_atoms from first snapshot
    n_atoms = len(pos_matrix[0])

    # initialize density matrix
    density_matrix = np.zeros((n_train, n_atoms), object)

    for i in range(n_train):
        for j in range(n_atoms):
            density_matrix[i][j] = np.exp(-(np.abs(d-pos_matrix[i][j])**2) / (2*sigma**2))

    expanded_array = density_matrix.tolist() # list converting
    return expanded_array


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('trajectories', help='Location of trajectory file')
    parser.add_argument('--n_grid_points', help='Number of grid points', type=int, default=50)
    parser.add_argument('--sigma', help='Kernel parameter for density smoothing [1 Å]',
                        type=int, default=0.1)
    parser.add_argument('--n_train', help='Number of training examples', type=int, default=1)
    parser.add_argument('--n_type', help='Types of atoms', type=int, default=1)
    args = parser.parse_args()

    # Load atoms objects from ASE trajectory file
    filename = args.trajectories
    atom_objects = load_trajectory(filename)

    # Create matrix of size (n_atoms, 3) with the coordinates of each atom
    pos = atoms_to_array(atom_objects, args.n_train)

    # sort atoms by atoms mean distance to the origin
    putting_sequence = mean_sorting(pos, args.n_train)

    # Make the gaussian expansion
    expanded_array = gaussian_expansion(putting_sequence, args.n_grid_points,
                                        args.sigma, args.n_train, args.n_type)

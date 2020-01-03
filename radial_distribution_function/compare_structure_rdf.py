import numpy as np
import matplotlib.pyplot as plt
import ase
import asap3
import argparse
from asap3.analysis.rdf import RadialDistributionFunction
from ase.visualize import view
from ase import Atoms
from ase.io import Trajectory

def load_trajectory(filename):
    """load trajectory

    :param filename:
    :return list of ASE atom objects
    """
    atom_objects = Trajectory(filename)
    return atom_objects

def renormalize_gen(generated_structure):
        """renormalize_gen
        Take the coordinates from the generated structure and place them within
        0 to 10 Å.
        :param structure: the generated structre
        :return: normalized structure list
        """
        n_generations= len(generated_structure) # obtain number of generations

        # renormalize
        n_grid_points = 50
        grid = np.linspace(0, 10, n_grid_points + 1)
        d = grid[:-1]+(grid[1]-grid[0])/2 # taking mid value of linspace
        for i in range(n_generations):
            for j in range(len(generated_structure[0])):
                for k in range(3):
                    generated_structure[i][j][k] = d[int(generated_structure[i][j][k])] + min(0.1, max(-0.1, np.random.uniform(-0.1,0.1))) # added noise for smooting
        return generated_structure

def rdf(atoms, rng, bins):
    """rdf
    Calculate radial distribution functions:
    measures the probability of finding an atom at distance r given that there
    is an atom at position 0; it is thus essentially a histogram of interatomic
    distances - and is calculated as such

    :param structure: the atoms object being analyzed
    :param rng: max distance up to which RDF is calculated, in Ångström
    :param bins: number of bins used for histogram, i.e. resolution of the RDF
    :return n_generations: the number of generations
    :return rdf vector: values of the RDF
    """
    # number of generations
    n_generations = len(atoms)
    # initialize trajectory object
    new_structure = []

    for i in range(n_generations):
        new_structure.append(Atoms(symbols='Si30', pbc=True, cell=[10.0, 10.0, 10.0])) # initialize objects
        new_structure[i].set_positions(atoms[i])
    # Determine RDF
    RDFobj = None
    for atoms in new_structure:
        if RDFobj is None:
            RDFobj = RadialDistributionFunction(atoms, rng, bins)
        else:
            RDFobj.atoms = atoms  # Fool RDFobj to use the new atoms
        RDFobj.update()           # Collect data
    # obtain RDF
    rdf = RDFobj.get_rdf()
    return n_generations, rdf

def rdf_org_traj(atoms, rng, bins, used_snaps):
    """rdf
    Calculate radial distribution functions:
    measures the probability of finding an atom at distance r given that there
    is an atom at position 0; it is thus essentially a histogram of interatomic
    distances - and is calculated as such

    :param structure: the ASE trajectory file
    :param rng: max distance up to which RDF is calculated, in Ångström
    :param bins: number of bins used for histogram, i.e. resolution of the RDF
    :param used_snaps: the number of used snapshots from the trajectory file
    :return rdf vector: values of the RDF
    """
    # truncation of trajectory
    trunc_traj = []
    for i in range(len(atoms)):
        trunc_traj.append(atoms[i])

    # Determine RDF
    RDFobj = None
    for atoms in trunc_traj:
        if RDFobj is None:
            RDFobj = RadialDistributionFunction(atoms, rng, bins)
        else:
            RDFobj.atoms = atoms  # Fool RDFobj to use the new atoms
        RDFobj.update()           # Collect data
    # obtain RDF
    rdf = RDFobj.get_rdf()
    return rdf

def rand_struct(side_len, n_generations):
    """rand_struct
    Create a random strucutre of same size as original and created structure

    :param axis_len: the side length of the box, default: 10 Ångstrøm
    :param n_generations: number of generations created
    """
    n_atoms = 30 # number of atoms in the box
    rand_struct = side_len * np.random.rand(n_generations, n_atoms, 3)
    return rand_struct

def comp_plot(struct_1, struct_2, struct_3, rng, bins):
    """comp_plot
    Compare RDF for structures by plotting

    :param structure_1,2,3: generated structure, original structure, random structure
    """
    x = (np.arange(bins) + 0.5) * rng / bins  # compute distances
    plt.plot(x, struct_1, label="Generated structure")
    plt.plot(x, struct_2, '--', label="Original structure")
    plt.plot(x, struct_3, label="Random structure")
    plt.xlabel("Distances [Å]")
    plt.ylabel("RDF")
    plt.legend()
    plt.grid(True)
    plt.show()
    return None

# Measure functions
def rmse(gen_struct_rdf, org_struct_rdf):
    """rmse
    root mean square error function for RDF
    """
    return np.sqrt(((gen_struct_rdf - org_struct_rdf) ** 2).mean())

def mse(gen_struct_rdf, org_struct_rdf):
    """mse
    mean square error function for RDF
    """
    return (np.square(gen_struct_rdf - org_struct_rdf)).mean()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('gen_struct', help='Location of generated structure file')
    parser.add_argument('org_struct', help='Location of original structure file')
    parser.add_argument('--rng', help='Maximum distance up to which RDF is calcuated in Ångstrøm', type=float, default=10.0)
    parser.add_argument('--bins', help='Resolution of the RDF', type=int, default=100)
    args = parser.parse_args()

    # Load our structure file and analyze it
    filename = args.gen_struct
    generated_structure = np.load(filename)
    generated_structure = renormalize_gen(generated_structure)
    n_generations, gen_struct_rdf = rdf(generated_structure, args.rng, args.bins)

    # Load the original trajectory file and analyze it
    original_structure = load_trajectory(args.org_struct)
    org_struct_rdf = rdf_org_traj(original_structure, args.rng, args.bins, n_generations)

    # Create a totally random structure for comparison
    random_structure = rand_struct(args.rng, n_generations)
    _, rand_struct_rdf = rdf(random_structure, args.rng, args.bins)

    # Compare RDF for structures
    comp_plot(gen_struct_rdf, org_struct_rdf, rand_struct_rdf, args.rng, args.bins)
    error_rdf = mse(gen_struct_rdf, org_struct_rdf)
    print("RDF reconstruction error:", error_rdf)
    
    


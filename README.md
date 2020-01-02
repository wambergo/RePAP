# RePAP <img src="https://github.com/wambergo/Atomic-structure-generation-with-recurrent-neural-networks/blob/master/other/math/various%20figures/repap_logo.jpg" align="right" />

## Overview
RePAP provides a atom structure generation model using recurrent probabilistic axis projections.

## Prerequisites
- Python 3.5+
- NumPy 1.10+: ```pip3 install numpy```
- SciPy 0.17+: ```pip3 install scipy```
- PyTorch: refer to [PyTorch website](https://pytorch.org/get-started/locally/)
- ASE: refer to [ASE website](https://wiki.fysik.dtu.dk/ase/install.html)

## Create densities
The probabilistic axis projection framework consists of representing each atom (in each snapshot
of the molecular training trajectory) as a discretized Gaussian density cloud that is projected onto
each of the three axes x, y, and z.

With 30 atoms in total, this corresponds to 90 axis projections (i.e. 90 atom coordinates).
The 30 atoms in a snapshot are ordered w.r.t. mean distance to origin over the entire trajectory, from smallest to largest distance. The 90 atom coordinates are then ordered as
x0, y0, z0, x1, y1, z1, ..., x29, y29, z29. This is the sequential order fed to the recurrent network.

With a 50-dimensional grid along each axis, the density of the i’th atom coordinate at the k’th
element of the grid vector is

![alt text](https://github.com/wambergo/Atomic-structure-generation-with-recurrent-neural-networks/blob/master/other/math/gauss_exp.png)

### Usage
(applies to all create_densites_xxx.py files)

In order to take an ASE atom object and make the gaussian expansion run the following in the terminal
```
python3 create_densities_xxx.py "/location_of_ASE_traj_file.traj"
```
Optional arguments:
- The number of grid points ```n_grid_points```, default=50
- Kernel parameter for density smoothing ```sigma```, default=0.1
- Number of training examples ```n_train```, default=1
- Types of atoms ```n_type```, default=1

## Model performance (radial_distribution_function)
As a measure of the quality of the generated samples, we consider the radial distribution function (RDF) which measures the average number density of atoms at a distance r.

We define the RDF reconstruction error as

![alt text](https://github.com/wambergo/Atomic-structure-generation-with-recurrent-neural-networks/blob/master/other/math/rdf_error.png) 

and is used to guide model development.

### Compare radial distribution functions
The compare_structure_rdf.py script can be used to compare the obtained radial distibution functions from our generation. You run the following in the terminal
```
python3 compare_structure_rdf.py "/location_of_generated_structure_converted_to_numpy.npy" "/location_of_original_ASE_traj_file.traj"
```
Optional arguments:
- Maximum distance up to which RDF is calcuated in Ångstrøm ```rng```, default=10
- Kernel parameter for density smoothing ```sigma```, default=0.1
- Resolution of the RDF ```bins```, default=100

### Example of comparison
Here is a comparison of the radial distribution function of one of our generated structures and the original structure 
![Screenshot](radial_distribution_function/plots/compare_rdfs.png)


## Google Colab

### Main setup
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1n1ZYoIAycpqcqS1h6HUPTPRlss9fMhUP#scrollTo=tCwAinf6qZDa)

### Testing setup
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1lCU9xOw1He0ArtYbdtE179HwEsGbrsul)

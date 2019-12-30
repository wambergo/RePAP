# Atomic-structure-generation-with-recurrent-neural-networks

## Create densities
The probabilistic axis projection framework consists of representing each atom (in each snapshot
of the molecular training trajectory) as a discretized Gaussian density cloud that is projected onto
each of the three axes x, y, and z.

With 30 atoms in total, this corresponds to 90 axis projections (i.e. 90 atom coordinates).
The 30 atoms in a snapshot are ordered w.r.t. mean distance to origin over the entire trajectory, from smallest to largest distance. The 90 atom coordinates are then ordered as
x0, y0, z0, x1, y1, z1, ..., x29, y29, z29. This is the sequential order fed to the recurrent network.

With a 50-dimensional grid along each axis, the density of the i’th atom coordinate at the k’th
element of the grid vector is



## First setup
https://colab.research.google.com/drive/1n1ZYoIAycpqcqS1h6HUPTPRlss9fMhUP#scrollTo=tCwAinf6qZDa

## Second setup
https://colab.research.google.com/drive/1lCU9xOw1He0ArtYbdtE179HwEsGbrsul

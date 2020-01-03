# RDF Plot examples <img src="https://github.com/wambergo/RePAP/blob/master/radial_distribution_function/rdf_figure.png" align="right" />

- Evaluate radial distribution function for original structure (.traj files): rdf_for_org_structure.py
- Evaluate radial distribution function for generated structure: rdf_for_new_structures.py
- Compare structures with: compare_structure_rdf.py

## RePAP Vanilla with 40000 snapshots used, sigma = 0.10 and layers with 256 units. Dim = 50
RDF reconstruction error: 351.32

![Screenshot](plots/compare_vanilla.png)

## RePAP+ Atom # with 40000 snapshots used, sigma = 0.10 and layers with 256 units. Dim = 139
RDF reconstruction error: 0.12050458005585839
![Screenshot](plots/compare_s01_40000_l1_256hidden.png)

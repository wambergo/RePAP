def mean_sorting(pos, n_train):
    """mean_sorting
    Take position array from the atoms_to_array function, and sort atoms
    by their mean position to the origin

    :param pos: (n_train, n_atoms, 3) array obtained from atoms_to_array
    :param n_train: number of snapshots used for sorting
    :return numpy array of dimension (n_train, n_atoms, 3)
    """

    mean_xyz = np.mean(pos[0:n_train], axis=0)
    mean_dist = np.linalg.norm(mean_xyz, axis=1)

    pos_mean_sort = []
    putting_sequence = []

    for i in range(len(pos)):
        pos_snap_i = pos[i]
        pos_mean_sort = pos_snap_i[mean_dist.argsort()]
        putting_sequence.append(pos_mean_sort)

    return putting_sequence

def original_sorting(pos):    
    putting_sequence = []
    
    for i in range(len(pos)):
        xyz = pos[i]
        dist = np.linalg.norm(xyz, axis = 1) # compute distance to origin
        pos_sorted = xyz[dist.argsort()] # sort the atoms by their distances
        putting_sequence.append(pos_sorted)
        
    return putting_sequence

def periodic_sorting(pos, unit_cell_len):
    putting_sequence = []
    corners8 = np.array([[0,0,0],
                        [unit_cell_len,0,0],
                        [0,unit_cell_len,0],
                        [0,0,unit_cell_len],
                        [unit_cell_len,unit_cell_len,0],
                        [unit_cell_len,0,unit_cell_len],
                        [0,unit_cell_len,unit_cell_len],
                        [unit_cell_len,unit_cell_len,unit_cell_len]])
    
    for i in range(len(pos)):
        xyz = pos[i]

        dist_vec=np.zeros((n_atoms,8))
        for j in range(8):
          dist_vec[:,j] = np.linalg.norm(xyz-corners8[j], axis = 1) # compute distance to origin
        
        dist = np.min(dist_vec,1)
        pos_sorted = xyz[dist.argsort()] # sort the atoms by their distances
        putting_sequence.append(pos_sorted)
        
    return putting_sequence

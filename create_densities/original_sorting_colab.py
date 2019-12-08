# new function for original sorting
def original_sorting(pos):

    putting_sequence = []

    for i in range(len(pos)):
        xyz = pos[i]
        dist = np.linalg.norm(xyz, axis = 1) # compute distance to origin
        pos_sorted = xyz[dist.argsort()] # sort the atoms by their distances
        print(4)
        putting_sequence.append(pos_sorted)

    return putting_sequence

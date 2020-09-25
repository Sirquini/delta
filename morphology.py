import numpy as np
from scipy import ndimage
from itertools import product, combinations

def complement(a, b):
    if a == 1 and b == 0:
        return 1
    return 0

vcomplement = np.vectorize(complement)

class MatrixLattice():
    def __init__(self, shape):
        """Create a lattice for a matrix with dimensions
        set by the pair of ints `shape`.
        """
        self.shape = shape
        self.bottom = np.zeros(shape, dtype=int)
        self.top = np.ones(shape, dtype=int)

    def __len__(self):
        """Returns the number of nodes in the lattice.
        """
        return self.shape[0] * self.shape[1]

    def entails(self, a, b):
        return np.all(np.greater_equal(a, b))

    def lub(self, iterable):
        """Least Upper Bound of `iterable`.
        """
        r = np.zeros(self.shape, dtype=int)
        for i in iterable:
            r = self.pair_lub((r, i))
        return r
    
    def pair_lub(self, pair):
        return np.logical_or(pair[0], pair[1]).astype(int)

    def glb(self, iterable):
        """Greatest Lower Bound of `iterable`.
        """
        r = np.ones(self.shape, dtype=int)
        for i in iterable:
            r = self.pair_glb((r,i))
        return r

    def pair_glb(self, pair):
        return np.multiply(pair[0], pair[1])

    def imply(self, a, b):
        """Returns a imply b.
        """
        return vcomplement(b, a)
    
    def join_irreducible_elements(self):
        """Returns a list of all the join-irreducible elements in the lattice.

        `x` is join-irreducible (or lub-irreducible) if:
        + `x != 0`
        + `a < x` and `b < x` imply `a lub b < x` for all `a` and `b` 
        """
        pairs = product(range(self.shape[0]), range(self.shape[1]))
        result = []
        for a, b in pairs:
            tmp = self.bottom.copy()
            tmp[a, b] = 1
            result.append(tmp)
        return result

def downset(element):
    mask = element.flatten()
    sequence = [0 for i in mask if i == 1]
    top = [1 for _ in sequence]
    result = []
    result.append(np.zeros(element.shape, dtype=int))
    while sequence != top:
        for pos, value in enumerate(sequence):
            if value + 1 <= 1:
                sequence[pos] = value + 1
                break
            else:
                sequence[pos] = 0
        tmp = np.zeros(len(mask), dtype=int)
        count = 0
        for pos, value in enumerate(mask):
            if value == 1:
                tmp[pos] = sequence[count]
                count += 1
        result.append(tmp.reshape(element.shape))
    return result


# ######################################
# Delta functions defined with the class
# MatrixLattice and local context
# ######################################

def partition_helper(lattice, functions, first, last, c, helper_cache, search_space):
    cached_result = helper_cache.get((tuple(map(tuple, c)), first, last-1))
    if cached_result is not None:
        return cached_result
    # No need to compute for the bottom
    if np.array_equal(lattice.bottom, c):
        return lattice.bottom
    fn_num = last - first
    if fn_num == 1:
        return ndimage.binary_dilation(c, structure=functions[first]).astype(c.dtype)
    else:
        mid_point = first + fn_num // 2
        result = lattice.glb(lattice.lub((partition_helper(lattice, functions, first, mid_point, a, helper_cache, search_space), partition_helper(lattice, functions, mid_point, last, lattice.imply(a, c), helper_cache, search_space))) for a in search_space if lattice.entails(c, a) )
        helper_cache[(tuple(map(tuple, c)), first, last-1)] = result
        return result

def delta5_plus(lattice, functions, image):
    """Calculate Delta* for a set of `dilation functions` over a `lattice` for an `image`
    partitioning the set of functions and using a look-up table.
    """
    n = len(functions)
    helper_cache = {}
    # Precalculate lattice search space
    search_space = downset(image)
    return partition_helper(lattice, functions, 0, n, image, helper_cache, search_space)

def apply_dilation(structure, image):
    """Wrapper for `ndimage.binary_dilation`."""
    return ndimage.binary_dilation(image, structure=structure).astype(image.dtype)

def delta_plus_jies(lattice, functions, image):
    """Calculate Delta+ for a set of dilation `functions` over a `lattice` for an `image`.
    
    Similar to `delta_partition`.

    This implementation takes advantage of the join irreducible
    elements to reduce the number of operations.

    Args:
      lattice: A Lattice instance.
      functions: A list of space functions.
    """
    if len(functions) == 1:
        return apply_dilation(functions[0], image)
    jie_s = lattice.join_irreducible_elements()
    return lattice.lub(lattice.glb(apply_dilation(fn, ji) for fn in functions) for ji in jie_s if lattice.entails(image, ji))

def run_dilations(entry, struct1, struct2):
    shape = entry.shape
    lattice = MatrixLattice(shape)
    dilation1 = ndimage.binary_dilation(entry, structure=struct1).astype(entry.dtype)
    dilation2 = ndimage.binary_dilation(entry, structure=struct2).astype(entry.dtype)
    dilation3 = ndimage.binary_dilation(entry, structure=np.logical_and(struct1, struct2)).astype(entry.dtype)
    print("\nDilation results:")
    print(dilation1)
    print(dilation2)
    print("Intersection d1, d2:")
    print(lattice.glb((dilation1, dilation2)))
    print("Dilation of intersecting s1, s2:")
    print(dilation3)
    print("Delta result:")
    print(delta_plus_jies(lattice, (struct1, struct2), entry))

def run_dilations2(entry, struct1, struct2):
    import matplotlib.pyplot as plt
    from time import perf_counter

    shape = entry.shape
    l = MatrixLattice(shape)
    dilation1 = ndimage.binary_dilation(entry, structure=struct1).astype(entry.dtype)
    dilation2 = ndimage.binary_dilation(entry, structure=struct2).astype(entry.dtype)
    dilation3 = ndimage.binary_dilation(entry, structure=np.logical_and(struct1, struct2)).astype(entry.dtype)
    print("\nDilation results:")
    print(dilation1)
    print(dilation2)
    print("Dilation of intersecting s1, s2:")
    print(dilation3)
    print("Old Delta result:")
    old_time = perf_counter()
    old_delta = delta5_plus(l, (struct1, struct2), entry)
    old_time = perf_counter() - old_time
    print(old_delta)
    print("New Delta result:")
    new_time = perf_counter()
    new_delta = delta_plus_jies(l, (struct1, struct2), entry)
    new_time = perf_counter() - new_time
    print(new_delta)

    _, axs = plt.subplots(nrows=3, ncols=3)
    axs[0, 0].set_title("Original image", pad=20)
    axs[0, 0].matshow(entry)

    axs[0, 1].set_title("Old Delta", pad=20)
    axs[0, 1].matshow(old_delta)

    axs[0, 2].set_title("New Delta", pad=20)
    axs[0, 2].matshow(new_delta)

    axs[1, 0].set_title("S1 AND S2", pad=20)
    axs[1, 0].matshow(np.logical_and(struct1, struct2))

    axs[2, 0].set_title("Expected", pad=20)
    axs[2, 0].matshow(dilation3)

    axs[1, 1].set_title("S1", pad=20)
    axs[1, 1].matshow(struct1)

    axs[1, 2].set_title("S2", pad=20)
    axs[1, 2].matshow(struct2)

    axs[2, 1].set_title("D1", pad=20)
    axs[2, 1].matshow(dilation1)
    axs[2, 1].set_xlabel(old_time)

    axs[2, 2].set_title("D2", pad=20)
    axs[2, 2].matshow(dilation2)
    axs[2, 2].set_xlabel(new_time)

    plt.show()

def counter_example():
    shape = (4, 5)
    entry = np.zeros(shape, dtype=int)
    entry[1, 1] = 1
    entry[1, 2] = 1

    entry2 = np.zeros(shape, dtype=int)
    entry2[2, 2] = 1
    entry2[2, 3] = 1

    union = np.logical_or(entry, entry2)
    
    struct1 = ndimage.generate_binary_structure(2, 2)
    struct2 = ndimage.generate_binary_structure(2, 1)

    np.logical_and(
        ndimage.binary_dilation(union, structure=struct1),
        ndimage.binary_dilation(union, structure=struct2)
    )

    np.logical_or(
        np.logical_and(
            ndimage.binary_dilation(entry, structure=struct1),
            ndimage.binary_dilation(entry, structure=struct2)),
        np.logical_and(
            ndimage.binary_dilation(entry2, structure=struct1),
            ndimage.binary_dilation(entry2, structure=struct2)
        )
    )

if __name__ == "__main__":
    shape = (10, 10)
    entry = np.zeros(shape, dtype=int)
    entry[2, 2] = 1
    entry[2, 6] = 1
    entry[6, 4] = 1
    entry[6, 5] = 1
    entry[7, 4] = 1
    entry[7, 5] = 1
    struct1 = ndimage.generate_binary_structure(2, 2)
    struct2 = ndimage.generate_binary_structure(2, 1)
    run_dilations(entry, struct1, struct2)

    struct1 = np.array([[True, False, False], [True, True, False], [True, True, True]])
    struct2 = np.array([[False, False, True], [False, True, True], [True, True, True]])
    run_dilations2(entry, struct1, struct2)



import random
from itertools import combinations, product
from delta import calculate_lubs, calculate_glbs

# ######################################
# Lattice object, functions, and methods
# ######################################

class Lattice:
    def __init__(self, implies_matrix):
        """ Create a lattice from a matrix of implications,
            where implies_matrix[a][b] == 1, means that a >= b.

            Also calculates the corresponding matrix of
            least upper bounds and greatest lower bounds.
        """
        self.lattice = implies_matrix
        self.lubs = calculate_lubs(implies_matrix)
        self.glbs = calculate_glbs(implies_matrix)
        # Lazy attributes
        self._space_functions = None

    def __len__(self):
        """ Returns the number of nodes in the lattice.
        """
        return len(self.lattice)

    def is_lattice(self):
        """ Returns True if the lattice is valid.
        """
        N = len(self)
        # List of pairs, leaving them as an iterator allows us
        # to short-circuit and not generate all the pairs if
        # the lattice is not valid.
        lt = combinations(range(N),2)
        return all(self.lubs[a][b] != -1 for (a, b) in lt)
    
    def is_fn_distributive(self, fn, test_pairs):
        """ Checks if a given function `fn` is distributive.
        """
        fn = (0,) + fn
        return all(fn[self.lubs[t0][t1]] == self.lubs[fn[t0]][fn[t1]] for t0, t1 in test_pairs)

    def lub(self, iterable):
        """ Least Upper Bound of `iterable`.
        """
        r = 0
        for i in iterable:
            r = self.lubs[r][i]
        return r

    def glb(self, iterable):
        """ Greatest Lower Bound of `iterable`.
        """
        r = len(self) - 1
        for i in iterable:
            r = self.glbs[r][i]
        return r

    def space_functions(self):
        """ Return the list of space functions, based on the lattice.

            The actual `space_functions` are only generated once, feel free
            to call this method multiple times.
        """
        if self._space_functions is None:
           self._space_functions = self._generate_space_functions()
        return self._space_functions

    def _generate_space_functions(self):
        """ Generate a list of space functions, based on the lattice.
        """
        N = len(self)
        test_pairs = tuple(combinations(range(N), 2))
        return [(0,) + fn for fn in product(range(N), repeat=N-1) if self.is_fn_distributive(fn, test_pairs)]

# ######################################

# ######################################
# Helper functions for papers' algorithm
# ######################################

# This may be the wrong way of calculating S. ¯\_(ツ)_/¯
def S(i):
    """ Probably the log_2 of i aproximated to the ceilling, but allways >= 2.
    """
    j = 2
    while j * j < i:
        j += 1
    return j

def find_max(i):
    global M, Q, joins    
    k = 0
    for j in range(i):
        if all(a == j or joins[a][j] != a for a in range(i)):
            M[k] = j
            k += 1

    a = random.randrange(k) + 1
    
    for j in range(k):
        Q[j] = False
    
    s = 0
    while s < a:
        j = random.randrange(k)
        if Q[j]:
            s -= 1
        else:
            Q[j] = True
        s += 1
    
    return k

def work(i):
    global M, Q, joins

    N = len(joins)
    if i == N - 1:
        for j in range(N):
            for l in range(N):
                if joins[j][l] == -1:
                    joins[j][l] = N - 1
        return
    
    q = S(N-i)
    # Added to remove the papers' UB
    u = 0

    if i == 1:
        u = 1
        M[0] = 0
        Q[0] = True
    elif random.randrange(q) == 0:
        u = find_max(i)

    for j in range(u):
        if Q[j]:
            joins[M[j]][i] = joins[i][M[j]] = i
    
    w = True
    while w:
        w = False
        for j in range(i):
            if joins[j][i] == i:
                for s in range(i):
                    if joins[s][j] == j and joins[s][i] != i:
                        w = True
                        joins[s][i] = joins[i][s] = i
        for j in range(i):
            if joins[j][i] == i:
                for l in range(i):
                    if joins[l][i] == i:
                        s = joins[j][l]
                        if s != -1 and joins[s][i] != i:
                            w = True
                            joins[s][i] = joins[i][s] = i
    for j in range(i):
        if joins[j][i] == i:
            for l in range(i):
                if joins[l][i] == i and joins[j][l] == -1:
                    joins[j][l] = joins[l][j] = i

# ######################################

# #######################################
# Helper functions for lattice conversion
# #######################################

def lattice_from_joins(joins):
    """ Converts a join matrix of a lattice to an implies matrix
    """
    N = len(joins)

    lattice = [[0 for i in range(N)] for j in range(N)]

    for i in range(N):
        for j in range(i + 1):
            c = joins[i][j]
            lattice[c][i] = lattice[c][j] = 1
    
    return lattice

def explode(lc_all):
    r"""
    Helper function for lattice_from_covers.
    Returns a list of all elements below each element `i` of the list of
    covers `lc_all`.
    """
    n = len(lc_all)
    result = [set(i) for i in lc_all]
    exploded = [False for _ in range(n)]
    for i in range(n):
        exploded[i] = True
        covers = result[i].copy()
        while covers:
            cover = covers.pop()
            if not exploded[cover]:
                covers.update(result[cover])
                exploded[cover] = True
            result[i].update(result[cover])
    return result
            

def lattice_from_covers(lc_all):
    """ Converts a list of lower covers of a lattice to an implies matrix
    """
    N = len(lc_all)
    exploded = explode(lc_all)
    lattice = [[0]*N for _ in range(N)]
    print(exploded)
    for i in range(N):
        lattice[i][i] = 1
        for j in exploded[i]:
            lattice[i][j] = 1

    return lattice

# #######################################

def print_table(table):
    for row in table:
        print(*row)

# #######################################
# Functions for random lattice generation
# #######################################

def random_lattice(n, p):
    r"""
    Return a random lattice as a list of covers.
    
    Algorithm taken from:
    https://github.com/sagemath/sage/blob/master/src/sage/combinat/posets/poset_examples.py
    
    We add elements one by one. Check that adding a maximal
    element `e` to a meet-semilattice `L` with maximal elements
    `M` will create a semilattice by checking that there is a
    meet for `e, m` for all `m \in M`. We do that by keeping
    track of the meet matrix and the list of maximal elements.
    """
    from math import sqrt, floor

    n = n - 1
    meets = [[None]*n for _ in range(n)]
    meets[0][0] = 0
    maxs = set([0])
    lower_covers = [[]]

    for i in range(1, n):
        a = i - 1 - floor(i * sqrt(random.random()))
        lc_list = [a]
        maxs.discard(a)
        max_meets = {m:meets[m][a] for m in maxs}

        while random.random() < p and 0 not in lc_list:
            # Make number of coatoms closer to number of atoms.
            a = i - 1 - floor(i * sqrt(random.random()))

            # Check for antichain
            if any(meets[a][lc] in [a, lc] for lc in lc_list):
                continue
            
            # Check for unique meet with any maximal element and `a`
            for m in maxs:
                meet_m = meets[m][a]
                if meets[meet_m][max_meets[m]] not in [meet_m, max_meets[m]]:
                    break
            
            else: # We found a lower cover for i.
                lc_list.append(a)
                for m in maxs:
                    max_meets[m] = max(max_meets[m], meets[m][a])
                maxs.discard(a)
        
        # Compute a new row and column to the meet matrix
        meets[i][i] = i
        for lc in lc_list:
            meets[i][lc] = meets[lc][i] = lc
        for e in range(i):
            meets[i][e] = meets[e][i] = max(meets[e][lc] for lc in lc_list)

        maxs.add(i)
        lower_covers.append(lc_list)
    lower_covers.append(list(maxs)) # Add the top element.
    return lower_covers
 
def gen_lattice(n):
    assert(n > 0)

    # The C description uses global variable, will probably
    # move this to a context local varible shared between functions.
    global M, Q, joins

    # Initialize the join table and helper arrays
    joins = [[(i if i == j else -1) for i in range(n)] for j in range(n)]
    M = [0 for _ in range(n)]
    Q = [False for _ in range(n)]

    for i in range(1, n):
        work(i)
    
    # Fix the bottom connections
    for i in range(1,n - 1):
        joins[0][i] = joins[i][0] = i

    return joins

# #######################################

def process_file(path, gen_functions=False):
    """
    Process the file located in `path`, which should contain
    a list of join_matrices and convert them to a dictionary
    of implieas_matrices so it can be used by `delta.py`.

    Also allows the preprocessing of all the space functions
    associated with each lattice in the list, and saves the
    results in files with the corresponding key in the dictionary
    as the name of the file with the prefix "sf_".

    INPUTS:
    - `path` -- The location of the file relative to the script.
    - `gen_functions` -- If True, generates the corresponding
      space functions for each processed lattice and saves them
      in a file named with the equivalent key of the resulting
      dictionary, prefixed with "sf_" and ending in ".in".
      If False (default), do nothing else.
    """
    import os

    # Read list of matrices from a file.
    try:
        with open(path) as f:
            matrices = eval(f.read())
            print("[i] Reading input matrices from file")
    except IOError:
        print("[w] File not found `{}`, aborting processing...".format(path))
        return []


    # We have a list of matrices to process.
    # Prepare the dictionary for conversion.
    print("[i] Converting matrices")
    results = {hash(tuple(matrix)):lattice_from_joins(matrix) for submatrices in matrices for matrix in submatrices}

    # Now we create the corresponding space
    # functions for each matrix and save it
    # in a file.
    if gen_functions:
        for key, value in results.items():
            fns_file = "sf_{}.in".format(key)
            # Check if the file already exists
            if os.path.isfile(fns_file):
                print("[i] File `{}` already exist. Skipping.".format(fns_file))
            else:
                print("[i] Generating space functions for `{}` ({} nodes)".format(key, len(value)))
                lattice = Lattice(value)
                space_functions = lattice.space_functions()
                # Save the space functions to a file
                with open(fns_file, "w") as f:
                    f.write(repr(space_functions))
                    print("[i] Saved space functions in file `{}`".format(fns_file))

    return results

if __name__ == "__main__":

    process_file("distributive_lattices.py", True)


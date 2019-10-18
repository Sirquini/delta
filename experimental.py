import multiprocessing as mp
from itertools import combinations, product
from graphviz import Digraph
from lattice import covers_from_lattice, lattice_from_covers, calculate_lubs, calculate_glbs
from generation import is_distributive
from time import perf_counter

# ######################################
# Lattice object, functions, and methods
# ######################################

class NxLattice:
    def __init__(self, implies_matrix, parallel=False):
        """ Create a lattice from a matrix of implications,
            where implies_matrix[a][b] == 1, means that a >= b.

            Also calculates the corresponding matrix of
            least upper bounds and greatest lower bounds.
        """
        self.lattice = implies_matrix
        self.lubs = calculate_lubs(implies_matrix)
        self.glbs = calculate_glbs(implies_matrix)
        self.parallel = parallel # Added for parallel execution
        # Lazy attributes
        self._space_functions = None
        self._impls = None

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
        return all(fn[self.lubs[t0][t1]] == self.lubs[fn[t0]][fn[t1]] for t0, t1 in test_pairs)

    def is_fn_distributive_packed(self, pack):
        """ Equivalent to is_fn_distributive but all arguments are `pack`ed into one.
            It also returns the element checked and the result.
            Usefull for parallel computation using map, imap and equivalents.
        """
        fn, test_pairs = pack
        return (fn, self.is_fn_distributive(fn, test_pairs))

    def is_fn_monotonic(self, fn):
        """ Checks if a given function `fn` is monotonic.
        """
        n = len(self)
        return all(self.lattice[fn[a]][fn[b]] == 1 for a in range(n) for b in range(n) if self.lattice[a][b] == 1)

    def monotonic_functions(self):
        """ Generate all bottom-preserving monotonic functions of `lattice`.
        """
        n = len(self)
        return [(0,) + fn for fn in product(range(n), repeat=n-1) if self.is_fn_monotonic((0,) + fn)]

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

    def imply(self, a, b):
        """ Returns a imply b.
        """
        return self.impls[a][b]

    def atoms(self):
        """ Returns a list of all the atoms.
            
            + `y` is an atom if `0` is covered by `y`
            + `x` is covered by `y` if `x < y` and `x <= z < y` implies `z = x`
        """
        n = len(self)
        return [i for i in range(n) if all(i != 0 and (i == j or j == 0 or self.lattice[i][j] == 0) for j in range(n))]

    def join_irreducible_elements(self):
        """ Returns a list of all the join-irreducible elements in the lattice.

            `x` is join-irreducible (or lub-irreducible) if:
            + `x != 0`
            + `a < x` and `b < x` imply `a lub b < x` for all `a` and `b` 
        """
        n = len(self)
        test_pairs = tuple(combinations(range(n), 2))
        return [x for x in range(1, n) if all((x == a or x == b or x != self.lubs[a][b] for a, b in test_pairs))]

    @property
    def space_functions(self):
        """ Return the list of space functions, based on the lattice.

            The actual `space_functions` are only generated once, feel free
            to call this method multiple times.
        """
        if self._space_functions is None:
           self._space_functions = self._generate_space_functions()
        return self._space_functions

    @property
    def impls(self):
        """ Return the matrix of implications, based on the lattice.

            The actual `implys` are only generated once, feel free
            to call this method multiple times.
        """
        if self._impls is None:
           self._impls = self._calculate_implications()
        return self._impls

    def diagram(self, space_function=None):
        """ Returns the graphviz Digraph representation of the lattice for
            further manipulation or DOT language representation.
        """
        graph = Digraph("Lattice",edge_attr={"arrowhead": "none"})
        for i in range(len(self)):
            graph.node(str(i))
        for pos, nodes in enumerate(covers_from_lattice(self.lattice)):
            for node in nodes:
                graph.edge(str(pos), str(node))
        if space_function is not None:
            graph.attr("edge", arrowhead="normal", color="blue", constraint="false")
            for pos, val in enumerate(space_function):
                graph.edge(str(pos), str(val))
        return graph

    def _generate_space_functions(self):
        """ Generate a list of space functions, based on the lattice.
        """
        N = len(self)
        test_pairs = tuple(combinations(range(N), 2))
        if self.parallel:
            pool = mp.Pool(mp.cpu_count())
            result = [fn for fn, keep in pool.imap_unordered(self.is_fn_distributive_packed, (((0,) + f, test_pairs) for f in product(range(N), repeat=N-1)), 1000) if keep]
            pool.close()
            return result
        else:
            return [(0,) + fn for fn in product(range(N), repeat=N-1) if self.is_fn_distributive((0,) + fn, test_pairs)]

    def _calculate_implications(self):
        """ Calculate the matrix of implications for the lattice.
        """
        N = len(self)
        return [[self._pair_implication(i, j) for j in range(N)] for i in range(N)]
    
    def _pair_implication(self, a, b):
        """ Calculate the greatest lower bound of the pair (`a`, `b`)
            in the lattice
        """
        # a -> b ::= glb_{i <= b} { i | a lub i >= b }
        return self.glb((i for i in range(len(self)) if self.lattice[b][i] == 1 and self.lattice[self.lubs[a][i]][b] == 1))

# ######################################

if __name__ == "__main__":
    matrix = lattice_from_covers([[], [0], [0], [1,2]])
    normal_lattice = NxLattice(matrix)
    print("Is ditributive", is_distributive(normal_lattice))
    parallel_lattice = NxLattice(matrix, True)
    start = perf_counter()
    sf = len(normal_lattice.space_functions)
    finish = perf_counter() - start
    print("Normal lattice space functions ({}) time: {}".format(sf, finish))
    start = perf_counter()
    sf = len(parallel_lattice.space_functions)
    finish = perf_counter() - start
    print("Parallel lattice space functions ({}) time: {}".format(sf, finish))
    start = perf_counter()
    monotonic = parallel_lattice.monotonic_functions()
    finish = perf_counter() - start
    print("Parallel lattice monotonic functions ({}) time: {}".format(len(monotonic), finish))
    print(monotonic)

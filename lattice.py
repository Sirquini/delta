import os
import random
from itertools import combinations, product
from collections import deque
from graphviz import Digraph

# ######################################
#  Utility functions for handling files
# ######################################

def get_relative_path(file_path):
    """Returns the absolute path of a file relative to the script."""
    dirname = os.path.dirname(__file__)
    return os.path.join(dirname, file_path)

# ######################################

# ######################################
# Lattice object, functions, and methods
# ######################################

class Lattice:
    def __init__(self, implies_matrix):
        """Creates a lattice from a matrix of implications,
        where implies_matrix[a][b] == 1, means that a >= b.

        Also calculates the corresponding matrix of
        least upper bounds and greatest lower bounds.
        """
        self.lattice = implies_matrix
        self.lubs = calculate_lubs(implies_matrix)
        self.glbs = calculate_glbs(implies_matrix)
        # Lazy attributes
        self._space_functions = None
        self._impls = None
        self._join_irreducibles = None
        self._covers = None
        self._topological_order = None

    def __len__(self):
        """Returns the number of nodes in the lattice."""
        return len(self.lattice)
    
    @classmethod
    def from_covers(cls, covers):
        l = cls(lattice_from_covers(covers))
        l._covers = covers
        return l

    def is_lattice(self):
        """Returns True if the lattice is valid."""
        N = len(self)
        # List of pairs, leaving them as an iterator allows us
        # to short-circuit and not generate all the pairs if
        # the lattice is not valid.
        lt = combinations(range(N),2)
        return all(self.lattice[self.lubs[a][b]][a] == 1 and self.lattice[self.lubs[a][b]][b] == 1 for (a, b) in lt)
    
    def is_modular(self):
        """Returns True if the lattice is modular, False otherwise.

        A lattice is modular if, for all elements `a`, `b`, and `c`, it satisfies:
        
        + If  `a` >= `c`, then `a` glb (`b` lub `c`) = (`a` glb `b`) lub `c`.

        Warning: Not Tested!!
        """
        n = len(self)
        return all(self.glbs[a][self.lubs[b][c]] == self.lubs[self.glbs[a][b]][c] for c in range(n) for a in range(n) for b in range(n) if self.lattice[a][c] == 1)

    def is_fn_distributive(self, fn, test_pairs):
        """Checks if a given function `fn` is distributive."""
        return all(fn[self.lubs[t0][t1]] == self.lubs[fn[t0]][fn[t1]] for t0, t1 in test_pairs)

    @property
    def topological_order(self):
        """Returns a list with a linear order, from bottom to top, of the lattice.

        Comparable elements preserve their order. Incomparable elements are given
        an arbitrary order.

        The actual `topological_order` is only generated once, feel free
        to refer to this property multiple times.
        """
        if self._topological_order is None:
            self._topological_order = self._topological_sort()
        return self._topological_order

    def lub(self, iterable):
        """Least Upper Bound of `iterable`."""
        r = 0
        for i in iterable:
            r = self.lubs[r][i]
        return r

    def glb(self, iterable):
        """Greatest Lower Bound of `iterable`."""
        r = len(self) - 1
        for i in iterable:
            r = self.glbs[r][i]
        return r

    def imply(self, a, b):
        """Returns `a` imply `b`."""
        return self.impls[a][b]

    def atoms(self):
        """Returns a list of all the atoms.
            
        + `y` is an atom if `0` is covered by `y`
        + `x` is covered by `y` if `x < y` and `x <= z < y` implies `z = x`
        """
        n = len(self)
        return [i for i in range(n) if all(i != 0 and (i == j or j == 0 or self.lattice[i][j] == 0) for j in range(n))]

    @property
    def covers(self):
        """Returns a list of covers.
            
        + `y` is an atom if `0` is covered by `y`
        + `x` is covered by `y` if `x < y` and `x <= z < y` implies `z = x`

        The actual `covers` are only generated once, feel free
        to call this method multiple times.
        """
        if self._covers is None:
            self._covers = covers_from_lattice(self.lattice)
        return self._covers

    @property
    def join_irreducibles(self):
        """Returns a list of all the join-irreducible elements in the lattice.

        `x` is join-irreducible (or lub-irreducible) if:
        + `x != 0`
        + `a < x` and `b < x` imply `a lub b < x` for all `a` and `b`

        The actual `join_irreducibles` are only generated once, feel free
        to call this method multiple times.
        """
        if self._join_irreducibles is None:
            self._join_irreducibles = self.join_irreducible_elements()
        return self._join_irreducibles

    def join_irreducible_elements(self):
        """Returns a list of all the join-irreducible elements in the lattice.

        `x` is join-irreducible (or lub-irreducible) if:
        + `x != 0`
        + `a < x` and `b < x` imply `a lub b < x` for all `a` and `b` 
        """
        n = len(self)
        # Alternative if the list of covers is known
        if self._covers is not None:
            return [x for x, c in enumerate(self._covers) if len(c) == 1]
        else:
            test_pairs = tuple(combinations(range(n), 2))
            return [x for x in range(1, n) if all((x == a or x == b or x != self.lubs[a][b] for a, b in test_pairs))]

    @property
    def space_functions(self):
        """Returns the list of space functions valid for the lattice.

        The actual `space_functions` are only generated once, feel free
        to call this method multiple times.
        """
        if self._space_functions is None:
           self._space_functions = self._generate_space_functions()
        return self._space_functions

    @property
    def impls(self):
        """Returns the matrix of Heyting implications valid for the lattice.

        The actual `impls` are only generated once, feel free to call this
        method multiple times.
        """
        if self._impls is None:
           self._impls = self._calculate_implications()
        return self._impls

    def diagram(self, space_function=None):
        """Returns the graphviz Digraph representation of the lattice for
        further manipulation, or for DOT language representation.
        """
        graph = Digraph("Lattice",edge_attr={"arrowhead": "none"})
        for i in range(len(self)):
            graph.node(str(i))
        for pos, nodes in enumerate(self.covers):
            for node in nodes:
                graph.edge(str(pos), str(node))
        if space_function is not None:
            graph.attr("edge", arrowhead="normal", color="blue", constraint="false")
            for pos, val in enumerate(space_function):
                graph.edge(str(pos), str(val))
        return graph
    
    def space_functions_iter(self):
        """Returns an iterator of space functions, based on the lattice."""
        if self._space_functions is None:
            return iter(self._generate_space_functions_iter())
        return iter(self._space_functions)

    def _generate_space_functions_iter(self):
        """Returns an iterator of space functions, based on the lattice."""
        N = len(self)
        test_pairs = tuple(combinations(range(N), 2))
        for fn in product(range(N), repeat=N-1):
            result = (0,) + fn
            if self.is_fn_distributive(result, test_pairs):
                yield result


    def _generate_space_functions(self):
        """Generates a list of space functions, based on the lattice."""
        N = len(self)
        test_pairs = tuple(combinations(range(N), 2))
        return [(0,) + fn for fn in product(range(N), repeat=N-1) if self.is_fn_distributive((0,) + fn, test_pairs)]

    def _calculate_implications(self):
        """Calculates the matrix of implications for the lattice."""
        N = len(self)
        return [[self._pair_implication(i, j) for j in range(N)] for i in range(N)]
    
    def _topological_sort(self):
        """Calculates a topological ascending order for the lattice."""
        result = []
        n = len(self)
        
        if n == 1:
            return result
        work = deque((0,), n)

        # Covers are used implicitly
        edges = [len(cover) for cover in self.covers]
        while work:
            e = work.popleft()
            result.append(e)
            for node, covers in enumerate(self.covers):
                if e in covers:
                    edges[node] -= 1
                    if edges[node] == 0:
                        work.append(node)
        return tuple(result)
    
    def _pair_implication(self, a, b):
        """Calculates the Heyting implication of the pair (`a`, `b`)."""
        # a -> b ::= glb_{i <= b} { i | a lub i >= b }
        return self.glb((i for i in range(len(self)) if self.lattice[b][i] == 1 and self.lattice[self.lubs[a][i]][b] == 1))

# ########## Helper Functions ##########

def calculate_lubs(lattice):
    """Calculates the matrix of Least Upper Bounds for the `lattice`."""
    n = len(lattice)
    # Creates a decoding table
    encoding = {}
    for i in range(n):
        encoding[tuple(v[i] for v in lattice)] = i
    result = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(i+1):
            result[i][j] = pair_lub(i, j, lattice, encoding)
            result[j][i] = result[i][j]
    return result

def pair_lub(a, b, lattice, encoding):
    """Calculates the least upper bound of the pair (`a`, `b`)."""
    if lattice[a][b] == 1: return a
    elif lattice[b][a] == 1: return b
    else:
        n = len(lattice)
        entry = [row[a] * row[b] for row in lattice]
        return encoding.get(tuple(entry), n-1) # result or top

def calculate_glbs(lattice):
    """Calculates the matrix of Greatest Lower Bounds for the `lattice`."""
    n = len(lattice)
    # Creates a decoding table
    encoding = {}
    for i, v in enumerate(lattice):
        encoding[tuple(v)] = i
    result = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(i+1):
            result[i][j] = pair_glb(i, j, lattice, encoding)
            result[j][i] = result[i][j]
    return result

def pair_glb(a, b, lattice, encoding):
    """Calculates the greatest lower bound of the pair (`a`, `b`)."""
    if lattice[a][b] == 1: return b
    elif lattice[b][a] == 1: return a
    else:
        entry = [a * b for a, b in zip(lattice[a], lattice[b])]
        return encoding.get(tuple(entry), 0) # result or bottom

# ######################################

# ######################################
# Delta functions defined with the class
# Lattice and local context
# ######################################

class HelperCache:
    """Helper class for delta_ast_partition.

    Acts as a look-up table so data can be inserted/query.
    """

    def __init__(self, n, m):
        """Creates a HelperCache instance.

        Args:
            n: The size of the lattice.
            m: The number of space functions to evaluate.
        """
        self.cache = [[[None] * m for _ in range(m)] for _ in range(n)]
    
    def insert(self, c, first, last, result):
        self.cache[c][first][last] = result
    
    def get(self, c, first, last):
        return self.cache[c][first][last]

def partition_helper(lattice, functions, first, last, c, helper_cache):
    cached_result = helper_cache[c][first][last-1]
    if cached_result is not None:
        return cached_result
    if c == 0:
        return 0
    fn_num = last - first
    if fn_num == 1:
        return functions[first][c]
    else:
        n = len(lattice)
        mid_point = first + fn_num // 2
        result = lattice.glb(lattice.lub((partition_helper(lattice, functions, first, mid_point, a, helper_cache), partition_helper(lattice, functions, mid_point, last, lattice.imply(a, c), helper_cache))) for a in range(n) if lattice.lattice[c][a] == 1)
        helper_cache[c][first][last-1] = result
        return result

def delta_ast_partition(lattice, functions):
    """Calculates Delta* for a set of `functions` over a `lattice`
    partitioning the set of functions and using a look-up table.

    Args:
        lattice: A Lattice instance.
        functions: A list of space-functions.
    """
    n = len(functions)
    helper_cache = [[[None] * n for _ in range(n)] for _ in range(len(lattice))]
    return [partition_helper(lattice, functions, 0, n, c, helper_cache) for c in range(len(lattice))]

def delta_partition(lattice: Lattice, functions, jie_s=None):
    """Calculates Delta* for a set of `functions` over a `lattice`
    partitioning the set of functions and using a look-up table.

    This implementation takes advantage of the join irreducible
    elements to reduce the number of recursive calls.

    Args:
      lattice: A Lattice instance.
      functions: A list of space-functions.
      jie_s: A list of Join-Irreducible elements, if `None` it generates such list.
    """

    n = len(functions)
    helper_cache = [[[None] * n for _ in range(n)] for _ in range(len(lattice))]
    if jie_s is None:
        jie_s = lattice.join_irreducibles
    
    # Only call the recursive function for the join-irreducible elements
    result = [0 for _ in range(len(lattice))]
    for ji in jie_s:
        result[ji] = partition_helper(lattice, functions, 0, n, ji, helper_cache)
    for c in range(len(lattice)):
        if c not in jie_s:
            result[c] = lattice.lub(result[j] for j in jie_s if lattice.lattice[c][j] == 1)
    return result

class FooContext:
    """Helper class for delta_foo, containing:

    + (S) One sub-list of good pairs for each element in the lattice.
    + (C) All conflicting tuples from all elements in the lattice.
    + (R) A cross-referencing list of elements in the lattice.
    + (F) All pairs of elements in the lattice that lost support.
    """
    def __init__(self, n):
        """Creates a helper class for delta_foo, containing:

        + (S) One sub-list of good pairs for each element in the lattice.
        + (C) All conflicting tuples from all elements in the lattice.
        + (R) A cross-referencing list of elements in the lattice.
        + (F) All pairs of elements in the lattice that lost support.
        """
        self.good_pairs = [set() for _ in range(n)]
        self.conflicts = set()
        self.cross_references = [set() for _ in range(n)]
        self.falling_pairs = set()

    def process(self, lattice: Lattice, delta, u, v):
        """Sends the pair (u, v) to the set of conflicts (C), falling_pairs (F)
        or good_pairs (S), according to the property that holds for the pair.
        """
        w = lattice.lubs[u][v]
        if lattice.lattice[lattice.lubs[delta[u]][delta[v]]][delta[w]] != 1:
            self.conflicts.add(((u, v), w))
        elif lattice.lubs[delta[u]][delta[v]] == delta[w]:
            self.good_pairs[w].add((u, v))
            self.cross_references[u].add(v)
            self.cross_references[v].add(u)
        else:
            self.falling_pairs.add((u, v))
    
    def check_supports(self, lattice: Lattice, delta, u):
        """Identifies all pairs of the form (u, x) that lost their support
        because of a change in delta(u). It adds (u, x) to the appropiate set
        of conflicts (C), or falling_pairs (F).
        """
        for v in self.cross_references[u].copy():
            # v = self.cross_references[u].pop()
            w = lattice.lubs[u][v]
            if lattice.lubs[delta[u]][delta[v]] != delta[w]:
                self.good_pairs[w].discard((u, v))
                self.cross_references[u].discard(v)
                self.cross_references[v].discard(u)
                self.process(lattice, delta, u, v)

def delta_foo(lattice: Lattice, functions):
    """Calculates Delta using the Greatest Lower Bound between all the `functions`
    and then fixes the resulting function until it's a valid space-function.

    Args:
        lattice: A Lattice instance.
        functions: A list of space-functions.
    """
    n = len(lattice)
    # Here delta[c] = glb(fn_1[c], fn_2[c], ..., fn_n[c]), for fn_i in functions
    delta = [lattice.glb(i) for i in zip(*functions)]
    # Contains all delta_foo supporting structures (S, C, R, F)
    context = FooContext(n)
    # Calculate all initial conflicts in the candidate solution
    for u, v in combinations(range(n), 2):
        context.process(lattice, delta, u, v)

    while len(context.conflicts) != 0:
        (u, v), w = context.conflicts.pop()
        if lattice.lattice[lattice.lubs[delta[u]][delta[v]]][delta[w]] != 1:
            delta[w] = lattice.lub((delta[u], delta[v]))
            context.falling_pairs.update(context.good_pairs[w])
            context.good_pairs[w] = {(u, v)}

            context.check_supports(lattice, delta, w)
            context.cross_references[u].add(v)
            context.cross_references[v].add(u)   
        else:
            context.process(lattice, delta, u, v)
        
        while len(context.falling_pairs) != 0:
            x, y = context.falling_pairs.pop()
            z = lattice.lub((x, y))
            
            if lattice.lattice[delta[z]][lattice.lubs[delta[x]][delta[y]]] == 1:
                context.process(lattice, delta, x, y)
            else:
                if delta[x] != lattice.glb((delta[x], delta[z])):
                    delta[x] = lattice.glb((delta[x], delta[z]))
                    context.falling_pairs.update(context.good_pairs[x])
                    for u, v in context.good_pairs[x]:
                        context.cross_references[u].add(v)
                        context.cross_references[v].add(u)
                    context.good_pairs[x].clear()
                    context.check_supports(lattice, delta, x)
                
                if delta[y] != lattice.glb((delta[y], delta[z])):
                    delta[y] = lattice.glb((delta[y], delta[z]))
                    context.falling_pairs.update(context.good_pairs[y])
                    for u, v in context.good_pairs[y]:
                        context.cross_references[u].add(v)
                        context.cross_references[v].add(u)
                    context.good_pairs[y].clear()
                    context.check_supports(lattice, delta, y)

                if lattice.lub((delta[x], delta[y])) == delta[z]:
                    context.good_pairs[z].add((x, y))
                    context.cross_references[x].add(y)
                    context.cross_references[y].add(x)
                else:
                    context.conflicts.add(((x, y), z))
    return delta

def delta_foo_cvrs(lattice: Lattice, functions, covers=None):
    """Calculates Delta using the Greatest Lower Bound between all the `functions`
    and then fixes the resulting function until it's a valid space-function.

    This version of delta_foo only checks the pairs in the covers list, instead
    of all the pairs when initialicing its support structure.

    Args:
      lattice: A Lattice instance.
      functions: A list of space-functions.
      covers: A list of cover relations, if `None` it generates such list.
    """
    n = len(lattice)
    # Here delta[c] = glb(fn_1[c], fn_2[c], ..., fn_n[c]), for fn_i in functions
    delta = [lattice.glb(i) for i in zip(*functions)]
    # Contains all delta_foo supporting structures (S, C, R, F)
    context = FooContext(n)
    # Calculate the covers.
    if covers is None:
        covers = lattice.covers
    covers = [cvs + [node] for node, cvs in enumerate(covers)]
    # Calculate all initial conflicts in the candidate solution
    for cvs in covers:
        for a in range(len(cvs)):
            for b in range(a):
                context.process(lattice, delta, cvs[a], cvs[b])

    while len(context.conflicts) != 0:
        (u, v), w = context.conflicts.pop()
        if lattice.lattice[lattice.lubs[delta[u]][delta[v]]][delta[w]] != 1:
            delta[w] = lattice.lub((delta[u], delta[v]))
            context.falling_pairs.update(context.good_pairs[w])
            context.good_pairs[w] = {(u, v)}

            context.check_supports(lattice, delta, w)
            context.cross_references[u].add(v)
            context.cross_references[v].add(u)   
        else:
            context.process(lattice, delta, u, v)
        
        while len(context.falling_pairs) != 0:
            x, y = context.falling_pairs.pop()
            z = lattice.lub((x, y))
            
            if lattice.lattice[delta[z]][lattice.lubs[delta[x]][delta[y]]] == 1:
                context.process(lattice, delta, x, y)
            else:
                if delta[x] != lattice.glb((delta[x], delta[z])):
                    delta[x] = lattice.glb((delta[x], delta[z]))
                    context.falling_pairs.update(context.good_pairs[x])
                    for u, v in context.good_pairs[x]:
                        context.cross_references[u].add(v)
                        context.cross_references[v].add(u)
                    context.good_pairs[x].clear()
                    context.check_supports(lattice, delta, x)
                
                if delta[y] != lattice.glb((delta[y], delta[z])):
                    delta[y] = lattice.glb((delta[y], delta[z]))
                    context.falling_pairs.update(context.good_pairs[y])
                    for u, v in context.good_pairs[y]:
                        context.cross_references[u].add(v)
                        context.cross_references[v].add(u)
                    context.good_pairs[y].clear()
                    context.check_supports(lattice, delta, y)

                if lattice.lub((delta[x], delta[y])) == delta[z]:
                    context.good_pairs[z].add((x, y))
                    context.cross_references[x].add(y)
                    context.cross_references[y].add(x)
                else:
                    context.conflicts.add(((x, y), z))
    return delta

def delta_foo_jies(lattice, fns, jie_s=None):
    """Calculates Delta using the Greatest Lower Bound between all the `functions`
    and then fixes the resulting function until it's a valid space-function.

    This version of delta_foo only checks the pairs of join-irreducible elements,
    instead of all the pairs when initialicing its support structure.

    Args:
      lattice: A Lattice instance.
      functions: A list of space-functions.
      jie_s: A list of join-irreducible elements, if `None` it generates such list.
    """
    n = len(lattice)
    delta = [lattice.glb(i) for i in zip(*fns)]
    context = FooContext(n)
    processed = set()
    if jie_s is None:
        jie_s = lattice.join_irreducibles

    for a in range(len(jie_s)):
        for b in range(a):
            u, v = jie_s[a], jie_s[b]
            processed.add(u)
            processed.add(v)
            processed.add(lattice.lubs[u][v])
            context.process(lattice, delta, u, v)
    while len(context.conflicts) != 0:
        (u, v), w = context.conflicts.pop()
        if lattice.lattice[lattice.lubs[delta[u]][delta[v]]][delta[w]] != 1:
            delta[w] = lattice.lub((delta[u], delta[v]))
            context.falling_pairs.update(context.good_pairs[w])
            context.good_pairs[w] = {(u, v)}

            context.check_supports(lattice, delta, w)
            context.cross_references[u].add(v)
            context.cross_references[v].add(u)   
        else:
            context.process(lattice, delta, u, v)
        
        while len(context.falling_pairs) != 0:
            x, y = context.falling_pairs.pop()
            z = lattice.lub((x, y))
            
            if lattice.lattice[delta[z]][lattice.lubs[delta[x]][delta[y]]] == 1:
                context.process(lattice, delta, x, y)
            else:
                if delta[x] != lattice.glb((delta[x], delta[z])):
                    delta[x] = lattice.glb((delta[x], delta[z]))
                    context.falling_pairs.update(context.good_pairs[x])
                    for u, v in context.good_pairs[x]:
                        context.cross_references[u].add(v)
                        context.cross_references[v].add(u)
                    context.good_pairs[x].clear()
                    context.check_supports(lattice, delta, x)
                
                if delta[y] != lattice.glb((delta[y], delta[z])):
                    delta[y] = lattice.glb((delta[y], delta[z]))
                    context.falling_pairs.update(context.good_pairs[y])
                    for u, v in context.good_pairs[y]:
                        context.cross_references[u].add(v)
                        context.cross_references[v].add(u)
                    context.good_pairs[y].clear()
                    context.check_supports(lattice, delta, y)

                if lattice.lub((delta[x], delta[y])) == delta[z]:
                    context.good_pairs[z].add((x, y))
                    context.cross_references[x].add(y)
                    context.cross_references[y].add(x)
                else:
                    context.conflicts.add(((x, y), z))
    for i in range(n):
        if i not in processed:
            delta[i] = lattice.lub(delta[j] for j in jie_s if lattice.lattice[i][j] == 1)
    return delta

def delta_foo_cvrs_jies(lattice: Lattice, fns, covers=None, jie_s=None):
    """Calculates Delta using the Greatest Lower Bound between all the `functions`
    and then fixes the resulting function until it's a valid space-function.

    This version of delta_foo only checks the pairs of join-irreducibles from
    the same set of covers, instead of all the pairs, when initialicing its
    support structure.

    Args:
      lattice: A Lattice instance.
      functions: A list of space-functions.
      covers: A list of cover relations, if `None` it generates such list.
      jie_s: A list of join-irreducible elements, if `None` it generates such list.
    """
    n = len(lattice)
    delta = [lattice.glb(i) for i in zip(*fns)]
    context = FooContext(n)
    processed = set()
    
    # Calculate the covers.
    if covers is None:
        covers = lattice.covers
    covers = [cvs + [node] for node, cvs in enumerate(covers)]
    
    # Calculate the set of join-irredubles.
    if jie_s is None:
        jie_s = lattice.join_irreducibles
    
    # Process only pairs of join-irreducibles from the same set of covers.
    for a in range(len(jie_s)):
        for b in range(a):
            u, v = jie_s[a], jie_s[b]
            if any((u in cvrs) and (v in cvrs) for cvrs in covers):
                processed.add(u)
                processed.add(v)
                processed.add(lattice.lubs[u][v])
                context.process(lattice, delta, u, v)
    
    # Unchanged
    while len(context.conflicts) != 0:
        (u, v), w = context.conflicts.pop()
        if lattice.lattice[lattice.lubs[delta[u]][delta[v]]][delta[w]] != 1:
            delta[w] = lattice.lub((delta[u], delta[v]))
            context.falling_pairs.update(context.good_pairs[w])
            context.good_pairs[w] = {(u, v)}

            context.check_supports(lattice, delta, w)
            context.cross_references[u].add(v)
            context.cross_references[v].add(u)   
        else:
            context.process(lattice, delta, u, v)
        
        while len(context.falling_pairs) != 0:
            x, y = context.falling_pairs.pop()
            z = lattice.lub((x, y))
            
            if lattice.lattice[delta[z]][lattice.lubs[delta[x]][delta[y]]] == 1:
                context.process(lattice, delta, x, y)
            else:
                if delta[x] != lattice.glb((delta[x], delta[z])):
                    delta[x] = lattice.glb((delta[x], delta[z]))
                    context.falling_pairs.update(context.good_pairs[x])
                    for u, v in context.good_pairs[x]:
                        context.cross_references[u].add(v)
                        context.cross_references[v].add(u)
                    context.good_pairs[x].clear()
                    context.check_supports(lattice, delta, x)
                
                if delta[y] != lattice.glb((delta[y], delta[z])):
                    delta[y] = lattice.glb((delta[y], delta[z]))
                    context.falling_pairs.update(context.good_pairs[y])
                    for u, v in context.good_pairs[y]:
                        context.cross_references[u].add(v)
                        context.cross_references[v].add(u)
                    context.good_pairs[y].clear()
                    context.check_supports(lattice, delta, y)

                if lattice.lub((delta[x], delta[y])) == delta[z]:
                    context.good_pairs[z].add((x, y))
                    context.cross_references[x].add(y)
                    context.cross_references[y].add(x)
                else:
                    context.conflicts.add(((x, y), z))
    # Fix all the other elements
    for i in range(n):
        if i not in processed:
            delta[i] = lattice.lub(delta[j] for j in jie_s if lattice.lattice[i][j] == 1)
    return delta

# ######################################

# ######################################
# Helper functions for random lattice 
# generation for the papers' algorithm
# ######################################

# This may be the wrong way of calculating S. ¯\_(ツ)_/¯
def S(i):
    """Probably the ceilling of log_2 of i, but allways >= 2."""
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
    """Converts a join matrix of a lattice to an implies matrix."""
    N = len(joins)

    lattice = [[0 for i in range(N)] for j in range(N)]

    for i in range(N):
        for j in range(i + 1):
            c = joins[i][j]
            lattice[c][i] = lattice[c][j] = 1
    
    return lattice

def explode(lc_all):
    """Helper function for lattice_from_covers.

    Returns:
        A list of all elements below each element `i` of the list of covers
        `lc_all`.
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
    """Converts a list of lower covers of a lattice to an implies matrix."""
    N = len(lc_all)
    exploded = explode(lc_all)
    lattice = [[0]*N for _ in range(N)]
    for i in range(N):
        lattice[i][i] = 1
        for j in exploded[i]:
            lattice[i][j] = 1

    return lattice

def covers_from_lattice(matrix):
    """Converts an implies matrix into an equivalent list of covers."""
    n = len(matrix)
    result = [set() for _ in range(n)]
    for i, row in enumerate(matrix):
        for j, cell in enumerate(row):
            if cell == 1 and i != j:
                result[i].add(j)

    for i, lowers in enumerate(result):
        for value, elem in enumerate(result):
            if i in elem:
                result[value].difference_update(lowers)

    return [list(i) for i in result]

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

def powerset_lattice(n):
    """Generate the equivalent lattice for a powerset of 2^n elements.

    The generated lattice can be used with the `Lattice` class, or converted
    with the `covers_from_lattice` function.

    Args:
        n: The power of 2 for the number of elements. For example, powerset of
            `n = 3` produces the covers for a lattice of 2^3, or 8, elements. 
    """
    powerset = [set(elem) for i in range(n+1) for elem in combinations(range(1, n+1), i)]
    result = [[0 for _ in range(len(powerset))] for _ in range(len(powerset))]
    for i, row in enumerate(result):
        for j, _ in enumerate(row):
            if powerset[j].issubset(powerset[i]):
                result[i][j] = 1
    return result

# #######################################

def process_file(path, gen_functions=False):
    """Processes the file located in `path`, which should contain a list of
    join_matrices, and converts them to a dictionary of implieas_matrices so it
    can be used by `delta.py`.

    Also allows the preprocessing of all the space functions associated with
    each lattice in the list, and saves the results in files with the
    corresponding key in the dictionary as the name of the file with the prefix
    "sf_".

    Args:
        path: The location of the file relative to the script.
        gen_functions:
            - If True, generates the corresponding space functions for
            each processed lattice and saves them in a file named with the
            equivalent key of the resulting dictionary, prefixed with "sf_" and
            ending in ".in".
            - If False (default), do nothing else.
    """

    # Read list of matrices from a file.
    try:
        with open(get_relative_path(path)) as f:
            matrices = eval(f.read())
            print("[i] Reading input matrices from file")
    except IOError:
        print("[w] File not found `{}`, aborting processing...".format(path))
        return {}


    # We have a list of matrices to process.
    # Prepare the dictionary for conversion.
    print("[i] Converting matrices")
    results = {hash(tuple(matrix)):lattice_from_joins(matrix) for submatrices in matrices for matrix in submatrices}

    # Now we create the corresponding space
    # functions for each matrix and save it
    # in a file.
    if gen_functions:
        generated_dir = "generated"
        for key, value in results.items():
            fns_file_name = os.path.join(generated_dir, "sf_{}.in".format(key))
            fns_file_path = get_relative_path(fns_file_name)
            # Check if the file already exists
            if os.path.isfile(fns_file_path):
                print("[i] File `{}` already exist. Skipping.".format(fns_file_name))
            else:
                print("[i] Generating space functions for `{}` ({} nodes)".format(key, len(value)))
                lattice = Lattice(value)
                space_functions = lattice.space_functions
                # Save the space functions to a file
                with open(fns_file_path, "w") as f:
                    f.write(repr(space_functions))
                    print("[i] Saved space functions in file `{}`".format(fns_file_name))

    return results

def test_equality(expected, actual, name):
    message = "test {} ...".format(name)
    if expected != actual:
        print("{} \x1b[31mFAILED\x1b[0m".format(message))
        print("Expected:", expected)
        print("  Actual:", actual)
    else:
        print("{} \x1b[32mok\x1b[0m".format(message))

if __name__ == "__main__":

    # process_file("distributive_lattices.py", True)

    print("Running tests")
    # Testing covers_from_lattice
    expected = [[], [0], [0], [0], [1,2], [1,3], [2,3], [4,5,6]] 
    actual = covers_from_lattice(lattice_from_covers(expected))
    test_equality(expected, actual, "covers_from_lattice")

    # Testing Powerset generation
    expected = [[], [0], [0], [0], [1,2], [1,3], [2,3], [4,5,6]]
    actual = covers_from_lattice(powerset_lattice(3))
    test_equality(expected, actual, "Powerset(3) covers")

    expected = [[], [0], [0], [0], [0], [1,2], [1,3], [1,4], [2,3], [2,4], [3,4], [5,6,8], [5,7,9], [6,7,10], [8,9,10], [11,12,13,14]]
    actual = covers_from_lattice(powerset_lattice(4))
    test_equality(expected, actual, "Powerset(4) covers")

    space_function = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    diagram = Lattice(powerset_lattice(4)).diagram(space_function)
    diagram.render(directory="results",view=True)


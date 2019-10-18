import random
import os
import sys
from enum import Enum
from itertools import product
from itertools import combinations
from itertools import permutations
from time import perf_counter

def calculate_lubs(lattice):
    """Calculates the matrix of Lower Upper Bounds for the `lattice`."""
    n = len(lattice)
    result = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(i+1):
            result[i][j] = pair_lub(i, j, lattice)
            result[j][i] = result[i][j]
    return result

def calculate_glbs(lattice):
    """Calculates the matrix of Greatest Lower Bounds for the `lattice`."""
    n = len(lattice)
    result = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(i+1):
            result[i][j] = pair_glb(i, j, lattice)
            result[j][i] = result[i][j]
    return result

def pair_lub(a, b, lattice):
    """Calculates the least upper bound of the pair (`a`, `b`)."""
    if lattice[a][b] == 1: return a
    elif lattice[b][a] == 1: return b
    else:
        n = len(lattice)
        # Create the list of Upper Bounds of (a, b) 
        upper_bounds = (i for i in range(n) if lattice[i][a] == 1 and lattice[i][b] == 1)
        least = n - 1 # or False
        for ub in upper_bounds:
            if lattice[least][ub] == 1:
                least = ub
        return least

def pair_glb(a, b, lattice):
    """Calculates the greatest lower bound of the pair (`a`, `b`)."""
    if lattice[a][b] == 1: return b
    elif lattice[b][a] == 1: return a
    else:
        n = len(lattice)
        # Create the list of Lower Bounds of (a, b) 
        lower_bounds = (i for i in range(n) if lattice[a][i] == 1 and lattice[b][i] == 1)
        greatest = 0 # or True
        for lb in lower_bounds:
            if lattice[lb][greatest] == 1:
                greatest = lb
        return greatest

def lub(iterable):
    """Least Upper Bound of `iterable`.

    Using global variable `LUBs`
    """
    r = 0
    for i in iterable:
        r = LUBs[r][i]
    return r

def glb(iterable):
    """Greatest Lower Bound of `iterable`.

    Using global variable `GLBs`
    """
    r = len(GLBs) - 1
    for i in iterable:
        r = GLBs[r][i]
    return r

def calculate_implications(lattice):
    """Calculates the matrix of implications for the `lattice`."""
    n = len(lattice)
    return [[pair_implication(i, j, lattice) for j in range(n)] for i in range(n)]

def pair_implication(a, b, lattice):
    """Calculates the greatest lower bound of the pair (`a`, `b`).

    Implictly using the global variables `LUBs` and `GLBs`
    """
    # a -> b ::= glb_{i <= b} { i | a lub i >= b }
    return glb((i for i in range(len(lattice)) if lattice[b][i] == 1 and lattice[lub((a, i))][b] == 1))

def imply(a, b):
    """Returns `a` imply `b`.

    Using global variable `IMPLs`
    """
    return IMPLs[a][b]

def generate_functions(n):
    """Generates a list of space functions, based on a lattice of `n` elements.
    
    Implictly using the global variables `LUBs` and `GLBs`
    """
    test_pairs = tuple(combinations(range(n), 2))
    return [(0,) + fn for fn in product(range(n), repeat=n-1) if is_distributive(fn, test_pairs)]

def is_distributive(fn, test_pairs):
    """Checks if a given function `fn` is distributive.

    Implictly using the global variables `LUBs` and `GLBs`
    """
    fn = (0,) + fn
    return all(fn[lub((t0, t1))] == lub((fn[t0], fn[t1])) for t0, t1 in test_pairs)

def is_lattice(lattice):
    """Returns True if `lattice` is a valid lattice."""
    N = len(lattice)
	# List of pairs, leaving them as an iterator allows us
    # to short-circuit and not generate all the pairs if
    # the lattice is not valid.
    lt = combinations(range(N),2)
    return all(lub(pair) != -1 for pair in lt)

def leq_fn(lattice, fn1, fn2):
    """Returns True if the function `fn1` <= `fn2` for a given `lattice`.
    Otherwise False.
    """
    return all(lattice[fn2[i]][fn1[i]] == 1 for i in range(len(lattice)))

def max_fn(lattice, functions):
    """Returns the maximun function from a collection of `functions` for a
    given `lattice`.
    """
    greatest = None
    for fn in functions:
        if greatest is None:
            greatest = fn
        elif leq_fn(lattice, greatest, fn):
            greatest = fn
    return greatest


def apply_fns(functions, elements):
    """Returns a list of mapped elements applying the correponding function in
    `functions` to each element in `elements`.
    """
    return [functions[i][elem] for i, elem in enumerate(elements)]

def atoms(lattice):
    """Returns a list of all the atoms in `lattice`.
        
    + `y` is an atom if `0` is covered by `y`
    + `x` is covered by `y` if `x < y` and `x <= z < y` implies `z = x`
    """
    n = len(lattice)
    return [i for i in range(n) if all(i != 0 and (i == j or j == 0 or lattice[i][j] == 0) for j in range(n))]

def join_irreducible_elements(lattice):
    """Returns a list of all the join-irreducible elements in the `lattice`.

    `x` is join-irreducible (or lub-irreducible) if:
    + `x != 0`
    + `a < x` and `b < x` imply `a lub b < x` for all `a` and `b`

    Using the global variable `LUBs`
    """
    n = len(lattice)
    test_pairs = tuple(combinations(range(n), 2))
    return [x for x in range(1, n) if all((x == a or x == b or x != LUBs[a][b] for a, b in test_pairs))]


def random_space_function(lattice):
    """Generates a random space function for the given `lattice`.

    This mapping function is only guaranteed to be a valid space function
    if the `lattice` is isomorphic to a powerset lattice.

    Args:
        lattice: A Lattice instance.
    """
    n = len(lattice)
    all_atoms = lattice.atoms()
    # Map 0
    result = [0 for _ in range(n)]
    # Map all atoms to a random value
    for atom in all_atoms:
        result[atom] = random.randrange(0, n)
    # Map all other values
    for i in range(1, n):
        if i not in all_atoms:
            # f(c) = Vf(a_i) where Va_i = c and each a_i is an atom.
            result[i] = lattice.lub((result[atom] for atom in all_atoms if lattice.lattice[i][atom] == 1))
    return result

def decode_powerset_function(lattice, encoded_fn):
    """Expands an atoms-mapping to a space function.

    The atoms in `encoded_fn` are read in the same order as calling 
    `lattice.atoms()`.

    Args:
        lattice: A Lattice instance.
        encoded_fn: A list mapping atoms.
    """
    n = len(lattice)
    all_atoms = lattice.atoms()
    # Map 0
    result = [0 for _ in range(n)]
    # Map all atoms to the values of encoded_fn
    for pos, atom in enumerate(all_atoms):
        result[atom] = encoded_fn[pos]
    # Map all other values
    for i in range(1, n):
        if i not in all_atoms:
            # f(c) = Vf(a_i) where Va_i = c and each a_i is an atom.
            result[i] = lattice.lub((result[atom] for atom in all_atoms if lattice.lattice[i][atom] == 1))
    return result

def i_projection(lattice, function, c):
    """Returns the i-projection of `c` for the `lattice`.

    Defined as `lub({e | c >= function[e] and e in lattice})`

    Implicictly using the global variable `LUBs`
    """
    return lub(e for e in range(len(lattice)) if lattice[c][function[e]] == 1)

def i_join_projection(lattice, functions, c):
    """Returns the I-join projection of `c` for a group of `functions` over a
    `lattice`.

    Defined as `lub({i_projection(fn, c) | fn in functions})`

    Implictly using the global variable `LUBs`
    """
    return lub(i_projection(lattice, fn, c) for fn in functions)

def space_projection(lattice, functions):
    """Calculates the Space Projection, or I-join projection, for all elements
    of `lattice` with `functions`.

    Implictly using the global variable `LUBs`
    """
    return [i_join_projection(lattice, functions, c) for c in range(len(lattice))]

# ######################################
# Delta functions for lattice operations
# ######################################

def delta_ast(lattice, functions):
    """Calculates Delta* for a set of `functions` over a `lattice`.

    Implictly using the global variables `LUBs` and `GLBs`
    """
    n = len(lattice)
    fn_number = len(functions)
    delta = [0 for i in range(n)]
    for c in range(n):
        # All possible tuples to later calculate if tuple t > c
        tuples = product(range(n), repeat=fn_number)
        derives_c = (lub(apply_fns(functions, i)) for i in tuples if lattice[lub(i)][c] == 1)
        delta[c] = glb(derives_c)
    return delta

def delta_ast_imply(lattice, functions):
    """Calculate Delta* for a set of `functions` over a `lattice`.
        
    Only implemented with for 2 functions or less.

    Implictly using the global variables `LUBs`, `GLBs`, and `IMPLs`
    """
    n = len(lattice)
    # TODO: Define for cases with more than two functions
    if len(functions) == 1:
        return functions[0]
    elif len(functions) != 2:
        return None
    else:
        delta = [0 for _ in range(n)]
        fns = list(permutations(functions, len(functions)))
        for c in range(n):
            delta[c] = glb((lub((fn[0][i], fn[1][imply(i, c)])) for i in range(n) if lattice[c][i] == 1 for fn in fns))
        return delta

def partition_helper_v2(lattice, functions, first, last, c):
    """Helper function for delta*
        
    - No look-up table.
    - No imply

    Implictly using the global variables `LUBs` and `GLBs`
    """
    fn_num = last - first
    if fn_num == 1:
        return functions[first][c]
    else:
        n = len(lattice)
        mid_point = first + fn_num // 2
        result = glb(
            lub(
                (partition_helper_v2(lattice, functions, first, mid_point, a),
                 partition_helper_v2(lattice, functions, mid_point, last, b))
            ) for a, b in permutations(range(n), 2) if lattice[LUBs[a][b]][c] == 1
        )
        return result

def delta_ast_v2(lattice, functions):
    """Calculates Delta* for a set of `functions` over a `lattice`
    partitioning the set of functions.

    Similar to Delta* with partitioning but without using imply or a look-up
    table.

    Implictly using the global variables `LUBs` and `GLBs`
    """
    n = len(functions)
    return [partition_helper_v2(lattice, functions, 0, n, c) for c in range(len(lattice))]

def partition_helper_v3(lattice, functions, first, last, c):
    """Helper function for delta*
        
    - No look-up table.
    - With imply
    - Checks all elements, not only elements <= c

    Implictly using the global variables `LUBs`, `GLBs`, and `IMPLs`
    """
    fn_num = last - first
    if fn_num == 1:
        return functions[first][c]
    else:
        n = len(lattice)
        mid_point = first + fn_num // 2
        result = glb(lub((partition_helper_v3(lattice, functions, first, mid_point, a), partition_helper_v3(lattice, functions, mid_point, last, imply(a, c)))) for a in range(n))
        return result

def delta_ast_v3(lattice, functions):
    """Calculates Delta* for a set of `functions` over a `lattice`
    partitioning the set of functions.

    Implementation of Delta* without a look-up table.

    Implictly using the global variables `LUBs`, `GLBs`, and `IMPLs`
    """
    n = len(functions)
    return [partition_helper_v3(lattice, functions, 0, n, c) for c in range(len(lattice))]

def partition_helper_v4(lattice, functions, first, last, c):
    """Helper function for delta*
        
    - No look-up table.
    - With imply
    - Only checks elements <= c

    Implictly using the global variables `LUBs`, `GLBs`, and `IMPLs`
    """
    if c == 0:
        return 0
    fn_num = last - first
    if fn_num == 1:
        return functions[first][c]
    else:
        n = len(lattice)
        mid_point = first + fn_num // 2
        result = glb(lub((partition_helper_v4(lattice, functions, first, mid_point, a), partition_helper_v4(lattice, functions, mid_point, last, imply(a, c)))) for a in range(n) if lattice[c][a] == 1)
        return result

def delta_ast_v4(lattice, functions):
    """Calculates Delta* for a set of `functions` over a `lattice`
    partitioning the set of functions.

    Latest implementation of Delta* without a look-up table.

    Implictly using the global variables `LUBs`, `GLBs`, and `IMPLs`
    """
    n = len(functions)
    return [partition_helper_v4(lattice, functions, 0, n, c) for c in range(len(lattice))]

def partition_helper(lattice, functions, first, last, c):
    cached_result = HELPER_CACHE[c][first][last-1]
    if cached_result is not None:
        return cached_result
    # No need to compute for the bottom
    if c == 0:
        return 0
    fn_num = last - first
    if fn_num == 1:
        return functions[first][c]
    else:
        n = len(lattice)
        mid_point = first + fn_num // 2
        result = glb(lub((partition_helper(lattice, functions, first, mid_point, a), partition_helper(lattice, functions, mid_point, last, imply(a, c)))) for a in range(n) if lattice[c][a] == 1)
        HELPER_CACHE[c][first][last-1] = result
        return result

def delta_ast_partition(lattice, functions):
    """Calculates Delta* for a set of `functions` over a `lattice`
    partitioning the set of functions.

    Latest implementation of Delta* with a look-up table.

    Implictly using the global variables `LUBs`, `GLBs`, and `IMPLs`
    """
    n = len(functions)
    global HELPER_CACHE
    HELPER_CACHE = [[[None] * n for _ in range(n)] for _ in range(len(lattice))]
    return [partition_helper(lattice, functions, 0, n, c) for c in range(len(lattice))]

def overlapping_partition_helper(lattice, functions, first, last, c):
    """Helper function for delta*
        
    - With look-up table.
    - With imply
    - Only checks elements <= c
    - The partitions overlap

    Implictly using the global variables `LUBs`, `GLBs`, and `IMPLs`
    """
    cached_result = HELPER_CACHE[c][first][last-1]
    if cached_result is not None:
        return cached_result
    fn_num = last - first
    if fn_num == 1:
        return functions[first][c]
    else:
        n = len(lattice)
        mid_point = first + fn_num // 2
        left = right = mid_point
        # Adjust the mid-point so the partition overlaps
        if mid_point - first > 2:
            right += 1
        if last - mid_point > 2:
            left -= 1
        result = glb(lub((overlapping_partition_helper(lattice, functions, first, right, a), overlapping_partition_helper(lattice, functions, left, last, imply(a, c)))) for a in range(n) if lattice[c][a] == 1)
        HELPER_CACHE[c][first][last-1] = result
        return result

def delta_ast_overlapping_partition(lattice, functions):
    """Calculates Delta* for a set of `functions` over a `lattice`
    partitioning the set of functions.

    Using overlaping partitions, see `overlapping_partition_helper`, and a
    look-up table.

    Implictly using the global variables `LUBs`, `GLBs`, and `IMPLs`
    """
    n = len(functions)
    global HELPER_CACHE
    HELPER_CACHE = [[[None] * n for _ in range(n)] for _ in range(len(lattice))]
    return [overlapping_partition_helper(lattice, functions, 0, n, c) for c in range(len(lattice))]

def check_fn_with_pair(fn, pair):
    a, b = pair
    return LUBs[fn[a]][fn[b]] == fn[LUBs[a][b]]

def to_sets(lattice, delta, u, v, good_pairs, conflicts, cross_references, falling_pairs):
    """Sends pair (u, v) to the set of conflicts (C),  falling_pairs (F) or
    good_pairs (S) according to the property that holds for the pair.

    Uing the global variable `LUBs`
    """
    w = LUBs[u][v]
    if lattice[LUBs[delta[u]][delta[v]]][delta[w]] != 1:
        conflicts.add(((u, v), w))
    elif LUBs[delta[u]][delta[v]] == delta[w]:
        good_pairs[w].add((u, v))
        cross_references[u].add(v)
        cross_references[v].add(u)
    else:
        falling_pairs.add((u, v))

def check_supports(lattice, delta, u, good_pairs, conflicts, cross_references, falling_pairs):
    """Identifies all pairs of the form (u, x) that lost their support because
    of a change in delta(u). It adds (u, x) to the appropiate set of
    conflicts (C), or falling_pairs (F).

    Uing the global variable `LUBs`
    """
    for v in cross_references[u].copy():
        # v = cross_references[u].pop()
        w = LUBs[u][v]
        if LUBs[delta[u]][delta[v]] != delta[w]:
            good_pairs[w].discard((u, v))
            cross_references[u].discard(v)
            cross_references[v].discard(u)
            to_sets(lattice, delta, u, v, good_pairs, conflicts, cross_references, falling_pairs)

def delta_foo(lattice, functions):
    """Calculates Delta using Greatest Lower Bounds and then fixes the
    resulting function until its a valid space-function.

    This version uses supporting structure to speed-up the process. 

    Implictly using the global variables `LUBs` and `GLBs`
    """
    n = len(lattice)
    # Here delta[c] = glb(fn_1[c], fn_2[c], ..., fn_n[c]), for fn_i in functions
    delta = [glb(i) for i in zip(*functions)]
    # (S) One sub-list of good pairs for each element in the lattice
    good_pairs = [set() for _ in range(n)]
    # (C) All conflicting tuples form all elements in the lattice
    conflicts = set()
    # (R) Used for cross-refencing elements in the good_pairs list
    cross_references = [set() for _ in range(n)]
    # (F) All pairs of elements in the lattice that lost support
    falling_pairs = set()
    # Calculate all initial conflicts in the candidate solution
    for u, v in combinations(range(n), 2):
        to_sets(lattice, delta, u, v, good_pairs, conflicts, cross_references, falling_pairs)

    while len(conflicts) != 0:
        (u, v), w = conflicts.pop()
        if lattice[LUBs[delta[u]][delta[v]]][delta[w]] != 1:
            delta[w] = lub((delta[u], delta[v]))
            falling_pairs.update(good_pairs[w])
            good_pairs[w] = set([(u, v)])

            check_supports(lattice, delta, w, good_pairs, conflicts, cross_references, falling_pairs)
            cross_references[u].add(v)
            cross_references[v].add(u)   
        else:
            to_sets(lattice, delta, u, v, good_pairs, conflicts, cross_references, falling_pairs)
        
        while len(falling_pairs) != 0:
            x, y = falling_pairs.pop()
            z = lub((x, y))
            
            if lattice[delta[z]][LUBs[delta[x]][delta[y]]] == 1:
                to_sets(lattice, delta, x, y, good_pairs, conflicts, cross_references, falling_pairs)
            else:
                if delta[x] != glb((delta[x], delta[z])):
                    delta[x] = glb((delta[x], delta[z]))
                    falling_pairs.update(good_pairs[x])
                    for u, v in good_pairs[x]:
                        cross_references[u].add(v)
                        cross_references[v].add(u)
                    good_pairs[x].clear()
                    check_supports(lattice, delta, x, good_pairs, conflicts, cross_references, falling_pairs)
                
                if delta[y] != glb((delta[y], delta[z])):
                    delta[y] = glb((delta[y], delta[z]))
                    falling_pairs.update(good_pairs[y])
                    for u, v in good_pairs[y]:
                        cross_references[u].add(v)
                        cross_references[v].add(u)
                    good_pairs[y].clear()
                    check_supports(lattice, delta, y, good_pairs, conflicts, cross_references, falling_pairs)

                if lub((delta[x], delta[y])) == delta[z]:
                    good_pairs[z].add((x, y))
                    cross_references[x].add(y)
                    cross_references[y].add(x)
                else:
                    conflicts.add(((x, y), z))
    return delta

def delta_foo_b(lattice, functions):
    """Calculates Delta using Greatest Lower Bounds and then fixes the
    resulting function until its a valid space-function.

    Implictly using the global variables `LUBs` and `GLBs`.

    Alternate version of `delta_foo` used for A/B Testing, with a sanity
    check that imposes a performance penalty.
    """
    n = len(lattice)
    # Here candidate_function[c] = glb(fn_1[c], fn_2[c], ..., fn_n[c]), for fn_i in functions
    candidate_function = [glb(i) for i in zip(*functions)]
    # One sub-list of good pairs for each element in the lattice
    good_pairs = [set() for _ in range(n)]
    # All conflicting tuples form all elements in the lattice
    conflict_tuples = set()
    # All pairs of elements in the lattice that lost support
    falling_pairs = set()
    # Calculate all initial conflicts in the candidate solution
    for pair in combinations(range(n), 2):
        w = lub(pair)
        u, v = pair
        if check_fn_with_pair(candidate_function, pair):
            good_pairs[w].add(pair)
        elif lattice[lub((candidate_function[u], candidate_function[v]))][candidate_function[w]] != 1:
            conflict_tuples.add((pair, w))
        else:
            falling_pairs.add(pair)

    while len(conflict_tuples) != 0:
        (u, v), w = conflict_tuples.pop()
        candidate_update = lub((candidate_function[u], candidate_function[v]))
        candidate_function[w] = candidate_update
        falling_pairs.update(good_pairs[w])
        good_pairs[w] = set([(u, v)])

        while len(falling_pairs) != 0:
            x, y = falling_pairs.pop()
            z = lub((x, y))
            
            if candidate_function[x] != glb((candidate_function[x], candidate_function[z])):
                falling_pairs.update(good_pairs[x])
                good_pairs[x].clear()
            
            if candidate_function[y] != glb((candidate_function[y], candidate_function[z])):
                falling_pairs.update(good_pairs[y])
                good_pairs[y].clear()

            candidate_function[x] = glb((candidate_function[x], candidate_function[z]))
            candidate_function[y] = glb((candidate_function[y], candidate_function[z]))

            if lub((candidate_function[x], candidate_function[y])) == candidate_function[z]:
                good_pairs[z].add((x, y))
            else:
                conflict_tuples.add(((x, y), z))
        # Check that all good_pairs remain good. Sanity Check. Performance penalty.
        for sublist in good_pairs:
            for pair in sublist:
                if LUBs[candidate_function[pair[0]]][candidate_function[pair[1]]] != candidate_function[LUBs[pair[0]][pair[1]]]:
                    conflict_tuples.add((pair, LUBs[pair[0]][pair[1]]))
    return candidate_function

def probed_delta_foo(lattice, functions):
    """Calculates Delta using Greatest Lower Bounds and then fixes the
    resulting function until its a valid space-function.

    This function makes implicit use of the globals `LUBs` and `GLBs`.

    Alternate version of `delta_foo` with a sanity check that imposes
    a performance penalty. Returns a tuple with the result and the number
    of times the candidate delta function was updated, for testing.
    """
    # For testing
    updates = 0
    n = len(lattice)
    # Here candidate_function[c] = glb(fn_1[c], fn_2[c], ..., fn_n[c]), for fn_i in functions
    candidate_function = [glb(i) for i in zip(*functions)]
    # One sub-list of good pairs for each element in the lattice
    good_pairs = [set() for _ in range(n)]
    # All conflicting tuples form all elements in the lattice
    conflict_tuples = set()
    # All pairs of elements in the lattice that lost support
    falling_pairs = set()
    # Calculate all initial conflicts in the candidate solution
    for pair in combinations(range(n), 2):
        w = lub(pair)
        u, v = pair
        if check_fn_with_pair(candidate_function, pair):
            good_pairs[w].add(pair)
        elif lattice[lub((candidate_function[u], candidate_function[v]))][candidate_function[w]] != 1:
            conflict_tuples.add((pair, w))
        else:
            falling_pairs.add(pair)

    while len(conflict_tuples) != 0:
        (u, v), w = conflict_tuples.pop()
        candidate_update = lub((candidate_function[u], candidate_function[v]))
        if candidate_function[w] != candidate_update:
            updates += 1
        candidate_function[w] = candidate_update
        falling_pairs.update(good_pairs[w])
        good_pairs[w] = set([(u, v)])

        while len(falling_pairs) != 0:
            x, y = falling_pairs.pop()
            z = lub((x, y))
            
            if candidate_function[x] != glb((candidate_function[x], candidate_function[z])):
                falling_pairs.update(good_pairs[x])
                good_pairs[x].clear()
                updates += 1
            
            if candidate_function[y] != glb((candidate_function[y], candidate_function[z])):
                falling_pairs.update(good_pairs[y])
                good_pairs[y].clear()
                updates += 1

            candidate_function[x] = glb((candidate_function[x], candidate_function[z]))
            candidate_function[y] = glb((candidate_function[y], candidate_function[z]))

            if lub((candidate_function[x], candidate_function[y])) == candidate_function[z]:
                good_pairs[z].add((x, y))
            else:
                conflict_tuples.add(((x, y), z))
        # Check that all good_pairs remain good. Sanity Check. Performance penalty.
        for sublist in good_pairs:
            for pair in sublist:
                if LUBs[candidate_function[pair[0]]][candidate_function[pair[1]]] != candidate_function[LUBs[pair[0]][pair[1]]]:
                    conflict_tuples.add((pair, LUBs[pair[0]][pair[1]]))
    # print("Candidate Function updates:", updates)
    return (candidate_function, updates)

def delta_n(lattice, space_functions, functions):
    """ Calculate Delta for a set of `functions` over a `lattice`
        and a all possible `space_functions` for that lattice.
    """
    valid_functions = (fn for fn in space_functions if all(leq_fn(lattice, fn, elem) for elem in functions))
    return list(max_fn(lattice, (fn for fn in valid_functions)))

# ######################################

# ######################################
# Utility functions for lattice creation
# ######################################

def lattice_m(n):
    """Generates an M lattice with `n` elements."""
    l = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        l[i][0] = 1
        l[i][i] = 1
        l[n-1][i] = 1

    return l

def lattice_n5():
    """Generates the lattice n_5."""
    n = 5
    l = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        l[i][0] = 1
        l[i][i] = 1
        l[n-1][i] = 1
    l[3][1] = 1

    return l

def lattice_kite():
    """Generates a lattice with 4 elements."""
    return lattice_m(4)

def lattice_square():
    """Generates a square lattice with 9 elements."""
    n = 9
    l = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        l[i][0] = 1
        l[i][i] = 1
        l[n-1][i] = 1
    l[3][1] = 1
    l[4][1] = 1
    l[4][2] = 1
    l[5][2] = 1
    l[6][1] = 1
    l[6][2] = 1
    l[6][3] = 1
    l[6][4] = 1
    l[7][1] = 1
    l[7][2] = 1
    l[7][4] = 1
    l[7][5] = 1

    return l

def lattice_power_3():
    """Generates a powerset-of-3 lattice with 8 elements."""
    n = 8
    l = [[0 for i in range(n)] for i in range(n)]
    for i in range(n):
        l[i][0] = 1
        l[i][i] = 1
        l[n-1][i] = 1
    l[4][1] = 1
    l[4][2] = 1
    l[5][1] = 1
    l[5][3] = 1
    l[6][2] = 1
    l[6][3] = 1

    return l

# ######################################

# ######################################
# Utility functions for lattice checking
# ######################################

def print_fn_table_check(fn, verbose=False):
    """Prints a table with the values for each pair a b in a lattice with `n`
    nodes showing `a`, `b`, `lub(fn[a], fn[b])`, `fn[lub(a, b)]` for
    testing porpuses.
    
    This function makes implicit use of the global `LUBs`.
    """
    n = len(LUBs)

    if not verbose:
        print("[i] Showing only failing pairs")

    header = "a\t| b\t| fn[a] lub fn[b] | fn[a lub b]"
    print(header)
    print("-" * (len(header) + 10))

    valid = True
    for a, b in combinations(range(n), 2):
        fn_ab = fn[lub((a, b))]
        fn_a_fn_b = lub((fn[a], fn[b]))
        if verbose or fn_a_fn_b != fn_ab:
            print(a, b, fn_a_fn_b, fn_ab, sep='\t| ')
        if valid and fn_a_fn_b != fn_ab:
            valid = False
    # Just some fany terminal colors, only for Linux or macOS  systems
    if valid:
        color = "\x1b[32m"
    else:
        color = "\x1b[31m"

    print(fn, "Space Function:"+color, valid, "\x1b[0m\n")

# ######################################

# ######################################
#  Utility functions for handling files
# ######################################

def relative_path(path, *paths):
    """Returns the absolute path of a file relative to the script."""
    dirname = os.path.dirname(__file__)
    return os.path.join(dirname, path, *paths)

def eprint(*args, **kwargs):
    """Similar to `print` but uses `sys.stderr` instead of `sys.stdin`
    
    Also flushes the stream.
    """
    print(*args, file=sys.stderr, flush=True, **kwargs)

# ######################################

class TestResults:
    """Utility class to group and show algorithms' results."""
    
    def __init__(self, name, skip_times=False):
        # The name of the algorith tested
        self.name = name
        # Errors of the form (space_functions, bad_result, expected_result)
        self.errors = [] 
        self.avg_time = 0.0
        self.min_time = float("inf")
        self.max_time = 0.0
        # Used to indicate if timming values are not important
        self.skip_times = skip_times
    
    def update_times(self, new_time, max_iterations):
        """Updates the internal time counter with the `new_time`. Keeping
        track of the average (over `max_iterations`), minimum, and maximum
        time of the TestResult associated algorithm.
        """
        self.avg_time += new_time / max_iterations
        self.min_time = min(self.min_time, new_time)
        self.max_time = max(self.max_time, new_time)
    
    def print_times(self):
        """Prints the current average, minimum and maximum time."""
        print("\nTimes for", self.name)
        print("-- Average:", self.avg_time)
        print("-- Min:    ", self.min_time)
        print("-- Max:    ", self.max_time)
    
    def print_status(self):
        """Prints if the TestResults has failed any of the tests
        or passed all of them.
        """
        if len(self.errors) == 0:
            print("[i] All", self.name, "\x1b[32mpassed\x1b[0m")
        else:
            print("[w] Some", self.name, "\x1b[31mFAILED\x1b[0m")
    
    def print_errors(self):
        """If the TestResult presents any failed test cases, then it prints
        each failed test case.
        
        For each failed test case shows the space functions used, the
        expected (correct), and the actual (reported) result.
        """
        errors = len(self.errors)
        if errors > 0:
            print("\n[i] For", self.name, "showing", errors ,"errors:")
            for space_functions, actual, expected in self.errors:
                print("\n________For Space Functions________")
                for fn in space_functions:
                    print(fn)
                print("Expected:", expected)
                print("Actual:  ", actual)
                print("___________________________________")

class Delta(Enum):
    FOO = 0
    FOO_B = 1
    AST = 2
    AST_V2 = 3
    AST_V3 = 4
    AST_V4 = 5
    AST_LATEST = 6
    AST_OVERLAP = 7
    N = 8
    OTHER = 10

def run_test_case(fn, lattice, test_functions):
    """Runs `fn` with a `lattice` and an iterable of `test_functions`.
        
    Returns:
        A tuple with the execution_time and the actual result.
    """
    fn_time = perf_counter()
    result = fn(lattice, test_functions)
    fn_time = perf_counter() - fn_time
    return fn_time, result

def run(lattice, verbose = False, test_functions = None, n_tests = 100, n_functions = 2, fns_file = None, save_functions = False):
    """Runs all the delta algorithms against a `lattice` and `test_functions`,
        either random or explicit, and presents the results.

    Args:
        lattice:
            Matrix representing The lattice to test against.
        verbose:
            - If True, show each test case and the individual results,
            and a summary of failing test cases
            - If False (default), only show failing test cases (faster)
        test_functions:
            The list of space-functions to test against. By default, if no 
            `test_functions` are provided, generates random test_functions.
        n_tests:
            If `test_functions` is None, indicates the number of test to run
            with `n_functions` random test functions.
        n_functions:
            If `test_funtions` is None, indicates how many random test
            functions to use per test-case.
        fns_file:
            String with the file path to read and/or write calculated
            space-functions for `lattice`. 
        save_functions:
            - If True, write the generated space-function to `fns_file` except
            when the space-functions were read from the same file.
            - If False (default), do not write the generated space-functions
            anywhere.
        
    Returns:
        A dictionary with the acumulated results, including timings and
        failures (if any).
    """
    # Used to calculate elapsed time
    start_time = perf_counter()
    # Number of iterations to test (test-cases)
    n = 1

    preproc_time = perf_counter()
    # Matrix of least upper bounds and greatest lower bounds
    global LUBs, GLBs, IMPLs
    LUBs = calculate_lubs(lattice)
    GLBs = calculate_glbs(lattice)
    IMPLs = calculate_implications(lattice)
    preproc_time = perf_counter() - preproc_time

    # Valitdate the input lattice
    if not is_lattice(lattice):
        print("[e] Invalid lattice, aborting execution!")
        return {}

    func_gen_time = perf_counter()
    if fns_file is None:
        # Calculate space functions.
        space_functions = generate_functions(len(lattice))
    else:
        # Read space functions from a file
        try:
            # TODO: Validate that the space functions from the file correspond to the lattice
            with open(fns_file) as f:
                space_functions = eval(f.read())
                print("[i] Reading space functions from file")
                save_functions = False # Nothing changed no need to overwrite
        except IOError:
            print("[w] File not found, generating space functions...")
            # Calculate space functions.
            space_functions = generate_functions(len(lattice))
    
    validation_time = perf_counter()
    # By default use from 1 to `n`, inclusive, random test_functions
    # if none are provided, or check the test_functions
    # provided by the user.
    # Note: It was "temporarily" changed to signify the number of test-cases.
    if test_functions is None:
        n = n_tests
    else:
        valid_input = True
        for fn in test_functions:
            if fn not in space_functions:
                print("[E] Invalid Test Function found:", repr(fn))
                valid_input = False
        if not valid_input:
            print ("[E] Aborting execution!")
            return {}
    
    validation_time = perf_counter() - validation_time
    func_gen_time = perf_counter() - func_gen_time - validation_time

    print("Space functions:", len(space_functions))
    print("Space functions preprocessing time:   ", func_gen_time)
    print("LUBs, GLBs, IMPLYs preprocessing time:", preproc_time)
    print("________________________________________________________________\n")

    # Save the space functions to a file
    if save_functions and fns_file is not None:
        with open(fns_file, "w") as f:
            f.write(repr(space_functions))
            print("[i] Writing space functions to file `{}`".format(fns_file))

    # Used for showing the aggregate results at the end
    delta_results = {
        # Delta.AST: TestResults("Delta*"),
        Delta.FOO: TestResults("Delta_foo(V4)"),
        # Delta.FOO_B: TestResults("Delta_foo(V3)"),
        # Delta.AST_PART: TestResults("Delta*(OVERLAP)"),
        # Delta.AST_PART_B: TestResults("Delta*(PART+O2)"),
        Delta.N: TestResults("Delta_n(Preprocessed)", True)
    }

    for _ in range(n):
        # Get some space functions at random or use given ones
        if test_functions is None:
            sample_functions = random.sample(range(len(space_functions)), n_functions)
            sample_functions = [space_functions[j] for j in sample_functions]
        else:
            sample_functions = [fn for fn in test_functions if fn in space_functions]
        
        if verbose:
            print("__________Space Functions__________\n")
            for fn in sample_functions:
                print(fn)
            print("___________________________________\n")

        # fn_time, delta_ast_result = run_test_case(delta_ast, lattice, sample_functions)
        # delta_results[Delta.AST].update_times(fn_time, n)
        # if verbose:
        #     print("Delta*:          ", repr(delta_ast_result))
        #     print("-- Time:", fn_time, "\n")

        fn_time, delta_foo_result = run_test_case(delta_foo, lattice, sample_functions)
        delta_results[Delta.FOO].update_times(fn_time, n)
        if verbose:
            print("Delta_foo(V4):   ", repr(delta_foo_result))
            print("-- Time:", fn_time, "\n")
        
        # fn_time, delta_foo_b_result = run_test_case(delta_foo_b, lattice, sample_functions)
        # delta_results[Delta.FOO_B].update_times(fn_time, n)
        # if verbose:
        #     print("Delta_foo(V3):   ", repr(delta_foo_b_result))
        #     print("-- Time:", fn_time, "\n")
        #     print("Delta_0:         ", [glb(i) for i in zip(*sample_functions)], "\n")

        # fn_time, delta_ast_overlapping_result = run_test_case(delta_ast_overlapping_partition, lattice, sample_functions)
        # delta_results[Delta.AST_OVERLAP].update_times(fn_time, n)
        # if verbose:
        #     print("Delta*(OVERLAP): ", repr(delta_ast_overlapping_result))
        #     print("-- Time:", fn_time, "\n")

        # fn_time, delta_ast_part_result = run_test_case(delta_ast_partition, lattice, sample_functions)
        # delta_results[Delta.AST_LATEST].update_times(fn_time, n)
        # if verbose:
        #     print("Delta*(PART+O2): ", repr(delta_ast_part_result))
        #     print("-- Time:", fn_time, "\n")

        fn_time = perf_counter()
        delta_max_result = delta_n(lattice, space_functions, sample_functions)
        fn_time = perf_counter() - fn_time
        delta_results[Delta.N].update_times(fn_time, n)
        if verbose:
            print("Delta:           ", repr(delta_max_result))
            print("-- With preprocessing:   ", fn_time)
            print("-- Without preprocessing:", fn_time + func_gen_time)
            print("----------------------------------------------------------------\n")

        # If any of the algorithms failed to compute delta, add the error and context.
        # if delta_ast_result != delta_max_result:
        #     delta_results[Delta.AST].errors.append((sample_functions, delta_ast_result, delta_max_result))
        if delta_foo_result != delta_max_result:
            delta_results[Delta.FOO].errors.append((sample_functions, delta_foo_result, delta_max_result))
        # if delta_foo_b_result != delta_max_result:
        #     delta_results[Delta.FOO_B].errors.append((sample_functions, delta_foo_b_result, delta_max_result))
        # if delta_ast_overlapping_result != delta_max_result:
        #     delta_results[Delta.AST_OVERLAP].errors.append((sample_functions, delta_ast_overlapping_result, delta_max_result))
        # if delta_ast_part_result != delta_max_result:
        #     delta_results[Delta.AST_LATEST].errors.append((sample_functions, delta_ast_part_result, delta_max_result))

    print("Number of iterations:", n)

    # Show error test cases
    for result in delta_results.values():
        result.print_errors()

    # Show summary of results
    for result in delta_results.values():
        result.print_status()

    for result in delta_results.values():
        result.print_times()

    elapsed_time = perf_counter() - start_time
    print("\nTotal time:", elapsed_time)
    return {"result": delta_results.values(), "total": elapsed_time, "sf": len(space_functions), "sf_gen_time": func_gen_time, "preproc": preproc_time}

def write_test_results_csv(name, results):
    """Writes the test results in a CSV file."""
    import csv
    # Collect headers
    if results:
        headers = ["Lattice", "Nodes", "Space Functions", "Test Functions", "Elapsed Time", "Preprocessing Time"]
        for delta_result in results[0]["result"]:
            if not delta_result.skip_times:
                headers.append("{} Average".format(delta_result.name))
                headers.append("{} Min".format(delta_result.name))
                headers.append("{} Max".format(delta_result.name))
            headers.append("{} Errors".format(delta_result.name))
        with open(relative_path("results", name), 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            for result in results:
                row = [
                    result["lattice"],
                    result["nodes"],
                    result["sf"],
                    result["functions"],
                    result["total"],
                    result["preproc"]
                ]
                for delta_result in result["result"]:
                    if not delta_result.skip_times:
                        row.append(delta_result.avg_time)
                        row.append(delta_result.min_time)
                        row.append(delta_result.max_time)
                    row.append(len(delta_result.errors))
                writer.writerow(row)
        print("[i] Generated CSV file {}".format(name))

def run_full_tests():
    from lattice import process_file
    from datetime import datetime
    from generation import progress_bar
    # Read the lattice to test from
    lattices = process_file("distributive_lattices.py")
    results = []
    start_time = datetime.now().strftime("%Y-%m-%d-%H%M")
    count = 0
    for key, lattice in lattices.items():
        print("* Using lattice `{}` ({} nodes)".format(key, len(lattice)))
        count += 1
        progress_bar(count, len(lattices), 50)
        for i in range(2, 6):
            print("Test Functions:", i)
            result = run(lattice, n_functions=i, fns_file=relative_path("generated", "sf_{}.in".format(key)))
            result["lattice"] = key
            result["nodes"] = len(lattice)
            result["functions"] = i
            results.append(result)
        print("================================================================\n")
    eprint(" Done")
    write_test_results_csv("results-{}.csv".format(start_time), results)

def run_square():
    run(lattice_square(), n_tests=1000, n_functions=4, fns_file=relative_path("generated","sf_square.in"))

def run_failling_foo():
    """Test-cases in which delta_foo used to fail."""
    from lattice import process_file
    
    lattices = process_file("distributive_lattices.py")
    test_functions = [(0, 5, 3, 7, 3, 3, 7, 7),
        (0, 7, 7, 7, 0, 7, 7, 7),
        (0, 3, 4, 7, 6, 6, 7, 7)]

    run(lattices[4698136515449058355], test_functions=test_functions, fns_file=relative_path("generated", "sf_4698136515449058355.in"))

    test_functions = [(0, 0, 7, 7, 9, 9, 9, 9, 7, 9),
        (0, 3, 3, 3, 1, 3, 3, 3, 9, 9),
        (0, 7, 1, 7, 8, 8, 9, 9, 7, 9)]

    run(lattices[-286441500945297568], test_functions=test_functions, fns_file=relative_path("generated", "sf_-286441500945297568.in"))

    test_functions = [(0, 5, 6, 6, 8, 3, 8, 8, 8, 8),
        (0, 5, 4, 9, 9, 9, 9, 9, 9, 9),
        (0, 8, 5, 8, 8, 0, 5, 8, 8, 8),
        (0, 2, 8, 8, 8, 7, 8, 8, 8, 8),
        (0, 5, 4, 9, 9, 5, 9, 5, 9, 9)]

    run(lattices[8013726431884705816], test_functions=test_functions, fns_file=relative_path("generated", "sf_8013726431884705816.in"))
    
    test_functions = [(0, 8, 1, 8, 9, 5, 7, 8, 8, 9),
        (0, 4, 1, 4, 4, 9, 9, 9, 9, 9),
        (0, 6, 9, 9, 9, 3, 9, 8, 9, 9),
        (0, 7, 8, 8, 9, 1, 8, 7, 8, 9),
        (0, 1, 0, 1, 1, 9, 9, 9, 9, 9)]

    run(lattices[8013726431884705816], test_functions=test_functions, fns_file=relative_path("generated", "sf_8013726431884705816.in"))

    test_functions = [(0, 4, 8, 9, 9, 9, 9, 9, 9, 9),
        (0, 4, 7, 9, 9, 8, 8, 9, 9, 9),
        (0, 2, 4, 4, 9, 1, 4, 3, 4, 9),
        (0, 7, 4, 9, 9, 8, 9, 8, 9, 9),
        (0, 5, 2, 6, 8, 4, 4, 9, 9, 9)]

    run(lattices[8013726431884705816], test_functions=test_functions, fns_file=relative_path("generated", "sf_8013726431884705816.in"))

def run_random():
    """Equivalent to calling `run()` with a random lattice as its parameter."""
    from lattice import random_lattice, lattice_from_covers
    covers = random_lattice(8, 0.95)
    print("* Using lattice:", covers)
    lattice = lattice_from_covers(covers)
    run(lattice)

def run_lattice_implementations():
    from lattice import lattice_from_covers, Lattice
    covers = [[], [0], [0], [0], [2, 3, 1], [0], [2], [4, 5, 6]]
    print("* Using lattice:", covers)
    lattice_matrix = lattice_from_covers(covers)

    start_time = perf_counter()

    # Compare the LUBs, GLBs generation
    legacy_time = perf_counter()
    legacy_lubs = calculate_lubs(lattice_matrix)
    legacy_glbs = calculate_glbs(lattice_matrix)
    legacy_time = perf_counter() - legacy_time

    print("Legacy implementation LUBs, GLBs time:", legacy_time)

    # The new implementation encapsulates the LUBs and GLBs as part
    # of the Lattice object and are generated on initialization
    new_time = perf_counter()
    new_lattice = Lattice(lattice_matrix)
    new_time = perf_counter() - new_time

    print("New implementation LUBs, GLBs time:   ", new_time)

    # Compare results for sanity check
    if legacy_lubs != new_lattice.lubs or legacy_glbs != new_lattice.glbs:
        print("Some of the resulting matrices are \x1b[31mDIFFERENT\x1b[0m")
    else:
        print("All the resulting matrices are \x1b[32mequal\x1b[0m")

    # Needed for legacy implementation
    global LUBs, GLBs
    LUBs = legacy_lubs
    GLBs = legacy_glbs

    # Compare the IMPLs generation
    legacy_time = perf_counter()
    legacy_impls = calculate_implications(lattice_matrix)
    legacy_time = perf_counter() - legacy_time

    print("Legacy implementation IMPLs time:", legacy_time)

    # The new implementation encapsulates the LUBs and GLBs as part
    # of the Lattice object and are generated on initialization
    new_time = perf_counter()
    new_impls = new_lattice.impls
    new_time = perf_counter() - new_time

    print("New implementation IMPLs time:   ", new_time)

    # Compare results for sanity check
    if legacy_impls != new_impls:
        print("The resulting matrices are \x1b[31mDIFFERENT\x1b[0m")
    else:
        print("The resulting matrices are \x1b[32mequal\x1b[0m")

    # Compare space_functions generation
    legacy_time = perf_counter()
    legacy_space_fns = generate_functions(len(lattice_matrix))
    legacy_time = perf_counter() - legacy_time

    print("Legacy implementation space_functions time:", legacy_time)

    new_time = perf_counter()
    new_space_fns = new_lattice.space_functions
    new_time = perf_counter() - new_time

    print("New implementation space_functions time:   ", new_time)

    # Compare results for sanity check
    if legacy_space_fns != new_space_fns:
        print("Some of the resulting functions are \x1b[31mDIFFERENT\x1b[0m")
    else:
        print("All the resulting functions are \x1b[32mequal\x1b[0m")
    
    elapsed_time = perf_counter() - start_time
    print("\nElapsed time:", elapsed_time)

def run_random_space_functions():
    """Used for testing the random generation of space functions for a lattice.
    """
    from lattice import Lattice, lattice_from_covers

    covers = [[], [0], [0], [1], [1, 2], [2], [3, 4], [4], [4, 5], [6, 7], [6, 8], [7, 8], [9], [9, 10, 11], [11], [12, 13], [13, 14], [15, 16]]
    # covers = [[], [0], [0], [1], [1, 2], [2], [3], [3, 4], [4, 5], [5], [5], [6, 7, 8, 9, 10]]
    # covers = [[], [0], [0], [1, 2], [2], [3], [3, 4], [5,6]]
    # covers = [[], [0], [0], [0], [1,2], [1,3], [2,3], [4,5,6]]
    # covers = [[], [0], [1], [1], [2,3]]
    # covers = [[], [0], [0], [0], [1, 2, 3]]
    print("* Using lattice({}): {}".format(len(covers), covers))
    lattice_matrix = lattice_from_covers(covers)
    lattice = Lattice(lattice_matrix)

    print("\nAtoms for this lattice:", lattice.atoms())
    print("Join-irreducible elements for this lattice:", lattice.join_irreducible_elements())
    print("[i] Generating 100 random space functions.\n")
    space_fns = [random_space_function(lattice) for _ in range(100)]

    test_pairs = tuple(combinations(range(len(lattice)), 2))
    invalid_space_fns = [sf for sf in space_fns if not lattice.is_fn_distributive(sf, test_pairs)]

    if len(invalid_space_fns) != 0:
        print("Some ({}) of the random functions are \x1b[31mINVALID\x1b[0m".format(len(invalid_space_fns)))
        for i, v in enumerate(invalid_space_fns, 1):
            print("{} : {}".format(i, v))
    else:
        print("All the resulting functions are \x1b[32mvalid\x1b[0m")

def run_powerset(exponent = 10, verbose = False, test_functions = None, n_tests = 10, n_functions = 2):
    """Runs all the delta algorithms against a powerset lattice and
    `test_functions`, either random or explicit, similar to the run command,
    and present the results.

    Args:
        exponent:
            The number of base elements for the powerset, the total number of
            elements in the lattice is 2 ** n.
        verbose:
            - If `True`, show each test case, the individual results, and a
            summary of failing test cases.
            - If `False` (default), only show failing test cases.
        test_functions:
            The list of space-functions to test against. By default, if no 
            `test_functions` are provided, generates random test_functions.
        n_tests:
            If `test_functions` is None, indicates the number of test to run
            with `n_functions` random test functions.
        n_functions:
            If `test_funtions` is None, indicates how many random test
            functions to use per test-case.
        
    Returns:
        A dictionary with the acumulated results, including timings and
        failures (if any).
    """
    from lattice import Lattice, powerset_lattice
    # Used to calculate the elapsed time
    start_time = perf_counter()
    # Number of iterations to test (test-cases)
    n = 1
    # The actual number of elements in the lattice = 2^n
    elements = 2**exponent

    eprint("[i] Generating powerset lattice with {} nodes".format(elements))

    # Generate the powerset lattice and measure times
    gen_time = perf_counter()
    lattice_matrix = powerset_lattice(exponent)
    gen_time = perf_counter() - gen_time
    print("Lattice generation time:", gen_time)

    preproc_time = perf_counter()
    lattice = Lattice(lattice_matrix)
    # Generate the necessary globals
    global LUBs, GLBs, IMPLs
    LUBs = lattice.lubs
    GLBs = lattice.glbs
    IMPLs = lattice.impls
    preproc_time = perf_counter() - preproc_time

    if test_functions is None:
        n = n_tests
    else:
        valid_input = True
        test_pairs = tuple(combinations(range(len(lattice)), 2))
        for fn in test_functions:
            if not lattice.is_fn_distributive(fn, test_pairs):
                print("[E] Invalid Test Function found:", repr(fn))
                valid_input = False
        if not valid_input:
            print("[E] Aborting due to previous error!")
            return {}

    print("LUBs, GLBs, IMPLYs preprocessing time:", preproc_time)
    print("________________________________________________________________\n")
    
    # Used for showing the aggregate results at the end
    delta_results = {
        Delta.FOO: TestResults("Delta_foo(V4)"),
        # Delta.FOO_B: TestResults("Delta_foo(V3)"),
        # Delta.AST: TestResults("Delta+"),
        # Delta.AST_V2: TestResults("Delta++"),
        # Delta.AST_V3: TestResults("Delta+3"),
        # Delta.AST_V4: TestResults("Delta+4"),
        Delta.AST_LATEST: TestResults("Delta+5"),
    }
    deltas_are_equal = TestResults("assert_equal(Delta_foo(V4), Delta+5)", True)

    eprint("Running ", end='')
    for _ in range(n):
        # Get some space functions at random or use given ones
        if test_functions is None:
            sample_functions = [random_space_function(lattice) for _ in range(n_functions)]
        else:
            sample_functions = test_functions

        if verbose:
            print("__________Space Functions__________\n")
            for fn in sample_functions:
                print(fn)
            print("___________________________________\n")

        fn_time, delta_foo_result = run_test_case(delta_foo, lattice_matrix, sample_functions)
        delta_results[Delta.FOO].update_times(fn_time, n)
        if verbose:
            print("{}: {}".format(delta_results[Delta.FOO].name, repr(delta_foo_result)))
            print("-- Time: {}\n".format(fn_time))
        
        eprint(":", end='')

        # fn_time, delta_foo_b_result = run_test_case(delta_foo_b, lattice_matrix, sample_functions)
        # # fn_time, delta_foo_b_result = (0.0, [0])
        # delta_results[Delta.FOO_B].update_times(fn_time, n)
        # if verbose:
        #     print("{}: {}".format(delta_results[Delta.FOO_B].name, repr(delta_foo_b_result)))
        #     print("-- Time: {}\n".format(fn_time))
        
        # eprint(":", end='')

        # fn_time, delta_ast_result = run_test_case(delta_ast, lattice_matrix, sample_functions)
        # fn_time, delta_ast_result = (0.0, [0])
        # delta_results[Delta.AST].update_times(fn_time, n)
        # if verbose:
        #     print("{}: {}".format(delta_results[Delta.AST].name, repr(delta_ast_result)))
        #     print("-- Time: {}\n".format(fn_time))
        
        # eprint(":", end='')

        # fn_time, delta_ast_v2_result = run_test_case(delta_ast_v2, lattice_matrix, sample_functions)
        # fn_time, delta_ast_v2_result = (0.0, [0])
        # delta_results[Delta.AST_V2].update_times(fn_time, n)
        # if verbose:
        #     print("{}: {}".format(delta_results[Delta.AST_V2].name, repr(delta_ast_v2_result)))
        #     print("-- Time: {}\n".format(fn_time))
        
        # eprint(":", end='')

        # fn_time, delta_ast_v3_result = run_test_case(delta_ast_v3, lattice_matrix, sample_functions)
        # fn_time, delta_ast_v3_result = (0.0, [0])
        # delta_results[Delta.AST_V3].update_times(fn_time, n)
        # if verbose:
        #     print("{}: {}".format(delta_results[Delta.AST_V3].name, repr(delta_ast_v3_result)))
        #     print("-- Time: {}\n".format(fn_time))
        
        # eprint(":", end='')

        # fn_time, delta_ast_v4_result = run_test_case(delta_ast_v4, lattice_matrix, sample_functions)
        # fn_time, delta_ast_v4_result = (0.0, [0])
        # delta_results[Delta.AST_V4].update_times(fn_time, n)
        # if verbose:
        #     print("{}: {}".format(delta_results[Delta.AST_V4].name, repr(delta_ast_v4_result)))
        #     print("-- Time: {}\n".format(fn_time))
        
        # eprint(":", end='')

        fn_time, delta_ast_partition_result = run_test_case(delta_ast_partition, lattice_matrix, sample_functions)
        delta_results[Delta.AST_LATEST].update_times(fn_time, n)
        if verbose:
            print("{}: {}".format(delta_results[Delta.AST_LATEST].name, repr(delta_ast_partition_result)))
            print("-- Time: {}\n".format(fn_time))

        if delta_foo_result != delta_ast_partition_result:
            deltas_are_equal.errors.append((sample_functions, delta_foo_result, delta_ast_partition_result))

        eprint(".", end='')
    eprint("")

    print("Number of iterations:", n)

    deltas_are_equal.print_errors()
    deltas_are_equal.print_status()

    for result in delta_results.values():
        result.print_times()

    delta_results[Delta.OTHER] = deltas_are_equal

    elapsed_time = perf_counter() - start_time
    print("\nTotal time:", elapsed_time)
    return {"result": delta_results.values(), "total": elapsed_time, "sf": 0, "sf_gen_time": gen_time, "preproc": preproc_time}

def run_full_powerset_tests():
    """Equivalent to `run_full_tests()` but using `run_powerset()` instead of
    `run()`, and is used to time algorithms against powerset lattices.
    """
    from datetime import datetime
    results = []
    start_time = datetime.now().strftime("%Y-%m-%d-%H%M")
    for exponent in range(10, 11):
        nodes = 2**exponent
        print("* Using lattice `Powerset_{}` ({} nodes)".format(exponent, nodes))
        for i in [4, 8, 12, 16]: # TODO: Restore this to [4, 8, 12, 16]. Use only [4] for 4 and 5. Since Delta+ is O(mn^m), where m=n_functions 
            print("\nTest Functions:", i)
            result = run_powerset(exponent, n_tests=5, n_functions=i)
            result["lattice"] = "Powerset_{}".format(exponent)
            result["nodes"] = nodes
            result["functions"] = i
            results.append(result)
        print("================================================================\n")
    write_test_results_csv("results-{}.csv".format(start_time), results)

def run_space_projection():
    """Usecase for the space projection functions.

    + `i_projection()`
    + `i_join_projection()`
    + `space_projection()`
    """
    from lattice import Lattice, lattice_from_covers
    # 0 = bottom
    # 1 = a
    # 2 = b
    # 3 = top 
    covers = [[], [0], [0], [1, 2]]
    space_functions = [(0, 2, 1, 3), (0, 3, 2, 3)]
    lattice = Lattice(lattice_from_covers(covers))

    # Generate the necessary globals
    global LUBs, GLBs
    LUBs = lattice.lubs
    GLBs = lattice.glbs

    print("* Using lattice({}): {}".format(len(covers), covers))
    print("__________Space Functions__________\n")
    for fn in space_functions:
        print(fn)
    print("___________________________________\n")
    for c in range(len(lattice)):
        print("{}: {}".format(c, i_projection(lattice.lattice, space_functions[1], c)))
    print("Pi_1 a:", i_projection(lattice.lattice, space_functions[0], 2))
    print("Space Projection:", space_projection(lattice.lattice, space_functions))

def run_powerset_space_function_pairs():
    """Finds the pair of space functions such that the execution of delta_foo
    requires the maximum ammount of updates to the candidate_function.
    """
    from lattice import Lattice, covers_from_lattice, powerset_lattice
    # ID of the lattice
    key = "power_4"
    # Read the lattice to test from
    lattice = Lattice(powerset_lattice(4))
    
    # Globals used implicitly
    global LUBs, GLBs
    LUBs = lattice.lubs
    GLBs = lattice.glbs

    # Read space functions from a file
    try:
        with open(relative_path("generated", "sf_{}.in".format(key))) as f:
            space_functions = eval(f.read())
            print("[i] Reading space functions from file")
    except IOError:
        print("[w] File not found, generating space functions...")
        # Calculate space functions.
        space_functions = generate_functions(len(lattice))

    max_updates = 0
    results = []
    print("* Using lattice `{}` ({} nodes)".format(key, len(lattice)))
    print("* Covers:", covers_from_lattice(lattice.lattice))
    for pair in combinations(space_functions, 2):
        _, updates = probed_delta_foo(lattice.lattice, pair)
        if updates >= max_updates:
            results = pair
            max_updates = updates
    print("__________Space Functions__________\n")
    for fn in results:
        print(fn)
    print("___________________________________\n")
    print("Updates:", max_updates)
    print("================================================================\n")

def test_space_functions():
    from lattice import Lattice, powerset_lattice
    # covers = [[], [0], [0], [1, 2], [0], [1, 4], [2, 4], [3, 5, 6], [0], [1, 8], [2, 8], [3, 9, 10], [4, 8], [5, 9, 12], [6, 10, 12], [7, 11, 13, 14]]
    # lattice = Lattice(lattice_from_covers(covers))

    lattice = Lattice(powerset_lattice(5))

    encoded_space_functions = [
        (9, 30, 31, 27, 10),
        (10, 8, 0, 25, 9),
    ]

    # space_functions = [
    #     (0, 9, 14, 15, 15, 15, 15, 15, 11, 11, 15, 15, 15, 15, 15, 15),
    #     (0, 10, 9, 11, 0, 10, 9, 11, 12, 14, 13, 15, 12, 14, 13, 15)]

    # Globals used implicitly
    global LUBs, GLBs
    LUBs = lattice.lubs
    GLBs = lattice.glbs

    for sf in encoded_space_functions:
        sf = decode_powerset_function(lattice, sf)
        test_pairs = combinations(range(len(lattice)), 2)
        if lattice.is_fn_distributive(sf, test_pairs):
            print("\x1b[32mvalid\x1b[0m", sf)
        else:
            print("\x1b[31mINVALID\x1b[0m", sf)
    result, updates = probed_delta_foo(lattice.lattice, [decode_powerset_function(lattice, sf) for sf in encoded_space_functions])
    print("Updates:", updates)
    print("Result:", result)
    print("Delta_0:", [glb(i) for i in zip(*[decode_powerset_function(lattice, sf) for sf in encoded_space_functions])])

def run_space_functions_for_m():
    from lattice import Lattice    
    for i in range(9, 11):
        print("__\n* Using lattice M{} ({} nodes)".format(i - 2, i))
        lattice = Lattice(lattice_m(i))
        print("Space Functions:", len(lattice.space_functions))

def all_space_functions_diagram():
    from lattice import Lattice
    lattice = Lattice(lattice_m(6))
    total = len(lattice.space_functions)
    for pos, sf in enumerate(lattice.space_functions):
        file_name = "lattice_m{}_sf_{}.gv".format(pos - 2, pos)
        fns_file_path = relative_path("results", "diagrams", file_name)
        # Check if the file already exists
        if os.path.isfile(fns_file_path):
            print("[i] File `{}` already exist. Skipping.".format(file_name))
        else:
            print("[i] Generating space function {}/{}".format(pos + 1, total))
            # Save the diagram to a file
            lattice.diagram(sf).render(fns_file_path, cleanup=True)        

if __name__ == "__main__":
    # run_full_tests()
    # run_square()
    # run_random()
    # run_failling_foo()
    # run_lattice_implementations()
    # run_random_space_functions()
    # run_powerset(exponent=3, verbose=True, test_functions=[(0,4,5,6,7,7,7,7),(0,3,2,1,6,5,4,7)])
    # run_space_projection()
    run_full_powerset_tests()
    # run_powerset_space_function_pairs()
    # test_space_functions()
    # run_space_functions_for_m()
    # all_space_functions_diagram()

import random
import os
import sys
from enum import Enum
from itertools import product
from itertools import combinations
from itertools import permutations
from collections import deque
from time import perf_counter

from lattice import Lattice

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

def space_function(lattice):
    """Generates a random space function valid for the given `lattice`
    
    Args:
      lattice: A Lattice instance.
    
    Returns:
      A list representing a valid space function, only valid for the given
      `lattice` if it is distributive.
    """
    # Number of elements
    n = len(lattice)
    # The resulting space function mapped to None.
    result = [None for _ in range(n)]
    # The list of join-irreducible elements of the lattice
    ji_s = lattice.join_irreducible_elements()
    # Use a stack to know which Jis to update, preserving the monotonic property
    work = ji_s.copy()
    
    while work:
#         print("Working stack :", work)
        # Remove the last element from the work stack
        e = work.pop()
#         print(e, ": Checking")
        if result[e] is None:
            # Element unprocessed
            # Add dependencies
            dep = [j for j in ji_s if j != e and lattice.lattice[e][j] == 1]
#             print(e, ": Unprocessed")
            if len(dep) == 0:
                # No restriction
                result[e] = random.randrange(0,n)
#                 print(e, "Processed, Leaf")
            else:
                # check for unresolved dependencies
                dep_map = [result[d] for d in dep]
                if None in dep_map:
                    work.append(e)
                    work.extend([dep[i] for i, d in enumerate(dep_map) if d is None])
#                     print(e, ": Dependencies not met :", dep, dep_map)
                else:
                    # All dependencies covered, solve for e
                    result[e] = random.choice([i for i in range(n) if all(lattice.lattice[i][d] == 1 for d in dep_map)])
#                     print(e, ": Processed, all deps solved :", dep, dep_map)
    
    # Now map all the other elements
    for i in range(n):
        if i not in ji_s:
            result[i] = lattice.lub(result[j] for j in ji_s if lattice.lattice[i][j] == 1)
    return result

def decode_powerset_function(lattice, encoded_fn, atoms=None):
    """Expands an atoms-mapping to a space function.

    The atoms in `encoded_fn` are read in the same order as calling 
    `lattice.atoms()`.

    Args:
        lattice: A Lattice instance.
        encoded_fn: A list mapping atoms.
        atoms: The corresponding list of atoms, if already calculated.
    """
    n = len(lattice)
    if atoms is None:
        all_atoms = lattice.atoms()
    else:
        all_atoms = atoms
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

def all_space_functions(lattice):
    """Generates all space functions for the given `lattice`.

    This mapping function is only guaranteed generate valid space functions
    if the `lattice` is isomorphic to a powerset lattice.

    Args:
        lattice: A Lattice instance.
    """
    n = len(lattice)
    all_atoms = lattice.atoms()
    return [decode_powerset_function(lattice, fn, all_atoms) for fn in product(range(n), repeat=len(all_atoms))]

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
# Delta functions for Lattice instances
# ######################################

def delta_naive(lattice: Lattice, functions):
    """Calculates Delta+ for a set of `functions` over a `lattice`.

    Args:
      lattice: A Lattice instance.
      functions: A list of space functions.
    """
    n = len(lattice)
    return [lattice.glb(
        lattice.lub(apply_fns(functions, t)) for t in product(range(n), repeat=len(functions))
        if lattice.lattice[lattice.lub(t)][c] == 1
    ) for c in range(n)]

def delta_plus(lattice: Lattice, functions):
    """Calculates Delta+ for a set of `functions` over a `lattice`.

    Args:
      lattice: A Lattice instance.
      functions: A list of space functions.
    """
    fn_len = len(functions)
    if fn_len == 1:
        return functions[0]
    if fn_len == 2:
        n = len(lattice)
        return [lattice.glb(
            lattice.lub(apply_fns(functions, t)) for t in product(range(n), repeat=len(functions))
            if lattice.lattice[lattice.lub(t)][c] == 1
        ) for c in range(n)]
    return delta_plus(lattice, (delta_plus(lattice, functions[:2]), delta_plus(lattice, functions[2:])))

def delta_plus_imply(lattice: Lattice, functions):
    """Calculate Delta+ for a set of `functions` over a `lattice`.
        
    Similar to `delta_ast_partition`.

    - No look-up table.
    - With imply
    - Only checks elements <= c

    Args:
      lattice: A Lattice instance.
      functions: A list of space functions.
    """
    if len(functions) == 1:
        return functions[0]
    if len(functions) == 2:
        n = len(lattice)
        return [lattice.glb(
            lattice.lub(apply_fns(functions, (i, lattice.imply(i, c)))) for i in range(n) if lattice.lattice[c][i] == 1
        ) for c in range(n)]
    return delta_plus_imply(lattice, (delta_plus_imply(lattice, functions[:2]), delta_plus_imply(lattice, functions[2:])))

def delta_plus_jies(lattice: Lattice, functions, jie_s=None):
    """Calculate Delta+ for a set of `functions` over a `lattice`.
    
    Similar to `delta_partition`.

    This implementation takes advantage of the join irreducible
    elements to reduce the number of operations.

    Args:
      lattice: A Lattice instance.
      functions: A list of space functions.
      jie_s: A list of Join-Irreducible elements, if `None` it generates such list.
    """
    if len(functions) == 1:
        return functions[0]
    if jie_s is None:
        jie_s = lattice.join_irreducibles
    n = len(lattice)
    result = [None for _ in range(n)]
    for j in jie_s:
        result[j] = lattice.glb(fn[j] for fn in functions)
    return [lattice.lub(result[j] for j in jie_s if lattice.lattice[c][j] == 1) if v is None else v for c, v in enumerate(result)]

def dmeet_jies(lattice, functions, jie_s=None, covers=None):
    """ Like `delta_plus_jies`.

    This implementation takes advantage of the join irreducible
    elements to reduce the number of operations.

    Args:
      lattice: A Lattice instance.
      functions: A list of space functions.
      jie_s: A list of Join-Irreducible elements. If `None`, it generates such list.
      covers: A list of cover relations. If `None`, it generates such list.
    """
    fn_n = len(functions)
    if fn_n == 1:
        return functions[0]
    if covers is None:
        covers = lattice.covers
    if jie_s is None:
        jie_s = lattice.join_irreducibles
    n = len(lattice)
    result = [None for _ in range(n)]
    # Mark bottom
    result[0] = 0
    # Mark ji
    for j in jie_s:
        result[j] = lattice.glb(fn[j] for fn in functions)
    # Solve for the rest
    work = deque((c for c,v in enumerate(result) if v is None))
    while work:
        # Remove the last element from the stack
        e = work.pop()
        if result[e] is None:
            if result[covers[e][0]] is None or result[covers[e][1]] is None:
                work.append(e)
                if result[covers[e][0]] is None:
                    work.append(covers[e][0])
                if result[covers[e][1]] is None:
                    work.append(covers[e][1])
            else:
                result[e] = lattice.lub((result[covers[e][0]], result[covers[e][1]]))
    return result

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
    PLUS = 8
    N = 9
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

    # Global used by legacy delta_ast_ implementations
    global IMPLs
    IMPLs = legacy_impls 

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
    from lattice import Lattice

    covers = [[], [0], [0], [1], [1, 2], [2], [3, 4], [4], [4, 5], [6, 7], [6, 8], [7, 8], [9], [9, 10, 11], [11], [12, 13], [13, 14], [15, 16]]
    # covers = [[], [0], [0], [1], [1, 2], [2], [3], [3, 4], [4, 5], [5], [5], [6, 7, 8, 9, 10]]
    # covers = [[], [0], [0], [1, 2], [2], [3], [3, 4], [5,6]]
    # covers = [[], [0], [0], [0], [1,2], [1,3], [2,3], [4,5,6]]
    # covers = [[], [0], [1], [1], [2,3]]
    # covers = [[], [0], [0], [0], [1, 2, 3]]
    print("* Using lattice({}): {}".format(len(covers), covers))
    lattice = Lattice.from_covers(covers)

    print("\nAtoms for this lattice:", lattice.atoms())
    print("Join-irreducible elements for this lattice:", lattice.join_irreducibles)
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

def run_space_projection():
    """Usecase for the space projection functions.

    + `i_projection()`
    + `i_join_projection()`
    + `space_projection()`
    """
    from lattice import Lattice
    # 0 = bottom
    # 1 = a
    # 2 = b
    # 3 = top 
    covers = [[], [0], [0], [1, 2]]
    space_functions = [(0, 2, 1, 3), (0, 3, 2, 3)]
    lattice = Lattice.from_covers(covers)

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
    from lattice import Lattice, powerset_lattice
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
    print("* Covers:", lattice.covers)
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
    # run_lattice_implementations()
    # run_random_space_functions()
    # run_space_projection()
    # run_powerset_space_function_pairs()
    # test_space_functions()
    # run_space_functions_for_m()
    # all_space_functions_diagram()
    pass

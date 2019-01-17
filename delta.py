import random

from enum import Enum
from itertools import product
from itertools import combinations
from itertools import permutations
from time import time

def calculate_lubs(lattice):
    """ Calculate the matrix of Lower Upper Bounds for the `lattice`
    """
    n = len(lattice)
    result = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(i+1):
            result[i][j] = pair_lub(i, j, lattice)
            result[j][i] = result[i][j]
    return result

def calculate_glbs(lattice):
    """ Calculate the matrix of Greatest Lower Bounds for the `lattice`
    """
    n = len(lattice)
    result = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(i+1):
            result[i][j] = pair_glb(i, j, lattice)
            result[j][i] = result[i][j]
    return result

def pair_lub(a, b, lattice):
    """ Calculate the lower upper bound of the pair (`a`, `b`)
        given `lattice`
    """
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
    """ Calculate the greatest lower bound of the pair (`a`, `b`)
        given `lattice`
    """
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
    """ Least Upper Bound of `iterable`.
        Using global variable `LUBs`
    """
    r = 0
    for i in iterable:
        r = LUBs[r][i]
    return r

def glb(iterable):
    """ Greatest Lower Bound of `iterable`.
        Using global variable `GLBs`
    """
    r = len(GLBs) - 1
    for i in iterable:
        r = GLBs[r][i]
    return r

def calculate_implications(lattice):
    """ Calculate the matrix of implications for the `lattice`
    """
    n = len(lattice)
    result = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            result[i][j] = pair_implication(i, j, lattice)
    return result

def pair_implication(a, b, lattice):
    """ Calculate the greatest lower bound of the pair (`a`, `b`)
        given `lattice`
        Implictly using the global variables `LUBs` and `GLBs`
    """
    # a -> b ::= glb_{i <= b} { i | a lub i <= b }
    return glb((i for i in range(len(lattice)) if lattice[b][i] == 1 and lattice[lub((a, i))][b] == 1))

def imply(a, b):
    """ Returns a impply b.
        Using global variable `IMPLs`
    """
    return IMPLs[a][b]


def generate_functions(n):
    """ Generate an list of space functions, based on a lattice
        of `n` elements.
    """
    test_pairs = list(combinations(range(n), 2))
    return [(0,) + fn for fn in product(range(n), repeat=n-1) if is_distributive(fn, test_pairs)]

def is_distributive(fn, test_pairs):
    """ Checks if a given function `fn` is distributive.
    """
    fn = (0,) + fn
    for i in range(len(test_pairs)):
        t0 = test_pairs[i][0]
        t1 = test_pairs[i][1]
        if fn[lub([t0, t1])] != lub([fn[t0], fn[t1]]):
            return False
    return True

def leq_fn(lattice, fn1, fn2):
    """ Returns True if the function `fn1` < `fn2`
        for a given `lattice`. Otherwise False.
    """
    for i in range(len(lattice)):
        if lattice[fn2[i]][fn1[i]] != 1:
            return False
    return True

def max_fn(lattice, functions):
    """ Returns the maximun function from a collection of
        `functions` for a given `lattice`.
    """
    greatest = None
    for fn in functions:
        if greatest is None:
            greatest = fn
        elif leq_fn(lattice, greatest, fn):
            greatest = fn
    return greatest


def apply_fns(functions, elements):
    """ Returns a list of mapped elements applying the
        correponding function in `functions` to each element in `elements`
    """
    return [functions[i][elem] for i, elem in enumerate(elements)]

# ######################################
# Delta functions for lattice operations
# ######################################

def delta_ast(lattice, functions):
    """ Calculate Delta* for a set of `functions` over a `lattice`.
    """
    n = len(lattice)
    fn_number = len(functions)
    # All possible tuples to later calculate if tuple t > c
    tuples = list(product(range(n), repeat=fn_number))
    delta = [0 for i in range(n)]
    for c in range(n):
        derives_c = (lub(apply_fns(functions, i)) for i in tuples if lattice[lub(i)][c] == 1)
        delta[c] = glb(derives_c)
    return delta

def delta_ast_imply(lattice, functions):
    """ Calculate Delta* for a set of `functions` over a `lattice`.
        Using global variable `IMPLs`
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

def delta_ast2(lattice, functions):
    """ Calculate Delta* for a set of `functions` over a `lattice`.
    """
    fn_num = len(functions)
    if fn_num == 1:
        return functions[0]
    elif fn_num == 2:
        return delta_ast_imply(lattice, functions)
    else:
        return delta_ast2(lattice, [delta_ast(lattice, functions[:2])] + functions[2:])

def partition_helper(lattice, functions, c):
    # TODO: Define for len(functions) != log_2 k, where k is a natural number.
    n = len(lattice)
    fn_num = len(functions)
    if fn_num == 1:
        return functions[0][c]
    else:
        mid_point = fn_num // 2
        pairs = (pair for pair in product(range(n), repeat=2) if lattice[lub(pair)][c] == 1)
        return glb(lub((partition_helper(lattice, functions[:mid_point], a), partition_helper(lattice, functions[mid_point:], b))) for (a, b) in pairs)

def delta_ast_partition(lattice, functions):
    """ Calculate Delta* for a set of `functions` over a `lattice`
        partitioning the set of functions.
    """
    return [partition_helper(lattice, functions, c) for c in range(len(lattice))]


def check_fn_with_pair(fn, pair):
    a, b = pair
    return lub((fn[a], fn[b])) == fn[lub((a, b))]

def delta_foo_old(lattice, functions):
    """ Calculates Delta using Greatest Lower Bounds and then fixes the
        resulting function until its a valid space-function.
        This function makes implicit use of the globals `LUBs` and `GLBs`.
    """
    n = len(lattice)
    # Here candidate_function[c] = glb(fn_1[c], fn_2[c], ..., fn_n[c]), for fn_i in functions
    candidate_function = [glb(i) for i in zip(*functions)]
    # One sub-list of good pairs for each element in the lattice
    good_pairs = [[] for i in range(n)]
    # All conflicting tuples form all elements in the lattice
    conflict_tuples = set()
    # Calculate all initial conflicts in the candidate solution
    for pair in combinations(range(n), 2):
        w = lub(pair)
        if check_fn_with_pair(candidate_function, pair):
            good_pairs[w].append(pair)
        else:
            conflict_tuples.add((pair, w))

    while len(conflict_tuples) != 0:
        (u, v), w = conflict_tuples.pop()
        candidate_function[w] = lub((candidate_function[u], candidate_function[v]))
        
        for x, y in good_pairs[w]:
            candidate_function[x] = glb((candidate_function[x], candidate_function[w]))
            candidate_function[y] = glb((candidate_function[y], candidate_function[w]))
            
            conflict_tuples = conflict_tuples | set(product(good_pairs[x], (x,)))
            conflict_tuples = conflict_tuples | set(product(good_pairs[y], (y,)))

            good_pairs[x].clear()
            good_pairs[y].clear()
        
        conflict_tuples = conflict_tuples | set(product(good_pairs[w], (w,)))
        good_pairs[w] = [(u, v)]

        # Filter non-conflits from the set
        for pair, z in conflict_tuples.copy():
            if check_fn_with_pair(candidate_function, pair):
                good_pairs[z].append(pair)
                conflict_tuples.remove((pair, z))
    
    return candidate_function

def delta_foo(lattice, functions):
    """ Calculates Delta using Greatest Lower Bounds and then fixes the
        resulting function until its a valid space-function.
        This function makes implicit use of the globals `LUBs` and `GLBs`.
    """
    n = len(lattice)
    # Here candidate_function[c] = glb(fn_1[c], fn_2[c], ..., fn_n[c]), for fn_i in functions
    candidate_function = [glb(i) for i in zip(*functions)]
    # One sub-list of good pairs for each element in the lattice
    good_pairs = [[] for i in range(n)]
    # All conflicting tuples form all elements in the lattice
    conflict_tuples = set()
    # All pairs of elements in the lattice that lost support
    falling_pairs = set()
    # Calculate all initial conflicts in the candidate solution
    for pair in combinations(range(n), 2):
        w = lub(pair)
        if check_fn_with_pair(candidate_function, pair):
            good_pairs[w].append(pair)
        else:
            conflict_tuples.add((pair, w))

    while len(conflict_tuples) != 0:
        (u, v), w = conflict_tuples.pop()
        candidate_function[w] = lub((candidate_function[u], candidate_function[v]))
        falling_pairs = falling_pairs | set(good_pairs[w])
        good_pairs[w].clear()

        while len(falling_pairs) != 0:
            x, y = falling_pairs.pop()
            z = lub((x, y))
            
            old_mapped_x = candidate_function[x]
            old_mapped_y = candidate_function[y]

            candidate_function[x], candidate_function[y] = (glb((candidate_function[x], candidate_function[z])), glb((candidate_function[y], candidate_function[z])))
            
            if old_mapped_x != candidate_function[x]:
                falling_pairs = falling_pairs | set(good_pairs[x])
                good_pairs[x].clear()
            
            if old_mapped_y != candidate_function[y]:
                falling_pairs = falling_pairs | set(good_pairs[y])
                good_pairs[y].clear()
            
            if check_fn_with_pair(candidate_function, (x, y)):
                good_pairs[z].append((x, y))
            else:
                conflict_tuples.add(((x, y), z))
        good_pairs[w].append((u, v))
    return candidate_function

def delta_n(lattice, space_functions, functions):
    """ Calculate Delta for a set of `functions` over a `lattice`
        and a all possible `space_functions` for that lattice.
    """
    valid_functions = (fn for fn in space_functions if all(map(lambda elem: leq_fn(lattice, fn, elem), functions)))
    return list(max_fn(lattice, (fn for fn in valid_functions)))

# ######################################

# ######################################
# Utility functions for lattice creation
# ######################################

def lattice_m(n):
    """Generate an M lattice with `n` elements"""
    l = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        l[i][0] = 1
        l[i][i] = 1
        l[n-1][i] = 1

    return l

def lattice_n5():
    """Generate the lattice n5"""
    n = 5
    l = [[0 for i in range(n)] for j in range(n)]
    for i in range(n):
        l[i][0] = 1
        l[i][i] = 1
        l[n-1][i] = 1
    l[3][1] = 1

    return l

def lattice_kite():
    """Generate a lattice with 4 elements"""
    return lattice_m(4)

def lattice_square():
    """Generate a square lattice with 9 elements"""
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
    """Generate a power of 3 lattice with 8 elements"""
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
    """
    Prints a table with the values for each pair a b in a lattice with `n` nodes showing
    `a`, `b`, `lub(fn[a], fn[b])`, `fn[lub(a, b)]` for testing porpuses.

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

class TestResults:
    """Utility class to group and show algorithms' results"""
    
    def __init__(self, name):
        # The name of the algorith tested
        self.name = name
        # Errors of the form (space_functions, bad_result, expected_result)
        self.errors = [] 
        self.avg_time = 0.0
        self.min_time = float("inf")
        self.max_time = 0.0
    
    def update_times(self, new_time, max_iterations):
        """ Updates the internal time counter with the `new_time`.
            Keeping track of the average (over `max_iterations`),
            minimun, and maximun time of the TestResult associated
            algorithm.
        """
        self.avg_time += new_time / max_iterations
        self.min_time = min(self.min_time, new_time)
        self.max_time = max(self.max_time, new_time)
    
    def print_times(self):
        """ Prints the current average, minimun and maximun time.
        """
        print("\nTimes for", self.name)
        print("-- Average:", self.avg_time)
        print("-- Min:    ", self.min_time)
        print("-- Max:    ", self.max_time)
    
    def print_status(self):
        """ Prints if the TestResults has failed any of the tests
            or passed all tests.
        """
        if len(self.errors) == 0:
            print("[i] All", self.name, "\x1b[32mpassed\x1b[0m")
        else:
            print("[w] Some", self.name, "\x1b[31mFAILED\x1b[0m")
    
    def print_errors(self):
        """ If the TestResult presents any failed test cases, then
            it prints each failed test case, showing the space functions
            used, the expected (correct) and the actual (reported) results.
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
    AST = 0
    AST_IMPLY = 1
    AST2 = 2
    FOO_OLD = 3
    FOO = 4
    AST_PART = 5

def run_test_case(fn, lattice, test_functions):
    """ Runs `fn` with a `lattice` and an iterable of `test_functions`.
        Returns the execution_time and the result.
    """
    fn_time = time()
    result = fn(lattice, test_functions)
    fn_time = time() - fn_time
    return fn_time, result

def run(lattice, verbose = False, test_functions = None, fns_file = None, save_functions = False):
    """ Run all the delta algorithms against a `lattice`
        and test_functions, either random or explicit, and
        present the results.
    """
    # Used to calculate elapsed time
    start_time = time()
    # Number of iterations to test (test-cases)
    n = 1

    preproc_time = time()
    # Matrix of least upper bounds and greatest lower bounds
    global LUBs, GLBs, IMPLs
    LUBs = calculate_lubs(lattice)
    GLBs = calculate_glbs(lattice)
    IMPLs = calculate_implications(lattice)
    preproc_time = time() - preproc_time

    func_gen_time = time()

    if fns_file is None:
        # Calculate space functions.
        space_functions = generate_functions(len(lattice))
    else:
        # Read space functions from a file
        try:
            # TODO: Make the file location relative to the script
            # TODO: Validate that the space functions from the file correspond to the lattice
            with open(fns_file) as f:
                space_functions = eval(f.read())
                print("[i] Reading space functions from file")
                save_functions = False # Nothing changed no need to overwrite
        except IOError:
            print("[w] File not found, generating space functions...")
            # Calculate space functions.
            space_functions = generate_functions(len(lattice))

    validation_time = time()
    # By default use from 1 to `n`, inclusive, random test_functions
    # if none are provided, or check the test_functions
    # provided by the user.
    # Note: It was "temporarily" changed to signify the number of test-cases.
    if test_functions is None:
        n = 1000
    else:
        valid_input = True
        for fn in test_functions:
            if fn not in space_functions:
                print("[e] Invalid Test Function found:", repr(fn))
                valid_input = False
        if not valid_input:
            print ("[e] Aborting execution!")
            return
    
    validation_time = time() - validation_time
    func_gen_time = time() - func_gen_time - validation_time

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
    delta_results = [
        TestResults("Delta*"),
        TestResults("Delta*(IMPLYs)"),
        TestResults("Delta*(may fail)"),
        TestResults("Delta_foo(old)"),
        TestResults("Delta_foo"),
        TestResults("Delta*(PART)")
    ]

    for i in range(1, n+1):
        # i is the number of test-functions per test-case
        # TODO: Remove hard-coded value.
        i = 4
        # Get some space functions at random or use given ones
        if test_functions is None:
            sample_functions = random.sample(range(len(space_functions)), i)
            sample_functions = [space_functions[j] for j in sample_functions]
        else:
            sample_functions = [fn for fn in test_functions if fn in space_functions]
        
        if verbose:
            print("__________Space Functions__________\n")
            for fn in sample_functions:
                print(fn)
            print("___________________________________\n")

        
        fn_time, delta_ast_result = run_test_case(delta_ast, lattice, sample_functions)
        delta_results[Delta.AST.value].update_times(fn_time, n)
        if verbose:
            print("Delta*:          ", repr(delta_ast_result))
            print("-- Time:", fn_time, "\n")
        
        fn_time, delta_ast_imply_result = run_test_case(delta_ast_imply, lattice, sample_functions)
        delta_results[Delta.AST_IMPLY.value].update_times(fn_time, n)
        if verbose:
            print("Delta*(IMPLYs):  ", repr(delta_ast_imply_result))
            print("-- Time:   ", fn_time)
            print("-- Without preprocessing:", fn_time + preproc_time, "\n")

        fn_time, delta_ast2_result = run_test_case(delta_ast2, lattice, sample_functions)
        delta_results[Delta.AST2.value].update_times(fn_time, n)
        if verbose:
            print("Delta*(may fail):", repr(delta_ast2_result))
            print("-- Time:", fn_time, "\n")

        fn_time, delta_foo_old_result = run_test_case(delta_foo_old, lattice, sample_functions)
        delta_results[Delta.FOO_OLD.value].update_times(fn_time, n)
        if verbose:
            print("Delta_foo(old):  ", repr(delta_foo_old_result))
            print("-- Time:", fn_time, "\n")

        fn_time, delta_foo_result = run_test_case(delta_foo, lattice, sample_functions)
        delta_results[Delta.FOO.value].update_times(fn_time, n)
        if verbose:
            print("Delta_foo:       ", repr(delta_foo_result))
            print("-- Time:", fn_time, "\n")

            print("Delta_0:         ", [glb(i) for i in zip(*sample_functions)])

        fn_time, delta_ast_part_result = run_test_case(delta_ast_partition, lattice, sample_functions)
        delta_results[Delta.AST_PART.value].update_times(fn_time, n)
        if verbose:
            print("Delta*(PART):    ", repr(delta_ast_part_result))
            print("-- Time:", fn_time, "\n")

        fn_time = time()
        delta_max_result = delta_n(lattice, space_functions, sample_functions)
        fn_time = time() - fn_time
        if verbose:
            print("Delta:           ", repr(delta_max_result))
            print("-- With preprocessing:   ", fn_time)
            print("-- Without preprocessing:", fn_time + func_gen_time)
            print("----------------------------------------------------------------\n")

        # If any of the algorithms failed to compute delta, add the error and context.
        if delta_ast_result != delta_max_result:
            delta_results[Delta.AST.value].errors.append((sample_functions, delta_ast_result, delta_max_result))
        if delta_ast_imply_result is not None and delta_ast_imply_result != delta_max_result:
            delta_results[Delta.AST_IMPLY.value].errors.append((sample_functions, delta_ast_imply_result, delta_max_result))
        if delta_ast2_result != delta_max_result:
            delta_results[Delta.AST2.value].errors.append((sample_functions, delta_ast2_result, delta_max_result))
        if delta_foo_old_result != delta_max_result:
            delta_results[Delta.FOO_OLD.value].errors.append((sample_functions, delta_foo_old_result, delta_max_result))
        if delta_foo_result != delta_max_result:
            delta_results[Delta.FOO.value].errors.append((sample_functions, delta_foo_result, delta_max_result))
        if delta_ast_part_result != delta_max_result:
            delta_results[Delta.AST_PART.value].errors.append((sample_functions, delta_ast_part_result, delta_max_result))

    print("Number of iterations:", n)

    # Show error test cases
    for result in delta_results:
        result.print_errors()

    # Show summary of results
    for result in delta_results:
        result.print_status()

    for result in delta_results:
        result.print_times()

    elapsed_time = time() - start_time
    print("\nTotal time:", elapsed_time)

# Call with one argument to generate random test cases
# run(lattice, verbose = True, test_functions = None, fns_file = None, save_functions = False)
# 
# verbose:
#   - If True, show each test case and the individual results,
#     and a summary of failing test cases
#   - If False (default), only show failing test cases (faster)
# test_functions:
#     The list of space-functions to test against. By default
#     if no test_functions are provided, then generate random test_functions.
# fns_file:
#     String with the file path to read and/or write calculated
#     space-functions for `lattice`. 
# save_functions:
#   - If True, write the generated space-function to `fns_file` except when the
#     space-functions were read from the same file.
#   - If False (default), do not write the generated space-functions anywhere.

run(lattice_square(), fns_file="sf_square.in")

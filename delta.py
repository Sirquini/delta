import random
import os
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
    """ Generate a list of space functions, based on a lattice
        of `n` elements.
    """
    test_pairs = list(combinations(range(n), 2))
    return [(0,) + fn for fn in product(range(n), repeat=n-1) if is_distributive(fn, test_pairs)]

def is_distributive(fn, test_pairs):
    """ Checks if a given function `fn` is distributive.
    """
    fn = (0,) + fn
    return all(fn[lub([t0, t1])] == lub([fn[t0], fn[t1]]) for t0, t1 in test_pairs)

def is_lattice(lattice):
    """ Returns True if `lattice` is a valid lattice.
    """
    N = len(lattice)
	# List of pairs, leaving them as an iterator allows us
    # to short-circuit and not generate all the pairs if
    # the lattice is not valid.
    lt = combinations(range(N),2)
    return all(lub(pair) != -1 for pair in lt)

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
    n = len(lattice)
    fn_num = len(functions)
    if fn_num == 1:
        return functions[0][c]
    else:
        mid_point = fn_num // 2
        return glb(lub((partition_helper(lattice, functions[:mid_point], a), partition_helper(lattice, functions[mid_point:], imply(a, c)))) for a in range(n) if lattice[c][a] == 1)

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
    good_pairs = [[] for _ in range(n)]
    # All conflicting tuples form all elements in the lattice
    conflict_tuples = set()
    # All pairs of elements in the lattice that lost support
    falling_pairs = set()
    # Calculate all initial conflicts in the candidate solution
    for pair in combinations(range(n), 2):
        w = lub(pair)
        u, v = pair
        if check_fn_with_pair(candidate_function, pair):
            good_pairs[w].append(pair)
        elif lattice[lub((candidate_function[u], candidate_function[v]))][candidate_function[w]] != 1:
            conflict_tuples.add((pair, w))
        else:
            falling_pairs.add(pair)

    while len(conflict_tuples) != 0:
        (u, v), w = conflict_tuples.pop()
        candidate_function[w] = lub((candidate_function[u], candidate_function[v]))
        falling_pairs = falling_pairs | set(good_pairs[w])
        good_pairs[w] = [(u,v)]

        while len(falling_pairs) != 0:
            x, y = falling_pairs.pop()
            z = lub((x, y))
            
            if candidate_function[x] != glb((candidate_function[x], candidate_function[z])):
                falling_pairs = falling_pairs | set(good_pairs[x])
                good_pairs[x].clear()
            
            if candidate_function[y] != glb((candidate_function[y], candidate_function[z])):
                falling_pairs = falling_pairs | set(good_pairs[y])
                good_pairs[y].clear()

            candidate_function[x] = glb((candidate_function[x], candidate_function[z]))
            candidate_function[y] = glb((candidate_function[y], candidate_function[z]))

            if lub((candidate_function[x], candidate_function[y])) == candidate_function[z]:
                good_pairs[z].append((x, y))
            else:
                conflict_tuples.add(((x, y), z))
        # # Check that all good_pairs remain good. Sanity Check. Performance penalty.
        # for sublist in good_pairs:
        #     for bad_pair in filter(lambda p: LUBs[candidate_function[p[0]]][candidate_function[p[1]]] != candidate_function[LUBs[p[0]][p[1]]], sublist):
        #         conflict_tuples.add((bad_pair, LUBs[bad_pair[0]][bad_pair[1]]))
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

# ######################################
#  Utility functions for handling files
# ######################################

def relative_path(path, *paths):
    """ Returns the absolute path for a file relative to the script
    """
    dirname = os.path.dirname(__file__)
    return os.path.join(dirname, path, *paths)

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
    FOO = 1
    AST_PART = 2
    N = 3

def run_test_case(fn, lattice, test_functions):
    """ Runs `fn` with a `lattice` and an iterable of `test_functions`.
        Returns the execution_time and the result.
    """
    fn_time = time()
    result = fn(lattice, test_functions)
    fn_time = time() - fn_time
    return fn_time, result

def run(lattice, verbose = False, test_functions = None, n_tests = 100, n_functions = 2, fns_file = None, save_functions = False):
    """ Run all the delta algorithms against a `lattice` and `test_functions`,
        either random or explicit, and present the results.

        INPUTS:
        - lattice:
            Matrix representing The lattice to test against.
        - verbose:
            - If True, show each test case and the individual results,
            and a summary of failing test cases
            - If False (default), only show failing test cases (faster)
        - test_functions:
            The list of space-functions to test against. By default, if no 
            `test_functions` are provided, generates random test_functions.
        - n_tests:
            If `test_functions` is None, indicates the number of test to run
            with `n_functions` random test functions.
        - n_functions:
            If `test_funtions` is None, indicates how many random test
            functions to use per test-case.
        - fns_file:
            String with the file path to read and/or write calculated
            space-functions for `lattice`. 
        - save_functions:
            - If True, write the generated space-function to `fns_file` except
            when the space-functions were read from the same file.
            - If False (default), do not write the generated space-functions
            anywhere.
        
        OUPUTS:
            Returns a dictionary with the acumulated results, including timings
            and failures (if any).
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

    # Valitdate the input lattice
    if not is_lattice(lattice):
        print("[e] Invalid lattice, aborting execution!")
        return {}

    func_gen_time = time()
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
    
    validation_time = time()
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
    delta_results = {
        Delta.AST: TestResults("Delta*"),
        Delta.FOO: TestResults("Delta_foo"),
        Delta.AST_PART: TestResults("Delta*(PART)"),
        Delta.N: TestResults("Delta_n(Preprocessed)")
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

        fn_time, delta_ast_result = run_test_case(delta_ast, lattice, sample_functions)
        delta_results[Delta.AST].update_times(fn_time, n)
        if verbose:
            print("Delta*:          ", repr(delta_ast_result))
            print("-- Time:", fn_time, "\n")

        fn_time, delta_foo_result = run_test_case(delta_foo, lattice, sample_functions)
        delta_results[Delta.FOO].update_times(fn_time, n)
        if verbose:
            print("Delta_foo:       ", repr(delta_foo_result))
            print("-- Time:", fn_time, "\n")
            print("Delta_0:         ", [glb(i) for i in zip(*sample_functions)], "\n")

        fn_time, delta_ast_part_result = run_test_case(delta_ast_partition, lattice, sample_functions)
        delta_results[Delta.AST_PART].update_times(fn_time, n)
        if verbose:
            print("Delta*(PART):    ", repr(delta_ast_part_result))
            print("-- Time:", fn_time, "\n")

        fn_time = time()
        delta_max_result = delta_n(lattice, space_functions, sample_functions)
        fn_time = time() - fn_time
        delta_results[Delta.N].update_times(fn_time, n)
        if verbose:
            print("Delta:           ", repr(delta_max_result))
            print("-- With preprocessing:   ", fn_time)
            print("-- Without preprocessing:", fn_time + func_gen_time)
            print("----------------------------------------------------------------\n")

        # If any of the algorithms failed to compute delta, add the error and context.
        if delta_ast_result != delta_max_result:
            delta_results[Delta.AST].errors.append((sample_functions, delta_ast_result, delta_max_result))
        if delta_foo_result != delta_max_result:
            delta_results[Delta.FOO].errors.append((sample_functions, delta_foo_result, delta_max_result))
        if delta_ast_part_result != delta_max_result:
            delta_results[Delta.AST_PART].errors.append((sample_functions, delta_ast_part_result, delta_max_result))

    print("Number of iterations:", n)

    # Show error test cases
    for result in delta_results.values():
        result.print_errors()

    # Show summary of results
    for result in delta_results.values():
        result.print_status()

    for result in delta_results.values():
        result.print_times()

    elapsed_time = time() - start_time
    print("\nTotal time:", elapsed_time)
    return {"result": delta_results.values(), "total": elapsed_time, "sf": len(space_functions), "sf_gen_time": func_gen_time, "preproc": preproc_time}

def write_test_results_csv(name, results):
    """ Write the test results in a file for later usage
        in a CSV file.
    """
    import csv
    # Collect headers
    if results:
        headers = ["Lattice", "Nodes", "Space Functions", "Test Functions", "Elapsed Time", "Preprocessing Time"]
        for delta_result in results[0]["result"]:
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
                    row.append(delta_result.avg_time)
                    row.append(delta_result.min_time)
                    row.append(delta_result.max_time)
                    row.append(len(delta_result.errors))
                writer.writerow(row)
        print("[i] Generated CSV file {}".format(name))

def run_full_tests():
    from lattice import process_file
    from datetime import datetime
    # Read the lattice to test from
    lattices = process_file("distributive_lattices.py")
    results = []
    start_time = datetime.now().strftime("%Y-%m-%d-%H%M")
    for key, lattice in lattices.items():
        print("* Using lattice `{}` ({} nodes)".format(key, len(lattice)))
        for i in range(2, 6):
            print("Test Functions:", i)
            result = run(lattice, n_functions=i, fns_file=relative_path("generated", "sf_{}.in".format(key)))
            result["lattice"] = key
            result["nodes"] = len(lattice)
            result["functions"] = i
            results.append(result)
        print("================================================================\n")
    write_test_results_csv("results-{}.csv".format(start_time), results)

def run_square():
    run(lattice_square(), n_tests=1000, n_functions=4, fns_file=relative_path("generated","sf_square.in"))

def run_failling_foo():
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
    from lattice import random_lattice, lattice_from_covers
    covers = random_lattice(8, 0.95)
    print("* Using lattice: ", covers)
    lattice = lattice_from_covers(covers)
    run(lattice)

if __name__ == "__main__":
    # run_full_tests()
    # run_square()
    # run_random()
    run_failling_foo()

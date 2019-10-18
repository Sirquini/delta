import os
import sys
import random
import json
from time import perf_counter
from itertools import product, combinations, permutations

from lattice import Lattice, powerset_lattice, delta_foo, delta_ast_partition

# #######################################
# Experiments mainly over lattices,
# sapce-functions, Delta implementations,
# Lattice implementation, etc.
# #######################################

# #######################################
# Utility fumctions for file manipulation
# #######################################

def from_rel_path(path, *paths):
    """Returns the absolute path of a file relative to the script."""
    dirname = os.path.dirname(__file__)
    return os.path.join(dirname, path, *paths)

def eprint(*args, **kwargs):
    """Similar to `print` but uses `sys.stderr` instead of `sys.stdin`
    
    Also flushes the stream.
    """
    print(*args, file=sys.stderr, flush=True, **kwargs)

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
        with open(from_rel_path("results", name), 'w', newline='') as f:
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

# #######################################

# #######################################
# Utility functions for running tests.
# #######################################

from delta import TestResults, Delta, run_test_case, delta_n, random_space_function, all_space_functions
from generation import progress_bar

# #######################################

# #######################################
# Main Experiments: `run_` prefix
# #######################################

def run(l, verbose = False, test_functions = None, n_tests = 100, n_functions = 2, fns_file = None, save_functions = False):
    """Runs all the delta algorithms against a `lattice` and `test_functions`,
    either random or explicit, and presents the results.

    Args:
        l:
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
    # Avoid creating globals, Use the Lattice class.
    lattice = Lattice(l)
    lattice.impls # This value is lazy, we must call it to precompute it.
    preproc_time = perf_counter() - preproc_time

    # Valitdate the input lattice
    if not lattice.is_lattice():
        print("[E] Invalid lattice, aborting execution!")
        return {}

    func_gen_time = perf_counter()
    if fns_file is None:
        # Calculate space functions.
        space_functions = lattice.space_functions
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
            space_functions = lattice.space_functions
    
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
            print("[E] Aborting execution!")
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
        Delta.FOO: TestResults("DeltaGen+"),
        Delta.AST_LATEST: TestResults("DeltaPart3+"),
        Delta.N: TestResults("Brute-force")
    }

    for _ in range(n):
        # Get some space functions at random or use given ones
        if test_functions is None:
            sample_functions = random.sample(range(len(space_functions)), n_functions)
            sample_functions = [space_functions[j] for j in sample_functions]
        else:
            sample_functions = test_functions
        
        if verbose:
            print("__________Space Functions__________\n")
            for fn in sample_functions:
                print(fn)
            print("___________________________________\n")

        fn_time, delta_foo_result = run_test_case(delta_foo, lattice, sample_functions)
        delta_results[Delta.FOO].update_times(fn_time, n)
        if verbose:
            print("{}: {}".format(delta_results[Delta.FOO].name, repr(delta_foo_result)))
            print("-- Time: {}\n".format(fn_time))

        fn_time, delta_ast_part_result = run_test_case(delta_ast_partition, lattice, sample_functions)
        delta_results[Delta.AST_LATEST].update_times(fn_time, n)
        if verbose:
            print("{}: {}".format(delta_results[Delta.AST_LATEST].name, repr(delta_ast_part_result)))
            print("-- Time: {}\n".format(fn_time))

        fn_time = perf_counter()
        delta_max_result = delta_n(lattice.lattice, space_functions, sample_functions)
        fn_time = perf_counter() - fn_time
        delta_results[Delta.N].update_times(fn_time, n)
        if verbose:
            print("{}: {}".format(delta_results[Delta.N].name, repr(delta_max_result)))
            print("-- With preprocessing:   ", fn_time)
            print("-- Without preprocessing:", fn_time + func_gen_time)
            print("----------------------------------------------------------------\n")

        # If any of the algorithms failed to compute delta, add the error and context.
        if delta_foo_result != delta_max_result:
            delta_results[Delta.FOO].errors.append((sample_functions, delta_foo_result, delta_max_result))
        if delta_ast_part_result != delta_max_result:
            delta_results[Delta.AST_LATEST].errors.append((sample_functions, delta_ast_part_result, delta_max_result))

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

def run_full_tests(size_limit=0):
    """Run test with a pre-generated list of random lattices, with 
    pre-calculated space functions.

    Args:
      size_limit: Value, different than 0, it will limit the lattices from the
      pre-generated list to only those with size equal to `size_limit`
    """
    from lattice import process_file
    from datetime import datetime
    # Read the lattice to test from
    lattices = process_file("distributive_lattices.py")
    if size_limit > 0:
        lattices = {key: value for (key, value) in lattices.items() if len(value) == size_limit}
    results = []
    start_time = datetime.now().strftime("%Y-%m-%d-%H%M")
    count = 0
    for key, lattice in lattices.items():
        print("* Using lattice `{}` ({} nodes)".format(key, len(lattice)))
        count += 1
        progress_bar(count, len(lattices), 50)
        for i in range(4, 33, 4):
            print("Test Functions:", i)
            result = run(lattice, n_functions=i, fns_file=from_rel_path("generated", "sf_{}.in".format(key)))
            result["lattice"] = key
            result["nodes"] = len(lattice)
            result["functions"] = i
            results.append(result)
        print("================================================================\n")
    eprint(" Done.")
    write_test_results_csv("results-{}.csv".format(start_time), results)

def run_powerset(exponent = 10, verbose = False, test_functions = None, n_tests = 10, n_functions = 2, brute_force=False):
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
        brute_force:
            - If `True`, generates all the possible space functions and runs
            the brute-force delta calculation to compare. The default is 
            `False` since we usually don't want to use brute-force with big
            lattices.

    Returns:
        A dictionary with the acumulated results, including timings and
        failures (if any).
    """
    # Used to calculate the elapsed time
    start_time = perf_counter()
    # Number of iterations to test (test-cases)
    n = 1
    # The actual number of elements in the lattice = 2^n
    elements = 2**exponent

    eprint("[i] Generating powerset lattice with {} nodes ".format(elements), end='')

    # Generate the powerset lattice and measure times
    preproc_time = perf_counter()
    lattice = Lattice(powerset_lattice(exponent))
    lattice.impls # This value is lazy, we must call it to precompute it.
    preproc_time = perf_counter() - preproc_time

    eprint(".", end='')

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

    pb_max = 2*n # used for the progress_bar
    pb_c = 0 # current test index
    func_gen_time = 0
    if brute_force:
        pb_max = 3*n
        func_gen_time = perf_counter()
        if exponent == 5 and POWERSET_32_SF is not None:
            space_functions = POWERSET_32_SF
        else:
            space_functions = all_space_functions(lattice)
        func_gen_time = perf_counter() - func_gen_time
        eprint(".")
        print("Space functions:", len(space_functions))
        print("Space functions preprocessing time:   ", func_gen_time)
    else:
        eprint("")
    print("LUBs, GLBs, IMPLYs preprocessing time:", preproc_time)
    print("________________________________________________________________\n")

    # Used for showing the aggregate results at the end
    delta_results = {
        Delta.FOO: TestResults("DeltaGen+"),
        Delta.AST_LATEST: TestResults("DeltaPart3+")
    }

    if brute_force:
        delta_results[Delta.N] = TestResults("Brute-force")

    deltas_are_equal = TestResults("assert_equal(DeltaGen+, DeltaPart3+)", True)

    progress_bar(pb_c, pb_max, 50, "DeltaGen+")
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

        fn_time, delta_foo_result = run_test_case(delta_foo, lattice, sample_functions)
        delta_results[Delta.FOO].update_times(fn_time, n)
        if verbose:
            print("{}: {}".format(delta_results[Delta.FOO].name, repr(delta_foo_result)))
            print("-- Time: {}\n".format(fn_time))
        pb_c += 1
        progress_bar(pb_c, pb_max, 50, "DeltaP3+")

        fn_time, delta_ast_partition_result = run_test_case(delta_ast_partition, lattice, sample_functions)
        delta_results[Delta.AST_LATEST].update_times(fn_time, n)
        if verbose:
            print("{}: {}".format(delta_results[Delta.AST_LATEST].name, repr(delta_ast_partition_result)))
            print("-- Time: {}\n".format(fn_time))
        pb_c += 1
        progress_bar(pb_c, pb_max, 50, "BruteForce")

        if brute_force:
            fn_time = perf_counter()
            delta_max_result = delta_n(lattice.lattice, space_functions, sample_functions)
            fn_time = perf_counter() - fn_time
            delta_results[Delta.N].update_times(fn_time, n)
            if verbose:
                print("{}: {}".format(delta_results[Delta.N].name, repr(delta_max_result)))
                print("-- With preprocessing:   ", fn_time)
                print("-- Without preprocessing:", fn_time + func_gen_time)
                print("----------------------------------------------------------------\n")
            pb_c += 1
            progress_bar(pb_c, pb_max, 50, "DeltaGen+")

        if delta_foo_result != delta_ast_partition_result:
            deltas_are_equal.errors.append((sample_functions, delta_foo_result, delta_ast_partition_result))
    
    eprint(" Done.")
    print("Number of iterations:", n)

    deltas_are_equal.print_errors()
    deltas_are_equal.print_status()

    for result in delta_results.values():
        result.print_times()

    delta_results[Delta.OTHER] = deltas_are_equal

    elapsed_time = perf_counter() - start_time
    print("\nTotal time:", elapsed_time)
    return {"result": delta_results.values(), "total": elapsed_time, "sf": 0, "sf_gen_time": func_gen_time, "preproc": preproc_time}

def run_full_powerset_tests():
    """Equivalent to `run_full_tests()` but is used to time algorithms against
    powerset lattices.
    """
    from datetime import datetime
    results = []
    start_time = datetime.now().strftime("%Y-%m-%d-%H%M")
    for exponent in range(5, 6):
        nodes = 2**exponent
        print("* Using lattice `Powerset_{}` ({} nodes)".format(exponent, nodes))
        for i in [4, 8, 12,16, 24, 28, 32]: # TODO: Restore this to [4, 8, 12, 16]. Use only [4] for 4 and 5. Since Delta+ is O(mn^m), where m=n_functions 
            print("\nTest Functions:", i)
            result = run_powerset(exponent, n_tests=10, n_functions=i, brute_force=nodes<=32)
            result["lattice"] = "Powerset_{}".format(exponent)
            result["nodes"] = nodes
            result["functions"] = i
            results.append(result)
        print("================================================================\n")
    eprint(" Done.")
    write_test_results_csv("results-{}.csv".format(start_time), results)

# #######################################

def run_arbitrary_lattices(sizes):
    """Runs and times DeltaGen+ against arbitrary lattices and saves the
    results to a CSV file.

    Agrs:
      sizes: A list of lattice sizes to run.
    
    IO:
      Generates a CSV file with the results of the run, the name of the file
      follows the format `results-YEAR-MONTH-DAY-TIME.csv`
    """
    from datetime import datetime
    from time import perf_counter

    from experimental import NxLattice
    from lattice import random_lattice, lattice_from_covers

    results = []
    start_time = datetime.now().strftime("%Y-%m-%d-%H%M")
    

    progress_bar(0, 1, 50, "Generating")
    # Generate the random lattices, using experimental parallel computation.
    lattices = [[NxLattice(lattice_from_covers(random_lattice(size, 0.95)), size>7) for _ in range(10)] for size in sizes]
    progress_bar(1, 1, 50, "Generating")

    # used for counting
    n_lattices = len(lattices) * 10
    count = 0

    # Calculate space functions and IMPLYs
    progress_bar(count, n_lattices, 50, "Preprocc")
    for batch in lattices:
        for lattice in batch:
            lattice.impls
            lattice.space_functions
            count += 1
            progress_bar(count, n_lattices, 50, "Preprocc")

    eprint(" Done.")
    write_test_results_csv("results-{}.csv".format(start_time), results)

if __name__ == "__main__":
    # This is only temporary, since the number of space functions for powerset
    # with 32 nodes is 32^5
    # global POWERSET_32_SF
    POWERSET_32_SF = all_space_functions(Lattice(powerset_lattice(5)))
    # POWERSET_32_SF = None
    # run_full_tests(size_limit=10)
    run_full_powerset_tests()

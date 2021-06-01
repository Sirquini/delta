import os
import sys
import random
from time import perf_counter
from itertools import combinations

from lattice import Lattice, powerset_lattice, delta_foo, delta_ast_partition, delta_partition

# #######################################
# Experiments mainly over lattices,
# sapce-functions, Delta implementations,
# Lattice implementations, etc.
# #######################################

# #######################################
# Utility functions for file manipulation
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

from delta import TestResults, Delta, lattice_square, run_test_case, delta_n, random_space_function, all_space_functions, delta_plus_jies 
from generation import progress_bar

class DeltaTest:
    def __init__(self, name, fn, skip_times=False):
        self.fn = fn
        self.test_results = TestResults(name, skip_times)

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
        Delta.AST_LATEST: TestResults("DeltaJIEs"),
        Delta.OTHER: TestResults("DeltaPlus"),
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
        
        fn_time, delta_plus_result = run_test_case(delta_plus_jies, lattice, sample_functions)
        delta_results[Delta.OTHER].update_times(fn_time, n)
        if verbose:
            print("{}: {}".format(delta_results[Delta.OTHER].name, repr(delta_plus_result)))
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
        if delta_plus_result != delta_max_result:
            delta_results[Delta.AST_LATEST].errors.append((sample_functions, delta_plus_result, delta_max_result))

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

def run_random():
    """Equivalent to calling `run()` with a random lattice as its parameter."""
    from lattice import random_lattice, lattice_from_covers
    covers = random_lattice(8, 0.95)
    print("* Using lattice:", covers)
    lattice = lattice_from_covers(covers)
    run(lattice)

def run_deltas(l: Lattice, deltas, verbose = False, test_functions = None, n_tests = 10, n_functions = 2, brute_force=False):
    """Runs all the `deltas` algorithms against a given lattice and
    `test_functions`, either random or explicit, similar to the run command,
    and present the results.

    Args:
        l:
            Matrix representing the lattice to test against.
        deltas:
            The list of DeltaTests functions to time and test.
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
    # NOTE: This function is NOT finished. It is meant to be a generic alternative
    # to using `run` and `run_powerset`.
    pass

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

    pb_max = 4*n # used for the progress_bar
    pb_c = 0 # current test index
    func_gen_time = 0
    if brute_force:
        pb_max += n

        func_gen_time = perf_counter()
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
        Delta.AST: TestResults("DeltaJIe"),
        Delta.PLUS: TestResults("DeltaPlus"),
        Delta.FOO: TestResults("DeltaGen+"),
        Delta.AST_LATEST: TestResults("DeltaPart3+")
    }

    if brute_force:
        delta_results[Delta.N] = TestResults("Brute-force")

    deltas_are_equal = TestResults("assert_equal(DeltaPlus, DeltaJIe)", True)

    progress_bar(pb_c, pb_max, 50, "DeltaPlus")
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

        fn_time, delta_assoc_result = run_test_case(delta_plus_jies, lattice, sample_functions)
        delta_results[Delta.PLUS].update_times(fn_time, n)
        if verbose:
            print("{}: {}".format(delta_results[Delta.PLUS].name, repr(delta_assoc_result)))
            print("-- Time: {}\n".format(fn_time))
        pb_c += 1
        progress_bar(pb_c, pb_max, 50, "DeltaJIe")

        fn_time, delta_ast_result = run_test_case(delta_partition, lattice, sample_functions)
        delta_results[Delta.AST].update_times(fn_time, n)
        if verbose:
            print("{}: {}".format(delta_results[Delta.AST].name, repr(delta_ast_result)))
            print("-- Time: {}\n".format(fn_time))
        pb_c += 1
        progress_bar(pb_c, pb_max, 50, "DeltaGen+")

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

        if brute_force:
            progress_bar(pb_c, pb_max, 50, "BruteForce")
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
        progress_bar(pb_c, pb_max, 50, "DeltaPlus")

        if delta_assoc_result != delta_ast_result:
            deltas_are_equal.errors.append((sample_functions, delta_assoc_result, delta_ast_result))
    
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
            result = run_powerset(exponent, n_tests=10, n_functions=i, brute_force=nodes<32)
            result["lattice"] = "Powerset_{}".format(exponent)
            result["nodes"] = nodes
            result["functions"] = i
            results.append(result)
        print("================================================================\n")
    eprint(" Done.")
    write_test_results_csv("results-{}.csv".format(start_time), results)

# #######################################

def run_arbitrary_lattices(sizes, fixed_lattice=False):
    """Runs and times DeltaGen+ against arbitrary lattices and saves the
    results to a CSV file.

    Agrs:
      sizes: A list of lattice sizes to run/number of space functions.
      fixed_lattice: 
      - If `True`, `size` becomes the number of space functions per lattice of
      size 10.
      - If `False` (default), `size` becomes the sizes of lattices with a fixed
      number of space functions (4).
    
    IO:
      Generates a CSV file with the results of the run, the name of the file
      follows the format `results-YEAR-MONTH-DAY-TIME.csv`
    """
    from datetime import datetime
    from time import perf_counter

    from experimental import NxLattice
    from lattice import random_lattice, lattice_from_covers

    elapsed_time = perf_counter()
    start_time = datetime.now().strftime("%Y-%m-%d-%H%M")

    progress_bar(0, 1, 50, "Generating")
    # Generate the random lattices, using experimental parallel computation.
    if fixed_lattice:
        lattices = [NxLattice(lattice_from_covers(random_lattice(10, 0.95)), True)]
    else:
        lattices = [NxLattice(lattice_from_covers(random_lattice(size, 0.95)), size>7) for size in sizes]
    progress_bar(1, 1, 50, "Generating")

    # used for counting
    n_lattices = len(lattices)
    count = 0

    # Calculate space functions and IMPLYs
    sf_gen_times = []
    preproc_times = []
    progress_bar(count, n_lattices, 50, "Preproc")
    for lattice in lattices:
        preproc_time = perf_counter()
        lattice.impls
        preproc_times.append(perf_counter() - preproc_time)
        sf_gen_time = perf_counter()
        lattice.space_functions
        sf_gen_times.append(perf_counter() - sf_gen_time)
        count += 1
        progress_bar(count, n_lattices, 50, "Preproc")

    # Run each set of lattices, take times and save the results
    if fixed_lattice:
        n_lattices = len(sizes)
        items = sizes
        lattice = lattices[0]
    else:
        items = lattices
    count = 0
    time_results = [] # ([(delta_gen_time, brute_force_time), ...], test_time)
    progress_bar(count, n_lattices * 100, 50, "Running")
    for item in items:
        batch_results = []
        test_time = perf_counter()
        if fixed_lattice:
            size = item
        else:
            lattice = item
            size = 4
        # 100 runs per lattice
        for _ in range(100):
            # Generate som random space functions
            test_functions = random.sample(lattice.space_functions, size)

            # Run and time each algorithm
            delta_gen_time, _ = run_test_case(delta_foo, lattice, test_functions)
            delta_n_time = perf_counter()
            delta_n(lattice.lattice, lattice.space_functions, test_functions)
            delta_n_time = perf_counter() - delta_n_time

            # Save the results.
            batch_results.append((delta_gen_time, delta_n_time))
        time_results.append((batch_results, perf_counter() - test_time))
        count += 1
        progress_bar(count, n_lattices * 100, 50, "Running")

    results = []
    count = 0
    progress_bar(count, n_lattices, 50, "Mapping")
    # Format the data and save it to a CSV file for later use.
    for i, lattice in enumerate(lattices):
        if fixed_lattice:
            items = sizes
        else:
            items = [4]
        for k, size in enumerate(items):
            if fixed_lattice:
                i = k
                sf_gen_time = sf_gen_times[0]
                preproc_time = preproc_times[0]
            else:
                sf_gen_time = sf_gen_times[i]
                preproc_time = preproc_times[i]
            test_gen = TestResults("DeltaGen+")
            test_bf = TestResults("Brute-force")

            test_gen.avg_time = sum(times[0] for times in time_results[i][0]) / len(time_results[i][0])
            test_bf.avg_time = sum(times[1] for times in time_results[i][0]) / len(time_results[i][0])

            test_gen.min_time = min(times[0] for times in time_results[i][0])
            test_bf.min_time = min(times[1] for times in time_results[i][0])

            test_gen.max_time = max(times[0] for times in time_results[i][0])
            test_bf.max_time = max(times[1] for times in time_results[i][0])
            
            results.append({
                "result": (test_gen, test_bf),
                "total": time_results[i][1],
                "sf": len(lattice.space_functions),
                "sf_gen_time": sf_gen_time,
                "preproc": preproc_time,
                "lattice": f"id_{i}_n_{len(lattice)}",
                "nodes": len(lattice),
                "functions": size,
            })

            count += 1
            progress_bar(count, n_lattices, 50, "Mapping")

    elapsed_time = perf_counter() - elapsed_time
    eprint(" Done.\nRan for {} seconds.".format(elapsed_time))
    write_test_results_csv("results-{}.csv".format(start_time), results)

def run_square():
    run(lattice_square(), n_tests=1000, n_functions=4, fns_file=from_rel_path("generated","sf_square.in"))

def run_failling_foo():
    """Test-cases in which delta_foo used to fail."""
    from lattice import process_file
    
    lattices = process_file("distributive_lattices.py")
    test_functions = [(0, 5, 3, 7, 3, 3, 7, 7),
        (0, 7, 7, 7, 0, 7, 7, 7),
        (0, 3, 4, 7, 6, 6, 7, 7)]

    run(lattices[4698136515449058355], test_functions=test_functions, fns_file=from_rel_path("generated", "sf_4698136515449058355.in"))

    test_functions = [(0, 0, 7, 7, 9, 9, 9, 9, 7, 9),
        (0, 3, 3, 3, 1, 3, 3, 3, 9, 9),
        (0, 7, 1, 7, 8, 8, 9, 9, 7, 9)]

    run(lattices[-286441500945297568], test_functions=test_functions, fns_file=from_rel_path("generated", "sf_-286441500945297568.in"))

    test_functions = [(0, 5, 6, 6, 8, 3, 8, 8, 8, 8),
        (0, 5, 4, 9, 9, 9, 9, 9, 9, 9),
        (0, 8, 5, 8, 8, 0, 5, 8, 8, 8),
        (0, 2, 8, 8, 8, 7, 8, 8, 8, 8),
        (0, 5, 4, 9, 9, 5, 9, 5, 9, 9)]

    run(lattices[8013726431884705816], test_functions=test_functions, fns_file=from_rel_path("generated", "sf_8013726431884705816.in"))
    
    test_functions = [(0, 8, 1, 8, 9, 5, 7, 8, 8, 9),
        (0, 4, 1, 4, 4, 9, 9, 9, 9, 9),
        (0, 6, 9, 9, 9, 3, 9, 8, 9, 9),
        (0, 7, 8, 8, 9, 1, 8, 7, 8, 9),
        (0, 1, 0, 1, 1, 9, 9, 9, 9, 9)]

    run(lattices[8013726431884705816], test_functions=test_functions, fns_file=from_rel_path("generated", "sf_8013726431884705816.in"))

    test_functions = [(0, 4, 8, 9, 9, 9, 9, 9, 9, 9),
        (0, 4, 7, 9, 9, 8, 8, 9, 9, 9),
        (0, 2, 4, 4, 9, 1, 4, 3, 4, 9),
        (0, 7, 4, 9, 9, 8, 9, 8, 9, 9),
        (0, 5, 2, 6, 8, 4, 4, 9, 9, 9)]

    run(lattices[8013726431884705816], test_functions=test_functions, fns_file=from_rel_path("generated", "sf_8013726431884705816.in"))

if __name__ == "__main__":
    # run_full_tests(size_limit=10)
    # run_random()
    # run_powerset(exponent=3, verbose=True, test_functions=[(0,4,5,6,7,7,7,7),(0,3,2,1,6,5,4,7)])
    run_full_powerset_tests()
    # run_arbitrary_lattices([4, 5, 6, 7, 8, 9, 10])
    # run_arbitrary_lattices([4, 8, 12, 16, 24, 28, 32], fixed_lattice=True)

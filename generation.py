import sys
import random
from itertools import combinations, product
from lattice import Lattice, covers_from_lattice

def progress_bar(value, end, length=20):
    percent = value * 100 // end
    bar = '=' * (percent * length // 100)
    spaces = ' ' * (length - len(bar))
    print("\rProgress: [{}] {}/{}".format(bar + spaces, value, end), end='', file=sys.stderr, flush=True)

def print_distribution(dist, length=20, extra=None):
    top = max(dist.values())
    for d in sorted(dist):
        value = dist[d]
        bar = '-' * (value * length // top)
        spaces = ' ' * (length - len(bar))
        extra_info = ''
        if extra is not None:
            for entry in extra[d]:
                extra_info += "({} {})".format(entry, extra[d][entry])
        print("{} | {}({}){}".format(spaces + bar, d, value, extra_info))

def binomial(n, k):
    """ Calculates the binomial coefficient using the multiplicative formula.
    """
    if k < 0 or k > n:
        return 0
    result = 1
    for i in range(1, min(k, n - k) + 1):
        result = result * (n + 1 - i) // i
    return result

def is_distributive(lattice):
    """ Test if the `lattice` is a distributive lattice.

    Args:
        lattice: A Lattice instance.
    """
    n = len(lattice)
    pairs = combinations(range(n), 2)
    return all(lattice.glbs[lattice.lubs[a][b]][x] == lattice.lubs[lattice.glbs[a][x]][lattice.glbs[b][x]]
        for a, b in pairs for x in range(n))


def generate_subsets(n, p):
    """ Generates up to 2^`n` subsets with probability `p` of keeping a subset.
    """
    result = set()
    population = tuple(range(n))
    for _ in range(2**n):
        if random.random() < p:
            element = random.sample(population, random.randrange(0, n))
            result.add(frozenset(element))
    return result

def add_missing_subsets(subsets):
    """ Adds missing subsets so that it is closed under Union and Intersection.
    """
    queue = combinations(subsets, 2)
    while True:
        new_elements = []
        for a, b in queue:
            union = a | b
            intersection = a & b
            if union not in subsets:
                subsets.add(union)
                new_elements.append(union)
            if intersection not in subsets:
                subsets.add(intersection)
                new_elements.append(intersection)
        if len(new_elements) == 0:
            break
        queue = product(new_elements, subsets)
    return subsets


def distributive_lattice(n, p):
    """ Generates a random distributive lattice.

        The generated lattice can be used with the `Lattice` class, or
        converted with the `covers_from_lattice` function.

        INPUTS:
        + `n`: Number of elements of the population used to create the lattice.
        + `p`: Probability, from 0 to 1, for keeping a node.
    """
    assert(n>0)
    # Generate a random set of sets and close it under union and intersection.
    subsets = list(add_missing_subsets(generate_subsets(n, p)))

    if len(subsets) == 0:
        subsets = [frozenset()]

    # Find and reposition the top and bottom element
    # This is required because the bottom is asumed to be the first element,
    # and the top to be the last, this is an implementation detail that
    # may be removed in a revised version of the lattice class.
    n = len(subsets) - 1

    top_pos = subsets.index(max(subsets, key=len))
    subsets[top_pos], subsets[n] = subsets[n], subsets[top_pos]
    bottom_pos = subsets.index(min(subsets, key=len))
    subsets[bottom_pos], subsets[0] = subsets[0], subsets[bottom_pos]

    # Conver the set into a lattice's matrix representation
    result = [[0 for _ in range(len(subsets))] for _ in range(len(subsets))]
    for i, row in enumerate(result):
        for j, _ in enumerate(row):
            if subsets[j].issubset(subsets[i]):
                result[i][j] = 1
    return result

def test_distribution(n, p, samples):
    distribution = {}
    print("Running with n={} p={} and {} samples".format(n, p, samples))
    progress_bar(0, samples, 50)
    for i in range(1, samples + 1):
        if i % 1000 == 0:
            progress_bar(i, samples, 50)
        size = len(distributive_lattice(n, p))
        distribution[size] = distribution.get(size, 0) + 1
    print(" Done.", file=sys.stderr, flush=True)
    print_distribution(distribution, 50)
    return distribution


def test_distributive_lattices(n, p, limit, samples):
    from experimental import NxLattice
    distribution = {}
    extra = {}
    skipped = 0
    invalid = 0
    exceptions = 0
    print("Running with n={} p={} and a size limit of {}".format(n, p, limit))
    for i in range(1, samples + 1):
        progress_bar(i, samples, 40)
        random_lattice = distributive_lattice(n, p)
        size = len(random_lattice)
        if size > limit:
            skipped += 1
        else:
            distribution[size] = distribution.get(size, 0) + 1
            if size > 7:
                random_lattice = NxLattice(random_lattice, parallel=True)
            else:
                random_lattice = NxLattice(random_lattice, parallel=False)
            if is_distributive(random_lattice):
                ln_sf = binomial(2 * size - 2, size - 1)
                sf = len(random_lattice.space_functions)
                if size not in extra:
                    extra[size] = {"min": sf, "max": sf, "Ln": ln_sf}
                else:
                    if extra[size]["min"] > sf:
                        extra[size]["min"] = sf
                    if extra[size]["max"] < sf:
                        extra[size]["max"] = sf
                if ln_sf < sf:
                    exceptions += 1
                    print("\n* Lattice {} : {}".format(len(random_lattice), covers_from_lattice(random_lattice)))
                    print("+ Found {} (vs {} for Ln)".format(sf, ln_sf))
            else:
                invalid += 1
    print(" Done.", file=sys.stderr, flush=True)
    print("\nFinished with {} tests, {} skipped, {} invalid lattices".format(samples - skipped - invalid, skipped, invalid))
    print("Lattices with more SF than Ln:", exceptions)
    print("Distribution of Lattices by size:")
    print_distribution(distribution, 40, extra)

if __name__ == "__main__":
    from delta import lattice_square, lattice_m
    from lattice import test_equality
    square = Lattice(lattice_square())
    m3 = Lattice(lattice_m(5))
    print("Running pre-checks")
    # Testing binomial implementation
    test_equality(20, binomial(6, 3), "binomial 6 3 is 20")
    # Testing is_distributive
    test_equality(True, is_distributive(square), "is_distributive square is True")
    test_equality(False, is_distributive(m3), "is_distributive m_3 is False")

    test_distributive_lattices(5, 0.2, 9, 100)
    # test_distribution(4, 0.3, 100000)

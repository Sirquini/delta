import random
from delta import random_space_function
from lattice import Lattice, lattice_from_covers, random_lattice, powerset_lattice, delta_foo
from generation import distributive_lattice

# ###########
# Quick Start
# ###########

# Most of the utility classes and functions come from lattice.py, and is 
# complemented with helper functions from delta.py
# 
# In general everything related with lattices (creation, representation
# transformation, generation, etc) can be found in lattice.py

# Basic usage 
# ===========
# 
# In this example we want to create a random lattice with 7 nodes, obtain all
# the valid space functions for the lattice, select 3 at random and calculate
# the point-wise least upper bound of these three space functions. Finally,
# since the resulting lattice may not be distributive, we use `delta_foo` to
# obtain the corresponding Group Space Function Delta_I.

covers = random_lattice(7, 0.95) # 0.95 seems to generate good random lattices.
lattice = Lattice(lattice_from_covers(covers))
sf = random.sample(lattice.space_functions, 3)
elem = [lattice.lub(f) for f in zip(*sf)]
delta = delta_foo(lattice, sf)

# Explanation:
# - We obtain the `covers` representation of a random 7-node lattice.
# - Since the `Lattice` class uses a matrix representation, and not covers, 
#   we use the `lattice_from_covers` function to transform from one
#   representation to the other.
# - Now we create our `lattice` from the `Lattice` class, which gives us access
#   to basic lattice operations, like obtaining the least upper bound of a list
#   of values, Heyting implication, atoms, etc.
# - We use the `delta_foo` function from lattice.py that receives a `Lattice`
#   instance as it first parameter and doesn't require previous set-up.

# #######
# General
# #######

# Lattices are represented by a matrix of relations (or entailment)
# `matrix[a][b] = 1` means that the element `a` entails `b` (`a` is greater
# than or equal to `b`). A value of 0 means that `a` and `b` are not comparable.
# 
# Most functions from delta.py use the matrix representation of the lattice to
# make calculations. But most of the functions from lattice.py use a special
# class named `Lattice` that uses this matrix internally.

# ##################
# Lattice Generation
# ##################

# lattice.py contains some functions that help in the generation of random
# lattices. The main ones are:
#
# - `random_lattice`: Generates a random lattice, the resulting list needs
#   to be converted with the `lattice_from_covers` function to be used with
#   the `Lattice` class.
# - `powerset_lattice`: Generates a powerset lattice that can be used directly
#   with the `Lattice` class.
#
# generation.py contains othe generation functions that may be moved back to
# lattice.py or vice-versa:
#
# - `distributive_lattice`: Generates a random distributive lattice of variable
#   size that can be used direclty with the `Lattice` class.
#
# Finally, for some basic lattices you can use the utility functions from
# delta.py, all can be used directly with the class `Lattice`:
#
# - `lattice_m`: Generates a lattice M_n, where n is the size of the lattice.
# - `lattice_n5`: Generates the non-distributive lattice N_5.
# - `lattice_kite`: Generates the kite lattice.
# - `lattice_square`: Generates a square lattice (9 elements).
# - `lattice_power_3`: Generates the powerset lattice (8 elements).

# #######################
# Lattice Representations
# #######################

# There are two ways of representing a lattice, by its relations matrix or as a
# list of covers. Usually the covers representation is more human readable. As
# such, lattice.py provides some helper functions for converting from one
# representation to another:
#
# - `lattice_from_covers`: Converts a list of lower covers to a relations matrix
#   usable, for example, by the `Lattice` class.
# - `covers_from_lattice`: Converst an implies matrix of a lattice into its equivalent list of
#   covers. Useful when printing the lattice.
# - `lattice_from_joins`: Converts a matrix of joins (least upper bounds), to a
#   relations matrix.
#
# As an example, the following to lattices are equivalent, the kite lattice:

lattice_covers = [[], [0], [0], [1,2]]
lattice_matrix = [[1, 0, 0, 0],
                  [1, 1, 0, 0],
                  [1, 0, 1, 0],
                  [1, 1, 1, 1]]

# ###############
# Space Functions
# ###############

# Generating, validating, and counting are some of the common operations one
# may need to achieve when working with space functions. 
#
# # Generation
#
# The `Lattice` class from lattice.py provides the `space_function` property,
# a lazily evaluated attribute, that provides a list of all the space functions
# valid for the lattice. This combined with the random.sample functions is the
# general way of generating and selecting space functions.
#
# When working with Powerset lattices, it is recommended to use the
# `random_space_function` function from delta.py, which receives a `Lattice`
# instance as its parameter. Furthermore, given the interesting properties of
# the powerset lattice, it is also possible to define a space function only by
# the mapping of its atoms and you can use `decode_powerset_function` to do
# just that it returns a space function based on the atoms mapping.
#
# # Validation
#
# The `Lattice` class method `is_fn_distributive` check if a "function" is
# distributive over joins.

# #################
# Computing Delta_I
# #################

# As expected most of the functions related to computing Delta_I are in the
# delta.py file, but there's a catch. Due to legacy implementations, the
# functions defined in delta.py make use of the relations matrix directly,
# and require that certain global variables are defined before use:
#
# - `LUBs`: The matrix of least upper bounds.
# - `GLBs`: The matrix of greatest lower bounds.
# - `IMPLs`: The matrix of Heyting implications.
#
# If a functions requires any of this globals to be define it will appear
# in the function documentation.
#
# delta.py contains all versions of Delta*, delta_foo, and special version of
# delta_foo `probed_delta_foo` that also returns the number of times the
# candidate function `delta` was updated. Last but not least `delta_n`
# calculates Delta_I using brute-force.
#
# One trick is to create a function that defines this globals using the
# convinient `Lattice` class, like in the following example:

def setup_enviroment(lattice):
    """ Defines LUBs, GLBs, and IMPLs global from the Lattice `lattice` """
    global LUBs, GLBs, IMPLs
    LUBs = lattice.lubs
    GLBs = lattice.glbs
    IMPLs = lattice.impls

# However, two of the "delta" functions are also defined in lattice.py and make
# use of the `Lattice` instead of the relations matrix, and need no globals:
#
# - `delta_ast_partition`: For distributive lattices, may not be the latest
#   implementation, but uses a look-up table.
# - `delta_foo`: For general lattices. The supporting structures are defined
#   in the class `FooContext` instead of being separate.

# #########
# Utilities
# #########

# From delta.py:
#
# - `relative_path`: Pass the relative path of a file to the script. Behaves
#   like `os.path.join`. The name may change.
# - `eprint`: Like `print`, but to stderr and flushing the buffer.
#
# From lattice.py:
# 
# - `get_relative_path`: Like `relative_path`. Again, may change the name.
# - `test_equality`: Pretty prints an equality check, used for basic unit tests.
# - `process_file`: Used to process lattices generated by SAGEMath, convert them to
#   their relations matrix representation, (Optionally) calculate all space
#   functions and save the results on disk.
#
# From generation.py:
#
# - `is_distributive`: Checks if a `Lattice` instance is a distributive.
# - `binomial`: Calculates the binomial coefficient using the multiplicative
#   formula.
# - `progress_bar`: For long running processes, displays a progress bar to
#   `sys.stderr`.

# ############
# Experimental
# ############

# Since computing all the space functions of a lattice is an expesive process
# experimental.py provides the `NxLattice` class, with the same inteface as the
# `Lattice` class, but allows for parallel computation. Can be used where
# `Lattice` is used.

# ###########
# Experiments
# ###########

# Almost all the experiments ran before are located in delta.py and prefixed
# with the word `run`, e.g. `run_space_functions_for_m` and can act as
# examples.

# #############
# More Examples
# #############

# Generate a Powerset lattice, and obtain a random valid space function.
lattice = Lattice(powerset_lattice(4))
sf = random_space_function(lattice) # Use only with powersets.

# Create a PDF file with the graphical representation of the lattice.
# Optionally show space functions.
#
# * Requires Graphviz, try `pip3 install graphviz`
lattice.diagram(sf).render("my_lattice.pdf")

# Generate a Random Distributive lattice
lattice = Lattice(distributive_lattice(4, 0.3))

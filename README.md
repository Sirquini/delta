# Computing Distributed Information

Provides Python 3 implementation of lattices, under `lattice.py` and proposed
"delta" functions for computing distributed information for a group of agents,
see `delta.py`. In this context *space functions* are *join-endomorphisms*.

Most of the utility classes and functions come from `lattice.py`, and is 
complemented with helper functions from `delta.py`.

In general, everything related to lattices (creation, representation
transformation, generation, etc) can be found in `lattice.py`.

> The content of this repository is based on the [theory developed](https://doi.org/10.1016/j.jlamp.2021.100674)
> in [various](https://doi.org/10.1007/978-3-030-43520-2_16)
> [papers](https://doi.org/10.1007/978-3-030-88701-8_25).


## Dependencies

This project depends on [graphviz](https://graphviz.org/) for the visualization
of lattices, and space functions.

```
python3 -m pip install graphviz
```

Optionally, we recommend installing [NumPy](https://numpy.org/index.html)
for fast matrix operations, and [Jupyter Notebooks](https://jupyter.org/index.html)
for running experiments interactively.


```
python3 -m pip install numpy notebook
```

For quick examples, please take a look at `guide.py`. At a glance.

## Basic usage 

Let us create a random lattice with 7 nodes, obtain all
the valid space functions for the lattice, select 3 at random, and calculate
the point-wise least upper bound of these three space functions. Finally,
since the resulting lattice may not be distributive, we use `delta_foo` 
(also known as DGen) to obtain the corresponding Distributed Information Space Function.

```python
import random
import lattice as lat

covers = lat.random_lattice(7, 0.95) # 0.95 seems to generate good random lattices.
lattice = lat.Lattice.from_covers(covers)
sf = random.sample(lattice.space_functions, 3)
elem = [lattice.lub(f) for f in zip(*sf)]
delta = lat.delta_foo(lattice, sf)
```

### Explanation:

- We obtain the `covers` representation of a random 7-node lattice.
- Since the `Lattice` class uses a matrix representation, and not covers,
  we use the `from_covers` class method.
- This creates our `lattice` from the `Lattice` class, which gives us access
  to basic lattice operations, like obtaining the least upper bound of a list 
  of values, Heyting implication, atoms, etc.
- We use the `delta_foo`/DGen function from `lattice.py` that receives a `Lattice`
  instance as its first parameter and doesn't require previous set-up.

## Lattices

Lattices are represented by a matrix of relations (or entailments)
`matrix[a][b] = 1` means that the element `a` entails `b` (`a` is greater
than or equal to `b`). A value of 0 means that `a` and `b` are not comparable.

Most functions from `delta.py` use the matrix representation of the lattice to
make calculations. But most of the functions from `lattice.py` use a special
class named `Lattice` that uses this matrix internally.

### Lattice Generation

`lattice.py` contains some functions that help in the generation of random
lattices. The main ones are:

- `random_lattice`: Generates a random lattice, the resulting list needs
   to be converted `from_covers` to be used with the `Lattice` class.
- `powerset_lattice`: Generates a powerset lattice that can be used directly
  with the `Lattice` class.

`generation.py` contains othe generation functions that may be moved back to
`lattice.py` or vice-versa:

- `distributive_lattice`: Generates a random distributive lattice of variable
  size that can be used direclty with the `Lattice` class.

> From `lattice.py`, `random_lattice(size, prob)` generates a random arbitrary
lattice (based on [SageMath](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/posets/poset_examples.html#sage.combinat.posets.poset_examples.Posets.RandomLattice)).

Finally, for some basic lattices you can use the utility functions from
`delta.py`, all can be used directly with the class `Lattice`:

- `lattice_m`: Generates a lattice *M_n*, where *n* is the size of the lattice.
  The traditional *M_3* would be `lattice_m(5)`.
- `lattice_n5`: Generates the non-distributive lattice *N_5*.
- `lattice_kite`: Generates the kite lattice.
- `lattice_square`: Generates a square lattice (9 elements).
- `lattice_power_3`: Generates the powerset lattice (8 elements).


## Lattice Representations

There are two ways of representing a lattice, by its relations matrix or as a
list of covers. Usually the covers representation is more human readable. As
such, `lattice.py` provides some helper functions for converting from one
representation to the other:

- `lattice_from_covers`: Converts a list of lower covers to a relations matrix.
  If you want to use the `Lattice` class, make use of the constructor `Lattice.from_covers()`.
- `covers_from_lattice`: Converts an implies matrix of a lattice into its equivalent list of
  covers. Useful when printing the lattice.
- `lattice_from_joins`: Converts a matrix of joins (least upper bounds), to a
  relations matrix.

As an example, the following two lattices are equivalent to the kite lattice:

```python
lattice_covers = [[], [0], [0], [1,2]]
lattice_matrix = [[1, 0, 0, 0],
                  [1, 1, 0, 0],
                  [1, 0, 1, 0],
                  [1, 1, 1, 1]]
```

## Space Functions

Generating, validating, and counting are some of the common operations one
may need to achieve when working with space functions. 

### Generation

The `Lattice` class from `lattice.py` provides the `space_function` property,
a lazily evaluated attribute that provides a list of all the space functions
valid for the lattice. This combined with the `random.sample` functions is the
general way of generating and selecting space functions.

When working with Powerset lattices, it is recommended to use the
`random_space_function` function from `delta.py`, which receives a `Lattice`
instance as its parameter. Furthermore, given the interesting properties of
the powerset lattice, it is also possible to define a space function only by
the mapping of its atoms. You can use `decode_powerset_function` to do
just that it returns a space function based on the atoms mapping.

### Validation

The `Lattice` class method `is_fn_distributive` checks if a "function" is
distributive over joins.

## Distributed Information

As a note, the algorithm referred to DGen+ (DGen) is implemented in `lattice.py`
as `delta_foo`, and is mean to be used with arbitrary lattices. On the other
hand, in `delta.py`, `dmeet_jies` corresponds to DMeet, and is mean to be used
with distributive lattices.

As expected most of the functions related to computing DI are in the
`delta.py` file, but there's a catch. Due to legacy implementations, the
functions defined in `delta.py` prefixed with `delta_ast` or `delta_foo` make use
of the relations matrix directly,
and require certain global variables to be present before use:

- `LUBs`: The matrix of least upper bounds.
- `GLBs`: The matrix of greatest lower bounds.
- `IMPLs`: The matrix of Heyting implications.

If a functions requires any of this globals to be define it will appear
in the function documentation.

`delta.py` contains all versions of Delta*, delta_foo, and special version of
delta_foo `probed_delta_foo` that also returns the number of times the
candidate function `delta` was updated. Last but not least `delta_n`
calculates DI using brute-force.

> This cath does not apply for the `delta_plus` functions, and `dmeet_jies`.

> A version without global state for DGen, or `delta_foo`, is available in `lattice.py`.

## Utilities

From `delta.py`:

- `relative_path`: Pass the relative path of a file to the script. Behaves
  like `os.path.join`. The name may change.
- `eprint`: Like `print`, but to stderr and flushing the buffer.

From `lattice.py`:

- `get_relative_path`: Like `relative_path`. Again, may change the name.
- `test_equality`: Pretty prints an equality check, used for basic unit tests.
- `process_file`: Used to process lattices generated by SAGEMath, convert them to
  their relations matrix representation, (Optionally) calculate all space
  functions and save the results on disk.

From `generation.py`:

- `is_distributive`: Checks if a `Lattice` instance is a distributive lattice.
- `progress_bar`: For long running processes, displays a progress bar to
  `sys.stderr`. Also see `ProgressRange()` and `ProgressRange.from_iter()`,
  usable in for-loops.

## More Examples

Generate a Powerset lattice, and obtain a random valid space function.

```python
import lattice as lat
from delta import random_space_function

lattice = lat.Lattice(lat.powerset_lattice(4))
sf = random_space_function(lattice) # Use only with powersets.
```

Create a PDF file with the graphical representation of the lattice.
Optionally show space functions.
**Requires Graphviz**

```python
lattice.diagram(sf).render("my_lattice.pdf")
```

Generate a Random Distributive lattice

```python
lattice = Lattice(distributive_lattice(4, 0.3))
```

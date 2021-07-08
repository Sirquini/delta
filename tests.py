import unittest
import lattice as lat
import delta

class TestLatticeClass(unittest.TestCase):
    def setUp(self):
        self.lattice = lat.Lattice([
            [1, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 0, 0, 0, 0],
            [1, 0, 0, 1, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 0, 0, 0],
            [1, 0, 0, 0, 0, 1, 0, 0],
            [1, 0, 1, 0, 0, 0, 1, 0],
            [1, 1, 1, 1, 1, 1, 1, 1]
        ])

    def test_len(self):
        self.assertEqual(len(self.lattice), 8)

    def test_is_lattice(self):
        self.assertTrue(self.lattice.is_lattice())

        poset = lat.Lattice([
            [1, 0, 0 ,0],
            [0, 1 ,0, 0],
            [1, 1, 1, 0],
            [1, 1, 0, 1]
        ])
        self.assertFalse(poset.is_lattice())

    def test_least_upper_bounds(self):
        actual = self.lattice.lubs
        expected = [
            [0, 1, 2, 3, 4, 5, 6, 7],
            [1, 1, 4, 4, 4, 7, 7, 7],
            [2, 4, 2, 4, 4, 7, 6, 7],
            [3, 4, 4, 3, 4, 7, 7, 7],
            [4, 4, 4, 4, 4, 7, 7, 7],
            [5, 7, 7, 7, 7, 5, 7, 7],
            [6, 7, 6, 7, 7, 7, 6, 7],
            [7, 7, 7, 7, 7, 7, 7, 7]
        ]
        self.assertEqual(actual, expected)

    def test_greatest_lower_bounds(self):
        actual = self.lattice.glbs
        expected = [
            [0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 1, 0, 0, 1],
            [0, 0, 2, 0, 2, 0, 2, 2],
            [0, 0, 0, 3, 3, 0, 0, 3],
            [0, 1, 2, 3, 4, 0, 2, 4],
            [0, 0, 0, 0, 0, 5, 0, 5],
            [0, 0, 2, 0, 2, 0, 6, 6],
            [0, 1, 2, 3, 4, 5, 6, 7]
        ]
        self.assertEqual(actual, expected)

    def test_atoms(self):
        actual = self.lattice.atoms()
        expected = [1, 2, 3, 5]
        self.assertEqual(actual, expected)

    def test_implications(self):
        self.lattice = lat.Lattice([
            [1, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 0, 0, 0, 0],
            [1, 0, 0, 1, 0, 0, 0, 0],
            [1, 1, 1, 0, 1, 0, 0, 0],
            [1, 1, 0, 1, 0, 1, 0, 0],
            [1, 0, 1, 1, 0, 0, 1, 0],
            [1, 1, 1, 1, 1, 1, 1, 1]
        ])
        actual = self.lattice.impls
        expected = [
            [0, 1, 2, 3, 4, 5, 6, 7],
            [0, 0, 2, 3, 2, 3, 6, 6],
            [0, 1, 0, 3, 1, 5, 3, 5],
            [0, 1, 2, 0, 4, 1, 2, 4],
            [0, 0, 0, 3, 0, 3, 3, 3],
            [0, 0, 2, 0, 2, 0, 2, 2],
            [0, 1, 0, 0, 1, 1, 0, 1],
            [0, 0, 0, 0, 0, 0, 0, 0]
        ]
        self.assertEqual(actual, expected)

    def test__topological_sort(self):
        actual = self.lattice._topological_sort()
        expected = (0, 1, 2, 3, 5, 6, 4, 7)

    def test_topological_order(self):
        actual = self.lattice.topological_order
        expected = (0, 1, 2, 3, 5, 6, 4, 7)

    def test_join_irreducible_elements(self):
        actual = self.lattice.join_irreducible_elements()
        expected = [1, 2, 3, 5, 6]
        self.assertEqual(actual, expected)

    def test_join_irreducibles(self):
            actual = self.lattice.join_irreducibles
            expected = [1, 2, 3, 5, 6]
            self.assertEqual(actual, expected)
    
    def test_join_irreducibles_with_covers(self):
            # Trigers faster join-irreducibles generation
            self.lattice.covers

            actual = self.lattice.join_irreducibles
            expected = [1, 2, 3, 5, 6]
            self.assertEqual(actual, expected)

    def test_covers(self):
        actual = self.lattice.covers
        expected = [[], [0], [0], [0], [1, 2, 3], [0], [2], [4, 5, 6]]
        self.assertEqual(actual, expected)

        lattice = [
            [1, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 0, 0, 0, 0],
            [1, 0, 0, 0, 1, 0, 0, 0],
            [1, 1, 0, 0, 1, 1, 0, 0],
            [1, 0, 1, 0, 1, 0, 1, 0],
            [1, 1, 1, 1, 1, 1, 1, 1]
        ]
        actual = lat.Lattice(lattice).covers
        expected = [[], [0], [0], [1, 2], [0], [1, 4], [2, 4], [3, 5, 6]]

        self.assertEqual(actual, expected)
    
    def test_from_covers(self):
        lattice = lat.Lattice.from_covers([[], [0], [0], [0], [1, 2, 3], [0], [2], [4, 5, 6]])
        actual = lattice.lattice
        expected = self.lattice.lattice
        self.assertEqual(actual, expected)

    def test_space_functions(self):
        self.lattice = lat.Lattice([
            [1, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0],
            [1, 0, 1, 0, 0, 0],
            [1, 1, 0, 1, 0, 0],
            [1, 0, 1, 0, 1, 0],
            [0, 1, 1, 1, 1, 1]
        ])
        actual = self.lattice.space_functions
        expected = [
            (0, 0, 0, 0, 0, 0), (0, 0, 1, 0, 1, 1), (0, 0, 1, 1, 1, 1),
            (0, 0, 2, 0, 2, 2), (0, 0, 2, 2, 2, 2), (0, 0, 3, 0, 3, 3),
            (0, 0, 3, 1, 3, 3), (0, 0, 3, 3, 3, 3), (0, 0, 4, 0, 4, 4),
            (0, 0, 4, 2, 4, 4), (0, 0, 4, 4, 4, 4), (0, 0, 5, 0, 5, 5),
            (0, 0, 5, 1, 5, 5), (0, 0, 5, 2, 5, 5), (0, 0, 5, 3, 5, 5),
            (0, 0, 5, 4, 5, 5), (0, 0, 5, 5, 5, 5), (0, 1, 0, 1, 0, 1),
            (0, 1, 0, 1, 1, 1), (0, 1, 1, 1, 1, 1), (0, 1, 2, 1, 2, 5),
            (0, 1, 2, 1, 4, 5), (0, 1, 2, 1, 5, 5), (0, 1, 2, 3, 2, 5),
            (0, 1, 2, 3, 4, 5), (0, 1, 2, 3, 5, 5), (0, 1, 2, 5, 2, 5),
            (0, 1, 2, 5, 4, 5), (0, 1, 2, 5, 5, 5), (0, 1, 3, 1, 3, 3),
            (0, 1, 3, 3, 3, 3), (0, 1, 4, 1, 4, 5), (0, 1, 4, 1, 5, 5),
            (0, 1, 4, 3, 4, 5), (0, 1, 4, 3, 5, 5), (0, 1, 4, 5, 4, 5),
            (0, 1, 4, 5, 5, 5), (0, 1, 5, 1, 5, 5), (0, 1, 5, 3, 5, 5),
            (0, 1, 5, 5, 5, 5), (0, 2, 0, 2, 0, 2), (0, 2, 0, 2, 2, 2),
            (0, 2, 1, 2, 1, 5), (0, 2, 1, 2, 3, 5), (0, 2, 1, 2, 5, 5),
            (0, 2, 1, 4, 1, 5), (0, 2, 1, 4, 3, 5), (0, 2, 1, 4, 5, 5),
            (0, 2, 1, 5, 1, 5), (0, 2, 1, 5, 3, 5), (0, 2, 1, 5, 5, 5),
            (0, 2, 2, 2, 2, 2), (0, 2, 3, 2, 3, 5), (0, 2, 3, 2, 5, 5),
            (0, 2, 3, 4, 3, 5), (0, 2, 3, 4, 5, 5), (0, 2, 3, 5, 3, 5),
            (0, 2, 3, 5, 5, 5), (0, 2, 4, 2, 4, 4), (0, 2, 4, 4, 4, 4),
            (0, 2, 5, 2, 5, 5), (0, 2, 5, 4, 5, 5), (0, 2, 5, 5, 5, 5),
            (0, 3, 0, 3, 0, 3), (0, 3, 0, 3, 1, 3), (0, 3, 0, 3, 3, 3),
            (0, 3, 1, 3, 1, 3), (0, 3, 1, 3, 3, 3), (0, 3, 2, 3, 2, 5),
            (0, 3, 2, 3, 4, 5), (0, 3, 2, 3, 5, 5), (0, 3, 2, 5, 2, 5),
            (0, 3, 2, 5, 4, 5), (0, 3, 2, 5, 5, 5), (0, 3, 3, 3, 3, 3),
            (0, 3, 4, 3, 4, 5), (0, 3, 4, 3, 5, 5), (0, 3, 4, 5, 4, 5),
            (0, 3, 4, 5, 5, 5), (0, 3, 5, 3, 5, 5), (0, 3, 5, 5, 5, 5),
            (0, 4, 0, 4, 0, 4), (0, 4, 0, 4, 2, 4), (0, 4, 0, 4, 4, 4),
            (0, 4, 1, 4, 1, 5), (0, 4, 1, 4, 3, 5), (0, 4, 1, 4, 5, 5),
            (0, 4, 1, 5, 1, 5), (0, 4, 1, 5, 3, 5), (0, 4, 1, 5, 5, 5),
            (0, 4, 2, 4, 2, 4), (0, 4, 2, 4, 4, 4), (0, 4, 3, 4, 3, 5),
            (0, 4, 3, 4, 5, 5), (0, 4, 3, 5, 3, 5), (0, 4, 3, 5, 5, 5),
            (0, 4, 4, 4, 4, 4), (0, 4, 5, 4, 5, 5), (0, 4, 5, 5, 5, 5),
            (0, 5, 0, 5, 0, 5), (0, 5, 0, 5, 1, 5), (0, 5, 0, 5, 2, 5),
            (0, 5, 0, 5, 3, 5), (0, 5, 0, 5, 4, 5), (0, 5, 0, 5, 5, 5),
            (0, 5, 1, 5, 1, 5), (0, 5, 1, 5, 3, 5), (0, 5, 1, 5, 5, 5),
            (0, 5, 2, 5, 2, 5), (0, 5, 2, 5, 4, 5), (0, 5, 2, 5, 5, 5),
            (0, 5, 3, 5, 3, 5), (0, 5, 3, 5, 5, 5), (0, 5, 4, 5, 4, 5),
            (0, 5, 4, 5, 5, 5), (0, 5, 5, 5, 5, 5)]
        self.assertEqual(actual, expected)

class TestHelperFunctions(unittest.TestCase):
    def test_explode(self):
        covers = [[], [0], [0], [0], [1, 2, 3], [0], [2], [4, 5, 6]]
        actual = lat.explode(covers)
        expected = [set(), {0}, {0}, {0}, {0,1,2,3}, {0}, {0, 2}, {0, 1, 2 , 3, 4, 5, 6}]
        self.assertEqual(actual, expected)

        covers = [[], [0], [0], [1, 2], [0], [1, 4],  [2, 4], [3, 5, 6]]
        actual = lat.explode(covers)
        expected = [set(), {0}, {0}, {0, 1, 2}, {0}, {0, 1, 4}, {0, 2, 4}, {0, 1, 2 , 3, 4, 5, 6}]
        self.assertEqual(actual, expected)
    
    def test_lattice_from_covers(self):
        covers = [[], [0], [0], [0], [2, 3, 1], [0], [2], [4, 5, 6]]
        actual = lat.lattice_from_covers(covers)
        expected = [
            [1, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 0, 0, 0, 0],
            [1, 0, 0, 1, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 0, 0, 0],
            [1, 0, 0, 0, 0, 1, 0, 0],
            [1, 0, 1, 0, 0, 0, 1, 0],
            [1, 1, 1, 1, 1, 1, 1, 1]
        ]
        self.assertEqual(actual, expected)

        covers = [[], [0], [0], [1, 2], [0], [1, 4], [2, 4], [3, 5, 6]]
        actual = lat.lattice_from_covers(covers)
        expected = [
            [1, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 0, 0, 0, 0],
            [1, 0, 0, 0, 1, 0, 0, 0],
            [1, 1, 0, 0, 1, 1, 0, 0],
            [1, 0, 1, 0, 1, 0, 1, 0],
            [1, 1, 1, 1, 1, 1, 1, 1]
        ]
        self.assertEqual(actual, expected)
        
    def test_covers_from_lattice(self):
        lattice = [
            [1, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 0, 0, 0, 0],
            [1, 0, 0, 1, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 0, 0, 0],
            [1, 0, 0, 0, 0, 1, 0, 0],
            [1, 0, 1, 0, 0, 0, 1, 0],
            [1, 1, 1, 1, 1, 1, 1, 1]
        ]
        actual = lat.covers_from_lattice(lattice)
        expected = [[], [0], [0], [0], [1, 2, 3], [0], [2], [4, 5, 6]]
        self.assertEqual(actual, expected)

        lattice = [
            [1, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 0, 0, 0, 0],
            [1, 0, 0, 0, 1, 0, 0, 0],
            [1, 1, 0, 0, 1, 1, 0, 0],
            [1, 0, 1, 0, 1, 0, 1, 0],
            [1, 1, 1, 1, 1, 1, 1, 1]
        ]
        actual = lat.covers_from_lattice(lattice)
        expected = [[], [0], [0], [1, 2], [0], [1, 4], [2, 4], [3, 5, 6]]

        self.assertEqual(actual, expected)
    
    def test_lattice_from_joins(self):
        joins = [
            [0, 1, 2, 3, 4, 5, 6, 7],
            [1, 1, 4, 4, 4, 7, 7, 7],
            [2, 4, 2, 4, 4, 7, 6, 7],
            [3, 4, 4, 3, 4, 7, 7, 7],
            [4, 4, 4, 4, 4, 7, 7, 7],
            [5, 7, 7, 7, 7, 5, 7, 7],
            [6, 7, 6, 7, 7, 7, 6, 7],
            [7, 7, 7, 7, 7, 7, 7, 7]
        ]
        actual = lat.lattice_from_joins(joins)
        expected = [
            [1, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 0, 0, 0, 0],
            [1, 0, 0, 1, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 0, 0, 0],
            [1, 0, 0, 0, 0, 1, 0, 0],
            [1, 0, 1, 0, 0, 0, 1, 0],
            [1, 1, 1, 1, 1, 1, 1, 1]
        ]
        self.assertEqual(actual, expected)

    def test_random_space_function(self):
        lattice = lat.Lattice.from_covers([[], [0], [0], [1, 2], [0], [1, 4], [2, 4], [3, 5, 6]])
        for _ in range(100):
            self.assertIn(delta.space_function(lattice), lattice.space_functions)

class TestDeltaFunctions(unittest.TestCase):
    def setUp(self):
        self.lattice = lat.Lattice([
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            [1, 0, 1, 0, 0, 1, 1, 0, 0, 0],
            [1, 1, 0, 0, 0, 1, 0, 1, 0, 0],
            [1, 1, 1, 1, 0, 1, 1, 1, 1, 0],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        ])

        self.scenarios = [
            {
                "fns": [
                    [0, 8, 1, 8, 9, 5, 7, 8, 8, 9],
                    [0, 4, 1, 4, 4, 9, 9, 9, 9, 9],
                    [0, 6, 9, 9, 9, 3, 9, 8, 9, 9],
                    [0, 7, 8, 8, 9, 1, 8, 7, 8, 9],
                    [0, 1, 0, 1, 1, 9, 9, 9, 9, 9]
                ],
                "exp": [0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
            },
            {
                "fns": [
                    [0, 5, 6, 6, 8, 3, 8, 8, 8, 8],
                    [0, 5, 4, 9, 9, 9, 9, 9, 9, 9],
                    [0, 8, 5, 8, 8, 0, 5, 8, 8, 8],
                    [0, 2, 8, 8, 8, 7, 8, 8, 8, 8],
                    [0, 5, 4, 9, 9, 5, 9, 5, 9, 9]
                ],
                "exp": [0, 0, 0, 0, 8, 0, 0, 0, 0, 8]
            },
            {
                "fns": [
                    [0, 4, 8, 9, 9, 9, 9, 9, 9, 9],
                    [0, 4, 7, 9, 9, 8, 8, 9, 9, 9],
                    [0, 2, 4, 4, 9, 1, 4, 3, 4, 9],
                    [0, 7, 4, 9, 9, 8, 9, 8, 9, 9],
                    [0, 5, 2, 6, 8, 4, 4, 9, 9, 9]
                ],
                "exp": [0, 0, 0, 0, 8, 1, 1, 1, 1, 8]
            },
            {
                "fns": [
                    [0, 6, 4, 9, 9, 4, 4, 9, 9, 9],
                    [0, 2, 4, 4, 4, 6, 9, 6, 9, 9]
                ],
                "exp": [0, 2, 4, 4, 4, 2, 4, 2, 4, 4]},
            {
                "fns": [
                    [0, 9, 6, 9, 9, 5, 6, 9, 9, 9],
                    [0, 4, 2, 4, 4, 8, 8, 9, 9, 9]
                ],
                "exp": [0, 4, 2, 4, 4, 5, 6, 9, 9, 9]},
            {
                "fns": [
                    [0, 7, 7, 7, 7, 3, 8, 8, 8, 8],
                    [0, 5, 8, 8, 8, 9, 9, 9, 9, 9]
                ],
                "exp": [0, 5, 7, 7, 7, 3, 8, 8, 8, 8]},
            {
                "fns": [
                    [0, 2, 7, 8, 9, 7, 7, 8, 8, 9],
                    [0, 4, 9, 9, 9, 2, 9, 4, 9, 9]
                ],
                "exp": [0, 2, 7, 8, 9, 0, 7, 2, 8, 9]},
            {
                "fns": [
                    [0, 5, 9, 9, 9, 0, 9, 5, 9, 9],
                    [0, 1, 9, 9, 9, 0, 9, 1, 9, 9]
                ],
                "exp": [0, 0, 9, 9, 9, 0, 9, 0, 9, 9]
            }
        ]

    def test_delta_foo(self):
        for scenario in self.scenarios:
            actual = lat.delta_foo(self.lattice, scenario["fns"])
            self.assertEqual(actual, scenario["exp"])

    def test_delta_partition(self):
        for scenario in self.scenarios:
            actual = lat.delta_partition(self.lattice, scenario["fns"])
            self.assertEqual(actual, scenario["exp"])

        jies = self.lattice.join_irreducibles
        actual = lat.delta_partition(self.lattice, scenario["fns"], jies)
        self.assertEqual(actual, scenario["exp"])

    def test_delta_plus_jies(self):
        for scenario in self.scenarios:
            actual = delta.delta_plus_jies(self.lattice, scenario["fns"])
            self.assertEqual(actual, scenario["exp"])

        jies = self.lattice.join_irreducibles
        actual = delta.delta_plus_jies(self.lattice, scenario["fns"], jies)
        self.assertEqual(actual, scenario["exp"])
    
    def test_dmeet_jies(self):
        for scenario in self.scenarios:
            actual = delta.dmeet_jies(self.lattice, scenario["fns"])
            self.assertEqual(actual, scenario["exp"])

        jies = self.lattice.join_irreducibles
        actual = delta.dmeet_jies(self.lattice, scenario["fns"], jies)
        self.assertEqual(actual, scenario["exp"])

        covers = self.lattice.covers
        actual = delta.dmeet_jies(self.lattice, scenario["fns"], covers=covers)
        self.assertEqual(actual, scenario["exp"])

        actual = delta.dmeet_jies(self.lattice, scenario["fns"], jies, covers)
        self.assertEqual(actual, scenario["exp"])

if __name__ == "__main__":
    unittest.main()

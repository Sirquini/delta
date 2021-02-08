from poset import Poset
from generation import distributive_lattice as _distributive_lattice

import numpy as np
import time
import random, math, time
from itertools import combinations
from collections import deque, Counter


class LatticeCollection(list):
    
    def __getitem__(self, n):
        'List of all lattices of size n'
        if isinstance(n, int) and n>=len(self):
            while len(self) <= n:
                self.grow()
        return super().__getitem__(n)
    
    def grow(self):
        'Find and store lattices of size len(self)'
        if len(self) > 2:
            layer = list(self.grow_iter())
        else:
            layer = {
                0: [],
                1: [Poset.from_covers([[]])],
                2: [Poset.from_covers([[],[0]])],
            }[len(self)]
        self.append(layer)
    
    def grow_iter(self):
        'Generate lattices by adding one node and zero or more edges to those of self[-1]'
        initial = self[-1]
        found = []
        vis = set(initial)
        for u in initial:
            for v in u.iter_add_node():
                if v not in vis:
                    vis.add(v)
                    found.append(v)
                    yield v
        q = deque(list(found))
        while q:
            u = q.popleft()
            for v in u.iter_add_edge():
                if v not in vis:
                    vis.add(v)
                    found.append(v)
                    q.append(v)
                    yield v
    
    def flatten(self):
        return [V for l in self for V in l]

        
    @classmethod
    def up_to_n(cls, n, verbose=False):
        C = cls()
        C.grow()
        N = n
        for n in range(N):
            start = time.time()
            C.grow()
            if verbose:
                elapsed = time.time()-start
                print(f'[{elapsed:7.2f}s] n={n}, n_lattices={len(C[n])}', flush=True)
        return C
    
    @classmethod
    def flatten_up_to_n(cls, n):
        
        return


def test_count_f_lub_distributive(N=8):
    times = [[] for i in range(N+1)]
    try:
        for V in LatticeCollection.up_to_n(N).flatten():
            if V.is_distributive:
                t0 = time.time()
                a =  V.count_f_lub_distributive_v2()
                t1 = time.time()
                b = V.count_f_lub_bruteforce()
                t2 = time.time()
                times[V.n].append((t1-t0,t2-t1))
                assert a==b, (a,b,V.show())
    finally:
        print('Timings:')
        for n in range(N+1):
            if times[n]:
                d0,d1 = (np.array(times[n]).mean(axis=0)*1000).round(0)
                print(f'n={n:2d} | avg_v2={d0:5.0f}ms | avg_bf={d1:5.0f}ms')
    print(f'-----\nThe v2 algorithm is correct up to n={N} inclusive')
    return

def distributive_lattice(n, p):
    leq = np.array(_distributive_lattice(n, p), dtype=bool)
    leq.flags.writeable = False
    return Poset(leq)


def large_test():
    #https://oeis.org/A006982/list
    random.seed(0)
    found = {}
    for n in [2,3,4,5,6,7]:
        for _ in range(10 * n**2):
            p = random.random()
            V = distributive_lattice(n, p)
            found[V.hash] = V
    found = set(found.values())
    len(found)
    
    print('Conjecture to be tested:\n\n    n^log2(n) <= |F_lub(L)| <= comb(2n-2, n-1)')
    print('\n    for all lattices L, where n is its size (including bottom)')
    print('    and F_lub(L) is the set of endomorphisms that preserve least upper bounds')


    sizes = []
    for size, examples in sorted(Counter(V.n for V in found).items()):
        sizes.append(f'{size}(x{examples})' if examples>1 else f'{size}')
    print(f'\nSize of the examples to be tested ({len(found)} in total):')
    print('\n   n =', ', '.join(sizes))
    print('\n-----\n')

    start_main = time.time()
    try:
        for V in sorted(found, key=lambda V:V.n):
            lower = V.n**math.log2(V.n)
            upper = math.comb(2*V.n-2, V.n-1)
            print(f'n={V.n}')
            print(f'    Hash:       {hash(V)}')
            print(f'    Lattice:    {V}')
            start = time.time()
            middle = V.num_f_lub
            elapsed = round(time.time()-start, 2)
            print(f'    |F_lub(L)|: {middle}')
            print(f'    Elapsed:    {elapsed}s')
            print(f'    Conjecture: {round(lower,1)} <= {middle} <= {upper}')
            print(f'    Conjecture: {lower <= middle <= upper}')
            print(end='', flush=True)
            assert lower <= middle <= upper, (V.n,lower,middle,upper,V,V.show())
            print()
        print(f'No counter example was found.')
    finally:
        print(f'Total time elapsed: {round(time.time()-start_main, 2)}s')
    return found


LatticeCollection.up_to_n(8, verbose=True)
test_count_f_lub_distributive(7)
large_test()
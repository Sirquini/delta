from cached_property import cached_property
import pyhash
import numpy as np
from collections import deque
from itertools import product, groupby, chain
from functools import reduce
import time, sys



def product_list(*iterables, repeat=1, out=None):
    'same as itertools.product, but mutates the output instead of making tuples'
    dims = [list(it) for it in iterables] * repeat
    n = len(dims)
    if out is not None:
        assert len(out)==n, f'Incompatible output shape'
    out = [None]*n if out is None else out
    def backtrack(i):
        if i==n:
            yield out
        else:
            for x in dims[i]:
                out[i] = x
                yield from backtrack(i+1)
    yield from backtrack(0)


class Outfile:
    def __init__(self, outfile=None):
        self.outfile = outfile
        
    def __enter__(self):
        if self.outfile is not None:
            self.initial_stdout = sys.stdout
            sys.stdout = open(self.outfile, 'a')
    
    def __exit__(self, *args):
        if self.outfile is not None:
            sys.stdout.close()
            sys.stdout = self.initial_stdout


class PosetException(Exception):
    pass


class Poset:
    """
    Hashable object that represents an inmutable finite partial order.
    Hash is invariant under permutations.
    
    Requires external packages:
        - cached_property
        - pyhash
        - pydotplus (which needs graphviz 'dot' program)
    
    The main attributes (always present) are:
        - n: size of the poset. The elements of the poset are range(n)
        - leq: read only (inmutable) boolean nxn matrix. leq[i,j]==True iff i <= j
    
    All other attributes are lazy loaded and usually cached.
    
    Conventions:
        - child[i,j]==True iff j covers i (with no elements inbetween)
        - children[j] = [i for i in range(n) if leq[i,j]]
        - parents[i] = [j for j in range(n) if leq[i,j]]

        For lattices:
            - lub[i,j] is the least upper bound for i and j.
            - glb[i,j] is the greatest lower bound for i and j
    
    Why pyhash?
        Because it is stable (like hashlib) and fast (like hash).
        hashlib is not adequate because it adds an unnecessary computation footrint.
        hash(tuple(...)) is not adequate because it yields different results across
        several runs unless PYTHONHASHSEED is set prior to execution.
    
    Example:
        V = Poset.from_children([[],[0],[0],[1,2]])
        V.show()
        print(V.num_f_lub)
            for f in V.iter_f_lub_bruteforce():
            V.show(f)
    """
    
    def __init__(self, leq):
        'assumes that leq is indeeed a partial order'
        assert leq.dtype==bool, 'leq must be a boolean numpy array'
        assert leq.flags.writeable == False, 'leq must be read-only'
        n = leq.shape[0]
        assert tuple(leq.shape)==(n,n), f'leq must be squared {leq.shape}'
        self.n = n
        self.leq = leq
    
    # Representation methods
    
    @cached_property
    def child(self):
        'nxn boolean matrix. out[i,j] iff j covers i (with no elements inbetween)'
        return self.__class__.leq_to_child(self.leq)
    
    @classmethod
    def leq_to_child(cls, leq):
        'Compute child (aka cover) relation from the poset relation'
        n = len(leq)
        lt = leq.copy()
        lt[np.diag_indices_from(lt)] = False
        any_inbetween = np.matmul(lt, lt)
        return lt & ~any_inbetween
    
    def __repr__(self):
        s = repr(self.children)
        return f'P({s[1:-1]})'
    
    def show(self, f=None, as_edges=False, save=None, labels=None):
        'Use graphviz to display or save self (or the endomorphism f if given)'
        g = self.graphviz(f, as_edges, labels)
        png = g.create_png()
        
        if save is None:
            from IPython.display import display
            from IPython.display import Image
            img = Image(png)
            display(img)
        else:
            with open(save, 'wb') as f:
                f.write(png)
        return
    
    def graphviz(self, f=None, as_edges=False, labels=None):
        'Graphviz representation of self (or f if given)'
        n = self.n
        child = self.child
        extra_edges = None
        if labels is None:
            labels = range(n)
        if f is not None:
            n = self.n
            if as_edges:
                extra_edges = [(i,f[i]) for i in range(n)]
            else:
                labels = ['' for i in range(n)]
                for i, l in groupby(range(n), f.__getitem__):
                    labels[i] = ','.join(map(str,l))
        return self._graphviz(labels, extra_edges)
    
    def _graphviz(self, labels, extra_edges):
        n = self.n
        child = self.child
        
        from pydotplus import graph_from_edges
        from pydotplus.graphviz import Node, Edge

        g = graph_from_edges([], directed=True)
        g.set_rankdir('BT')
        for i in range(n):
            style = {}
            g.add_node(Node(i, label=f'"{labels[i]}"', **style))
        for i in range(n):
            for j in range(n):
                if child[i,j]:
                    g.add_edge(Edge(i,j))
        if extra_edges is not None:
            for i,j in extra_edges:
                style = {'color':'blue'}
                g.add_edge(Edge(i, j, **style))
        return g
    
    def throw(self, message):
        print(message)
        self.show()
        print('Covers:', self)
        print('Relation matrix:')
        print(self.leq.astype(int))
        raise PosetException(message)
     
    @cached_property
    def children(self):
        '''(aka. covers): top-down adjoint list (j in G[i] iff i covers j)'''
        n = self.n
        child = self.child
        return [[j for j in range(n) if child[j,i]] for i in range(n)]
     
    @cached_property
    def parents(self):
        '''bottom-up adjoint list (j in G[i] iff j covers i)'''
        n = self.n
        child = self.child
        return [[j for j in range(n) if child[i,j]] for i in range(n)]
    
    
    # Interface methods
    
    @classmethod
    def from_parents(cls, parents):
        n = len(parents)
        children = [[] for i in range(n)]
        for ch in range(n):
            for pa in parents[ch]:
                children[pa].append(ch)
        return cls.from_children(children)
    
    @classmethod
    def from_children(cls, children):
        n = len(children)
        child = np.zeros((n,n), dtype=bool)
        for pa in range(n):
            for ch in children[pa]:
                child[ch,pa] = True
        child.flags.writeable = False
        dist = cls.child_to_dist(child)
        dist.flags.writeable = False
        leq = dist < n
        leq.flags.writeable = False
        poset = cls(leq)
        poset.is_partial_order(leq) or poset.throw('Not a partial order')
        poset.__dict__['child'] = child
        poset.__dict__['dist'] = dist
        return poset
    
    @classmethod
    def is_partial_order(cls, rel):
        "Check if the given relation is transitive, reflexive and antysimetric"
        if not rel[np.diag_indices_from(rel)].all():
            return False # reflexivity
        if (rel&rel.T).sum() > len(rel):
            return False # antysimmetry
        rel2 = np.matmul(rel, rel)
        if ((~rel)&rel2).any():
            return False # transitivity
        return True
    
    @cached_property
    def dist(self):
        return self.__class__.child_to_dist(self.child)
    
    @classmethod
    def child_to_dist(cls, child):
        'Compute all pairs shortest distances using Floyd-Warshall algorithm'
        dist = child.astype(np.uint64)
        n = len(dist)
        dist[dist==0] = n
        dist[np.diag_indices_from(dist)] = 0
        for k in range(n):
            np.minimum(dist, dist[:,k,None] + dist[None,k,:], out=dist)
        dist.flags.writeable = False
        return dist
    
    # Graph structure methods
    
    def subgraph(self, domain):
        n = self.n
        m = len(domain)
        assert len(set(domain))==m<=n, f'Invalid domain: {domain}'
        leq = self.leq
        sub = np.zeros((m,m), dtype=bool)
        for i in range(m):
            for j in range(m):
                sub[i,j] = leq[domain[i], domain[j]]
        sub.flags.writeable = False
        return self.__class__(sub)
    
    @cached_property
    def toposort(self):
        n = self.n
        G = self.parents
        child = self.child
        indeg = [child[:,i].sum() for i in range(n)]
        topo = []
        q = deque([i for i in range(n) if indeg[i]==0])
        while q:
            u = q.popleft()
            topo.append(u)
            for v in G[u]:
                indeg[v] -= 1
                if indeg[v]==0:
                    q.append(v)
        len(topo)==n or self.throw('There is a cycle')
        return topo
    

    @cached_property
    def independent_components(self):
        'Graph components if all edges were bidirectional'
        n = self.n
        cmp = self.leq|self.leq.T
        G = [[j for j in range(n) if cmp[i,j]] for i in range(n)]
        color = np.ones(n, dtype=int)*-1
        def component(i):
            q = deque([i])
            found = []
            while q:
                u = q.popleft()
                for v in G[u]:
                    if color[v]!=color[u]:
                        color[v] = color[u]
                        q.append(v)
                found.append(u)
            return found
        comps = []
        for i in range(n):
            if color[i]==-1:
                color[i] = len(comps)
                comps.append(component(i))
        return comps
    
    # Lattice methods
    
    def assert_lattice(self):
        if self.n > 0:
            self.lub
            self.bottom
    
    @cached_property
    def lub(self):
        n = self.n
        leq = self.leq
        lub_id = {tuple(leq[i,:]):i for i in range(n)}
        lub = np.zeros((n,n), int)
        for i in range(n):
            for j in range(n):
                above = tuple(leq[i,:] & leq[j,:])
                above in lub_id or self._throw_lattice(i,j)
                lub[i,j] = lub_id[above]
        lub.flags.writeable = False
        return lub
    
    def _throw_lattice(self, i, j):
        'Throw explaining why self is not a lattice by looking at i and j'
        n = self.n
        leq = self.leq
        above = [k for k in range(n) if leq[i,k] and leq[j,k]]
        below = [k for k in range(n) if leq[k,i] and leq[k,j]]
        above or self.throw(f'Not a lattice: {i} lub {j} => (no common ancestor)')
        below or self.throw(f'Not a lattice: {i} glb {j} => (no common descendant)')
        lub = min(above, key=lambda k: sum(leq[:,k]))
        glb = max(below, key=lambda k: sum(leq[:,k]))
        for x in above:
            leq[lub,x] or self.throw(f'Not a lattice: {i} lub {j} => {lub} or {x}')
        for x in below:
            leq[x,glb] or self.throw(f'Not a lattice: {i} glb {j} => {glb} or {x}')
    
    @cached_property
    def bottom(self):
        '''bottom element of the Poset. Throws if not present'''
        n = self.n
        nleq = self.leq.sum(axis=0)
        zeros = [i for i in range(n) if nleq[i]==1]
        zeros or self.throw(f'No bottom found')
        len(zeros)==1 or self.throw(f'Multiple bottoms found: {zeros}')
        return zeros[0]
    
    @cached_property
    def top(self):
        '''top element of the Poset. Throws if not present'''
        n = self.n
        nleq = self.leq.sum(axis=0)
        zeros = [i for i in range(n) if nleq[i]==n]
        zeros or self.throw(f'No top found')
        len(zeros)==1 or self.throw(f'Multiple bottoms found: {zeros}')
        return zeros[0]
    
    @cached_property
    def irreducibles(self):
        n = self.n
        children = self.children
        return [i for i in range(n) if len(children[i])==1]
    
    @cached_property
    def glb(self):
        n = self.n
        geq = self.leq.T
        glb_id = {tuple(geq[i,:]):i for i in range(n)}
        glb = np.zeros((n,n), int)
        for i in range(n):
            for j in range(n):
                below = tuple(geq[i,:] & geq[j,:])
                below in glb_id or self._throw_lattice(i,j)
                glb[i,j] = glb_id[below]
        glb.flags.writeable = False
        return glb
    
    
    # Hashing and isomorphisms
    
    _hasher = pyhash.xx_64(seed=0)
    @classmethod
    def hasher(cls, ints):
        'Fast hash that is consistent across runs independently of PYTHONHASHSEED'
        return cls._hasher(str(ints)[1:-1])>>1 # Prevent uint64->int64 overflow
    
    def hash_perm_invariant(self, mat):
        HASH = self.__class__.hasher
        h = lambda l: HASH(sorted(l))
        a = [HASH((h(mat[:,i]), h(mat[i,:]))) for i in range(self.n)]
        return np.array(a, dtype=int)
    
    @cached_property
    def hash_elems(self):
        mat = self.leq.astype(np.int64)
        with np.errstate(over='ignore'):
            H = self.hash_perm_invariant(mat)
            for repeat in range(2):
                mat += np.matmul(H[:,None], H[None,:])
                H = self.hash_perm_invariant(mat)
        return H
    
    @cached_property
    def hash(self):
        return self.__class__.hasher(sorted(self.hash_elems))
    
    def __hash__(self):
        return self.hash
    
    def __eq__(self, other):
        'Equality up to isomorphism, i.e. up to relabeling'
        N_NO_HASH_COLLISIONS_TESTED = 10
        if self.n == other.n <= N_NO_HASH_COLLISIONS_TESTED:
            eq = hash(self) == hash(other)
        else:
            eq = self.find_isomorphism(other) is not None
        return eq
    
    def find_isomorphism(self, other):
        # Quick check:
        if self.n!=other.n or hash(self)!=hash(other):
            return None
        
        # Filter out some functions:
        n = self.n
        Ah = self.hash_elems
        Bh = other.hash_elems
        
        matches = [[j for j in range(n) if Ah[i]==Bh[j]] for i in range(n)]
        remaining = product_list(*matches)
        
        # Find isomorphism among remaining functions
        A = self.leq
        B = other.leq
        
        def is_isomorphism(f):
            return all(A[i,j]==B[f[i],f[j]] for i in range(n) for j in range(n))
        
        return next((f for f in remaining if is_isomorphism(f)), None)
    
    
    def relabel(self, f, inverse=False):
        'Relabelled copy of self that i is to self as f[i] to out'
        'If inverse==True, then f[i] is to self as i to out'
        n = self.n
        assert len(f)==n and sorted(set(f))==list(range(n)), f'Invalid permutation {f}'
        if inverse:
            inv = [0]*n
            for i in range(n):
                inv[f[i]] = i
            f = inv
        leq = self.leq
        out = np.zeros_like(leq)
        for i in range(n):
            for j in range(n):
                out[f[i],f[j]] = leq[i,j]
        out.flags.writeable = False
        return self.__class__(out)
    
    @cached_property
    def canonical_mapping(self):
        n = self.n
        leq = self.leq
        h = self.hash_elems
        def key(i):
            return (dist[0,i], leq[:,i].sum())
        g = sorted(range(n), key=lambda i: (leq[:,i].sum(), h[i], ))
        f = list(range(n))
        for i in range(n):
            f[g[i]] = i
        return f
    
    @cached_property
    def is_canonical(self):
        is_sorted = lambda l: all(a<=b for a,b in zip(l,l[1:]))
        return is_sorted(self.canonical_mapping)
    
    @cached_property
    def canonical(self):
        'Relabelled copy'
        return self.relabel(self.canonical_mapping)

    
    # Methods for atomic changes (grow-by-one)
    
    @cached_property
    def forbidden_pairs(self):
        "Pairs (i,j) that break lub uniqueness or partial order structure if i<=j is assumed"
        n = self.n
        leq = self.leq
        joi = self.lub
        nocmp = ~(leq + leq.T)
        def f(a,b):
            if leq[b,a]: return True
            if leq[a,b]: return False
            X = [x for x in range(n) if leq[x, a]]
            Y = [y for y in range(n) if ~leq[b,y] and nocmp[y,a]]
            return any(nocmp[joi[x,y],joi[b,y]] for y in Y for x in X)
        fb = np.array([[f(i,j) for j in range(n)] for i in range(n)], dtype=bool)
        return fb
    
    def iter_add_edge(self):
        "Grow self by adding one edge"
        n = self.n
        leq = self.leq
        fb = self.forbidden_pairs
        vis = set()
        h = self.hash_elems
        for i,j in product_list(range(n), repeat=2):
            if not fb[i,j] and not leq[i,j] and not (h[i],h[j]) in vis:
                vis.add((h[i],h[j]))
                new_leq = leq + np.matmul(leq[:, i:i+1], leq[j:j+1, :])
                new_leq.flags.writeable = False
                yield self.__class__(new_leq)
        return
    
    def iter_add_node(self):
        "Grow self by adding one node"
        n = self.n
        leq = self.leq
        new_leq = np.zeros((n+1, n+1), bool)
        new_leq[:-1,:-1] = leq
        new_leq[n,n] = True
        fb = self.forbidden_pairs
        vis = set() # Don't repeat isomorphical connections
        h = self.hash_elems
        for i,j in product_list(range(n), repeat=2):
            if not fb[i,j] and not (h[i],h[j]) in vis:
                vis.add((h[i],h[j]))
                out = new_leq.copy()
                out[:-1, :-1] += np.matmul(leq[:, i:i+1], leq[j:j+1, :])
                out[n, :-1] = leq[j,:]
                out[:-1, n] = leq[:,i]
                out.flags.writeable = False
                yield self.__class__(out)
        return
    
    @classmethod
    def iter_all_latices(cls, max_size):
        q = deque([cls.from_children(x) for x in [[],[[]],[[],[0]]]])
        vis = set()
        while q:
            U = q.popleft()
            yield U.canonical
            it = U.iter_add_node() if U.n < max_size else iter([])
            for V in chain(U.iter_add_edge(), it):
                if V not in vis:
                    vis.add(V)
                    q.append(V)
    
    @classmethod
    def all_latices(cls, max_size):
        return list(cls.iter_all_latices(max_size))
    
    
    # Methods for all endomorphisms
    
    def iter_f_all(self):
        'all endomorphisms'
        return product_list(range(self.n), repeat=self.n)
    
    @cached_property
    def num_f_all(self):
        return self.n**self.n
        
    def iter_f_all_bottom(self):
        'all endomorphisms f with f[bottom]=bottom'
        n = self.n
        if n>0:
            options = [range(n) if i!=self.bottom else [i] for i in range(n)]
            for f in product_list(*options):
                yield f
        return
    
    @cached_property
    def num_f_all_bottom(self):
        return self.n**(self.n-1)
    
    
    # Methods for all monotonic endomorphisms
    
    def f_is_monotone(self, f, domain=None):
        'check if f is monotone over domain'
        n = self.n
        domain = range(n) if domain is None else domain
        leq = self.leq
        for i in domain:
            for j in domain:
                if leq[i,j] and not leq[f[i],f[j]]:
                    return False
        return True

    def iter_f_monotone_bruteforce(self):
        'all monotone functions'
        for f in self.iter_f_all():
            if self.f_is_monotone(f):
                yield f
        return

    def iter_f_monotone_bottom_bruteforce(self):
        'all monotone functions with f[bottom]=bottom'
        for f in self.iter_f_all_bottom():
            if self.f_is_monotone(f):
                yield f
        return

    def iter_f_monotone(self):
        'all monotone functions'
        f = [None]*self.n
        yield from self.iter_f_monotone_restricted(_f=f)
    
    def iter_f_lub_bottom_bruteforce(self):
        'all space functions. Throws if no bottom'
        for f in self.iter_f_monotone_bottom():
            if self.f_is_lub(f):
                yield f
        return
    
        
    def iter_f_monotone_restricted(self, domain=None, _f=None):
        'generate all monotone functions f : domain -> self, padding non-domain with None'
        n = self.n
        leq = self.leq
        geq_list = [[j for j in range(n) if leq[i,j]] for i in range(n)]
        f = [None for i in range(n)] if _f is None else _f
        topo, children = self._toposort_children(domain)
        yield from self._iter_f_monotone_restricted(f, topo, children, geq_list)
    
    def _iter_f_monotone_restricted(self, f, topo, children, geq_list):
        n = self.n
        m = len(topo)
        lub = self.lub
        _lub_f = (lambda acum,b: lub[acum,f[b]])
        lub_f = lambda elems: reduce(_lub_f, elems, self.bottom)
        def backtrack(i):
            'f[topo[j]] is fixed for all j<i. Backtrack f[topo[k]] for all k>=i, k<m'
            if i==m:
                yield f
            else:
                for k in geq_list[lub_f(children[i])]:
                    f[topo[i]] = k
                    yield from backtrack(i+1)
        yield from backtrack(0)
    
    def _toposort_children(self, domain):
        'Compute a toposort for domain and the children lists filtered for domain'
        'j in out.children[i] iff j in out.topo and j is children of out.topo[i]'
        n = self.n
        D = range(n) if domain is None else domain
        topo = [i for i in self.toposort if i in D]
        sub = self.subgraph(topo)
        children = [[topo[j] for j in l] for l in sub.children]
        return topo, children
    

    def iter_f_monotone_bottom(self):
        'all monotone functions with f[bottom]=bottom'
        if not self.n:
            return
        f = [None]*self.n
        f[self.bottom] = self.bottom
        domain = [i for i in range(self.n) if i!=self.bottom]
        yield from self.iter_f_monotone_restricted(domain=domain, _f=f)
    
    
    # Methods for monotonic endomorphisms over irreducibles
        
    @cached_property
    def irreducible_components(self):
        'components of join irreducibles in toposort order and children lists for each component'
        n = self.n
        if n<=1: # no join irreducibles at all
            return (0, [], [])
        irr = self.irreducibles
        sub = self.subgraph(irr)
        subcomps = sub.independent_components
        m = len(subcomps)
        irrcomps = [[irr[j] for j in subcomps[i]] for i in range(m)]
        m_topo, m_children = zip(*(self._toposort_children(irrcomps[i]) for i in range(m)))
        return m, m_topo, m_children
        
    def _interpolate_funcs(self, funcs, domain, iter_bot=False):
        'extend each f in funcs outside domain using f[j]=lub(f[i] if i<=j and i in domain)'
        n = self.n
        lub = self.lub
        leq = self.leq
        bot = self.bottom
        no_domain = [i for i in range(n) if i not in domain]
        dom_leq = [[i for i in domain if leq[i,j]] for j in range(n)]
        lub_f = (lambda a,b: lub[a,b])
        for f in funcs:
            for j in no_domain:
                f[j] = reduce(lub_f, (f[x] for x in dom_leq[j]), bot)
            yield f
    
    def iter_f_irreducibles_monotone_bottom(self):
        'all functions given by f[non_irr]=lub(f[irreducibles] below non_irr)'
        if self.n == 0:
            return
        n = self.n
        leq = self.leq
        geq_list = [[j for j in range(n) if leq[i,j]] for i in range(n)]
        m, m_topo, m_children = self.irreducible_components
        f = [None for i in range(n)]
        def backtrack(i):
            if i==m:
                yield f
            else:
                for _ in self._iter_f_monotone_restricted(f, m_topo[i], m_children[i], geq_list):
                    yield from backtrack(i+1)
        funcs = backtrack(0)
        yield from self._interpolate_funcs(funcs, self.irreducibles)
        
    
    def iter_f_irreducibles_monotone(self):
        'all functions given by f[non_irr]=lub(f[irreducibles] below non_irr) and'
        'f[bottom] = any below or equal to glb(f[irreducibles])'
        n = self.n
        if n==0:
            return
        glb = self.glb
        _glb_f = (lambda acum,b: glb[acum,f[b]])
        glb_f = lambda elems: reduce(_glb_f, elems, self.top)
        leq = self.leq
        below = [[i for i in range(n) if leq[i,j]] for j in range(n)]
        bottom = self.bottom
        irreducibles = self.irreducibles
        for f in self.iter_f_irreducibles_monotone_bottom():
            for i in below[glb_f(irreducibles)]:
                f[bottom] = i
                yield f
    
    
    # Methods for endomorphisms that preserve lub
    
    def f_is_lub_bottom(self, f, domain=None):
        'check if f preserves lub and f[bottom] is bottom. Throws if no bottom'
        n = self.n
        if n==0 or (domain is not None and len(domain)<=1):
            return True
        bot = self.bottom
        if f[bot] != bot or (domain is not None and bot not in domain):
            return False
        return self.f_is_lub(f, domain)
    
    def f_is_lub(self, f, domain=None):
        'check if f preserves lub'
        n = self.n
        domain = range(n) if domain is None else domain
        lub = self.lub
        for i in domain:
            for j in domain:
                if f[lub[i,j]]!=lub[f[i],f[j]]:
                    return False
        return True

    def iter_f_lub_bruteforce(self):
        'all join endomorphisms'
        for f in self.iter_f_monotone():
            if self.f_is_lub(f):
                yield f
        return
                
    def iter_f_lub(self):
        'all lub preserving functions'
        it = self.iter_f_irreducibles_monotone()
        if self.is_distributive:
            yield from it
        else:
            for f in it:
                if self.f_is_lub(f):
                    yield f
    
    def iter_f_lub_bottom(self):
        'all lub preserving functions with f[bottom]=bottom'
        it = self.iter_f_irreducibles_monotone_bottom()
        if self.is_distributive:
            yield from it
        else:
            for f in it:
                if self.f_is_lub(f):
                    yield f
    
    
    @cached_property
    def num_f_lub(self):
        return self.count_f_lub_bruteforce()
    
    def count_f_lub_bruteforce(self):
        return sum(1 for f in self.iter_f_lub())
    
    @cached_property
    def num_f_lub_bottom(self):
        return self.count_f_lub_bottom()
    
    def count_f_lub_bottom(self):
        if self.is_distributive:
            num = self.count_f_lub_bottom_distributive()
        else:
            num = self.count_f_lub_bottom_bruteforce()
        return num
    
    def count_f_lub_bottom_bruteforce(self):
        return sum(1 for f in self.iter_f_lub_bottom())
    
    
    # Optimizations for distributive lattices
    
    @cached_property
    def is_distributive(self):
        return self._distributive_error is None
    
    @cached_property
    def _distributive_error(self):
        'Find i, j, k that violate distributivity. None otherwise'
        n = self.n
        lub = self.lub
        glb = self.glb
        for i in range(n):
            diff = glb[i,lub] != lub[np.ix_(glb[i,:], glb[i,:])]
            if diff.any():
                j, k = next(zip(*np.where(diff)))
                return (
                    f'Non distributive lattice:\n'
                    f'{i} glb ({j} lub {k}) = {i} glb {lub[j,k]} = '
                    f'{glb[i,lub[j,k]]} != {lub[glb[i,j],glb[i,k]]} = '
                    f'{glb[i,j]} lub {glb[i,k]} = ({i} glb {j}) lub ({i} glb {k})'
                )
        return None
    
    def assert_distributive(self):
        self.is_distributive or self.throw(self._distributive_error)
    
    def iter_f_lub_bottom_distributive(self):
        'generate and interpolate all monotone functions over irreducibles'
        self.assert_distributive()
        yield from self.iter_f_irreducibles_monotone_bottom()

    def count_f_lub_bottom_distributive(self):
        self.assert_distributive()
        if self.n == 0:
            return 0
        n = self.n
        leq = self.leq
        geq_list = [[j for j in range(n) if leq[i,j]] for i in range(n)]
        m, m_topo, m_children = self.irreducible_components
        f = [None for i in range(n)]
        def num(i):
            'num of monotone functions restricted to domain k_topo[i]'
            it = self._iter_f_monotone_restricted(f, m_topo[i], m_children[i], geq_list)
            return sum(1 for _ in it)
        k_independent = [num(k) for k in range(m)]
        return reduce(lambda a,b: a*b, k_independent, 1)
    

    # Methods related with entropy

    @cached_property
    def num_lt(self):
        return self.num_leq - self.n
    @cached_property
    def num_leq(self):
        return self.leq.sum()
    @cached_property
    def num_child(self):
        return self.child.sum()
    @cached_property
    def num_incomparable(self):
        cmp = self.leq | self.leq.T
        return (~cmp).sum() // 2
    
    
    # Testing methods
    
    def _test_iters_diff(self, it1, it2):
        'Compute set1 = set(it1)-set(it2) and set2 = set(it2)-set(it1)'
        'Assumes that the iterators do not repeat elements'
        set1 = set()
        set2 = set()
        for x,y in zip(it1, it2):
            if x != y:
                if x in set2:
                    set2.discard(x)
                else:
                    set1.add(x)
                if y in set1:
                    set1.discard(y)
                else:
                    set2.add(y)
        for x in it1:
            set1.add(x)
        for y in it2:
            set2.add(y)
        return set1, set2
    
    def _test_iters(self, it1, it2):
        'Check if two iterators yield the same values'
        def timed(it, key):
            cnt = total = 0
            t = time.time()
            for i in it:
                total += time.time()-t
                yield i
                t = time.time()
                cnt += 1
            times[key] = total
            count[key] = cnt
        times = {}
        count = {}
        it1 = timed(it1, 0)
        it2 = timed(it2, 1)
        set1, set2 = self._test_iters_diff(it1, it2)
        same = not set1 and not set2
        reason = not same and (
            f'Iterators are different:\n'
            f'Found by 1 not by 2: {set1}\n'
            f'Found by 2 not by 1: {set2}'
        )
        self._test_summary(times, count, same, reason)
    
    def _test_counts(self, f1, f2):
        times = {}
        count = {}
        t = time.time()
        count[0] = f1()
        times[0] = time.time() - t
        t = time.time()
        count[1] = f2()
        times[1] = time.time() - t
        same = count[0]==count[1]
        reason = not same and (
            f'Methods are different:\n'
            f'Output of 1: {count[0]}\n'
            f'Output of 2: {count[1]}'
        )
        self._test_summary(times, count, same, reason)
        
    def _test_summary(self, times, count, same, reason):
        print(
            f'repr: {self}\n'
            f'hash: {hash(self)}\n'
            f'n: {self.n}\n'
            f'is_distributive: {self.is_distributive}\n'
            f'Time used by method 1: {round(times[0], 2)}s\n'
            f'Time used by method 2: {round(times[1], 2)}s\n'
            f'Elements found by method 1: {count[0]}\n'
            f'Elements found by method 2: {count[1]}\n'
            f'Same output: {same}\n'
        )
        if not same:
            self.throw(reason)
            
    def _test_assert_distributive(self):
        self.is_distributive or self.throw(
            f'The test can not be executed because the lattice is not distributive:\n'
            f'{self._distributive_error}'
        )
    
    def test_iter_f_monotone(self, outfile=None):
        it1 = map(tuple, self.iter_f_monotone())
        it2 = map(tuple, self.iter_f_monotone_bruteforce())
        with Outfile(outfile):
            self._test_iters(it1, it2)
    
    def test_iter_f_monotone_bottom(self, outfile=None):
        it1 = map(tuple, self.iter_f_monotone_bottom())
        it2 = map(tuple, self.iter_f_monotone_bottom_bruteforce())
        with Outfile(outfile):
            self._test_iters(it1, it2)
    
    def test_iter_f_lub_bottom(self, outfile=None):
        it1 = map(tuple, self.iter_f_lub_bottom())
        it2 = map(tuple, self.iter_f_lub_bottom_bruteforce())
        with Outfile(outfile):
            self._test_iters(it1, it2)
    
    def test_iter_f_lub(self, outfile=None):
        it1 = map(tuple, self.iter_f_lub())
        it2 = map(tuple, self.iter_f_lub_bruteforce())
        with Outfile(outfile):
            self._test_iters(it1, it2)
    
    def test_iter_f_lub_bottom_distributive(self, outfile=None):
        self._test_assert_distributive()
        it1 = map(tuple, self.iter_f_lub_bottom())
        it2 = map(tuple, self.iter_f_lub_bottom_distributive())
        with Outfile(outfile):
            self._test_iters(it1, it2)
    
    def test_count_f_lub_bottom_distributive(self, outfile=None):
        self._test_assert_distributive()
        f1 = lambda : self.count_f_lub_bottom_distributive()
        f2 = lambda : self.count_f_lub_bottom_bruteforce()
        with Outfile(outfile):
            self._test_counts(f1, f2)
    
    
    # Methods for serialization
    
    def to_literal(self):
        def get_dtype_string(dtype):
            import re
            aux = re.compile(r'(<class \'numpy\.(.*)\'>)|(<class \'(.*?)\'>)|(.*)')
            g = aux.match(str(dtype)).groups()
            s = g[1] or g[3] or g[4]
            assert dtype == np.dtype(s), f'Non invertible dtype: {dtype} != np.dtype(\'{s}\')'
            return s
        
        out = self.__dict__.copy()
        for key, value in out.items():
            if isinstance(value, np.ndarray):
                out[key] = {
                    'dtype': get_dtype_string(value.dtype),
                    'array': value.tolist()
                }
        return out
    
    @classmethod
    def from_literal(cls, lit):
        def read_numpy_array(dict_):
            arr = np.array(dict_['array'], dtype=dict_['dtype'])
            if arr.size == 0:
                arr = arr.reshape((0, 0))
            arr.flags.writeable = False
            return arr
        V = cls(read_numpy_array(lit.pop('leq')))
        for key,value in lit.items():
            if isinstance(value, dict) and 'dtype' in value and 'array' in value:
                value = read_numpy_array(value)
            V.__dict__[key] = value
        return V
    
    # Operations between lattices
    
    def stack_edge(self, other, above=False):
        'Stacks self below other, adding an edge between self.top and other.bottom'
        if above:
            return other.stack_edge(self)
        c = [[i for i in l] for l in self.children]
        c.extend([[self.n+i for i in l] for l in other.children])
        c[self.n+other.bottom].append(self.top)
        return self.__class__.from_children(c)
    
    def series(self, other, above=False):
        'Stacks self below other, replacing self.top with other.bottom'
        if above:
            return other.series(self)
        if other.bottom != 0:
            other = other.canonical
        c = [[i for i in l] for l in self.children]
        c.extend([[self.n-1+i for i in l] for l in other.children][1:])
        return self.__class__.from_children(c)
    
    def decompose_series(self):
        n = self.n
        leq = self.leq
        cmp = leq | leq.T
        nodes = sorted(range(n), key=lambda i: leq[:,i].sum())
        cuts = [i for i in nodes if cmp[i,:].sum()==n]
        subs = [nodes[i:j] for i,j in zip(cuts, cuts[1:])]
        return [self.subgraph(sub) for sub in subs]


def example_2002(cls=Poset):
    grid = [[],[0],[0],[1],[1,2],[2],[3,4],[4,5],[6,7]]
    grid.extend([[0],[0],[9,2],[10,1]])
    for i,j in [(3,9),(5,10),(6,11),(7,12)]:
        grid[i].append(j)
    return cls.from_children(grid)

def example_1990(cls=Poset):
    grid = [[],[0],[0],[1],[1,2],[2],[3,4],[4,5],[6,7]]
    children = [[j+9*(i>=9) for j in grid[i%9]] for i in range(18)]
    for i,j in [(9,4),(10,6),(11,7),(13,8)]:
        children[i].append(j)
    return cls.from_children(children)
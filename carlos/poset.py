from cached_property import cached_property
import pyhash
import numpy as np
from collections import deque
from itertools import product, groupby
from functools import reduce


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
        V = Poset.from_covers([[],[0],[0],[1,2]])
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
        return f'{self.children}'
    
    def show(self, f=None, as_edges=False, save=None):
        'Use graphviz to display or save self or the endomorphism f if given'
        n = self.n
        child = self.child
        extra_edges = None
        labels = range(n)
        if f is not None:
            n = self.n
            if as_edges:
                extra_edges = [(i,f[i]) for i in range(n)]
            else:
                labels = ['' for i in range(n)]
                for i, l in groupby(range(n), f.__getitem__):
                    labels[i] = ','.join(map(str,l))
        
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
    
    def throw(self, message):
        print(message)
        self.show()
        print('Relation matrix:')
        print(self.leq.astype(int))
        raise Exception(message)
     
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
        return cls.from_covers(children)
    
    @classmethod
    def from_covers(cls, children):
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
        poset.__dict__['dist'] = child
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
        dist = child.astype(int)
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
        remaining = product(*matches)
        
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
        for i,j in product(range(n), repeat=2):
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
        for i,j in product(range(n), repeat=2):
            if not fb[i,j] and not (h[i],h[j]) in vis:
                vis.add((h[i],h[j]))
                out = new_leq.copy()
                out[:-1, :-1] += np.matmul(leq[:, i:i+1], leq[j:j+1, :])
                out[n, :-1] = leq[j,:]
                out[:-1, n] = leq[:,i]
                out.flags.writeable = False
                yield self.__class__(out)
        return
    
    
    # Methods for endomorphisms (functions self->self)
    
    def iter_f_all(self):
        'all endomorphisms'
        return product(range(self.n), repeat=self.n)
    
    @cached_property
    def num_f_all(self):
        return self.n**self.n
        
    def iter_f_all_bottom(self):
        'all endomorphisms f with f(bottom)=bottom'
        n = self.n
        if n>0:
            options = [range(n) if i!=self.bottom else [i] for i in range(n)]
            for f in product(*options):
                yield f
        return
    
    @cached_property
    def num_f_all_bottom(self):
        return self.n**(self.n-1)
    
    def f_is_lub(self, f, domain=None):
        'check if f preserves lub and f[bottom] is bottom. Throws if no bottom'
        n = self.n
        if n==0 or (domain is not None and len(domain)<=1):
            return True
        bot = self.bottom
        if f[bot] != bot or (domain is not None and bot not in domain):
            return False
        domain = range(n) if domain is None else domain
        lub = self.lub
        for i in domain:
            for j in domain:
                if f[lub[i,j]]!=lub[f[i],f[j]]:
                    return False
        return True

    def iter_f_lub_bruteforce(self):
        'all join endomorphisms. Throws if no bottom.'
        for f in self.iter_f_all_bottom():
            if self.f_is_lub(f):
                yield f
        return
    
    def count_f_lub_bruteforce(self):
        return sum(1 for f in self.iter_f_lub_bruteforce())

    def f_is_monotone(self, f, domain=None):
        n = self.n
        domain = range(n) if domain is None else domain
        leq = self.leq
        for i in domain:
            for j in domain:
                if leq[i,j] and not leq[f[i],f[j]]:
                    return False
        return True
    
    def iter_f_monotone_restricted(self, domain=None, _f=None):
        'generate all monotone functions f : domain -> self, padding non-domain with None'
        D = range(n) if domain is None else domain
        n = self.n
        sub = self.subgraph(D)
        m = sub.n
        cov = sub.children
        inv = [None for i in range(n)]
        for i in range(m):
            inv[D[i]] = i
        topo = [inv[i] for i in self.toposort if i in D]
        leq = self.leq
        geq_list = [[j for j in range(n) if leq[i,j]] for i in range(n)]
        lub = self.lub
        f = [None for i in range(n)] if _f is None else _f
        def backtrack(i):
            'f[D[topo[j]]] is fixed for all j<i. Backtrack f[D[topo[k]]] for all k>=i, k<m'
            if i==m:
                yield f
            else:
                if not cov[topo[i]]:
                    options = range(n)
                else:
                    f_cov = (f[D[j]] for j in cov[topo[i]])
                    cov_lub = reduce((lambda a,b: lub[a,b]), f_cov)
                    options = geq_list[cov_lub]
                for k in options:
                    f[D[topo[i]]] = k
                    yield from backtrack(i+1)
        yield from backtrack(0)
    
    
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
                msg = f'Non distributive lattice: ' +\
                f'{i} glb ({j} lub {k}) = {i} glb {lub[j,k]} = ' +\
                f'{glb[i,lub[j,k]]} != {lub[glb[i,j],glb[i,k]]} = ' +\
                f'{glb[i,j]} lub {glb[i,k]} = ({i} glb {j}) lub ({i} glb {k})'
                return msg
        return None
    
    def assert_distributive(self):
        self.is_distributive or self.throw(self._distributive_error)
    
    def count_f_lub_distributive_v1(self):
        return sum(1 for f in self.iter_f_lub_distributive_v1())
    
    def iter_f_lub_distributive_v1(self):
        'generate and interpolate all monotone functions over irreducibles'
        self.assert_distributive()
        if not self.n:
            return
        irr = self.irreducibles
        funcs = self.iter_f_monotone_restricted(irr)
        yield from self._interpolate_funcs(funcs, irr)
    
    def _interpolate_funcs(self, funcs, domain):
        'extend each f in funcs outside domain using f[j]=lub(i if i<=j and i in domain)'
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
    
    def iter_f_lub_distributive_v2(self):
        'same as _v1 but treating different irreducible components separately'
        if not self.n:
            return
        f, aux = self._aux_iter_f_lub_distributive_v2()
        funcs = product(*aux)
        yield from self._interpolate_funcs(irr, funcs)
        
    def _aux_iter_f_lub_distributive_v2(self):
        'treat irreducibles components separately'
        self.assert_distributive()
        n = self.n
        irr = self.irreducibles
        sub = self.subgraph(irr)
        subcomps = sub.independent_components
        irrcomps = [[irr[i] for i in comp] for comp in subcomps]
        f = [None for i in range(n)]
        partial = lambda D: self.iter_f_monotone_restricted(D, _f=f)
        return f, [partial(domain) for domain in irrcomps]
    
    def count_f_lub_distributive_v2(self):
        if not self.n:
            return 0
        f, partials = self._aux_iter_f_lub_distributive_v2()
        independent = [sum(1 for _ in gen) for gen in partials]
        return reduce(lambda a,b: a*b, independent, 1)
    
    @cached_property
    def num_f_lub(self):
        if self.n >= 4 and self.is_distributive:
            num = self.count_f_lub_distributive_v2()
        else:
            num = self.count_f_lub_bruteforce()
        return num

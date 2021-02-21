import sqlite3
from ast import literal_eval


class KVStorage(dict):
    '''
    Holds a dictionary in disk and in memory simultaneously.
        Keys must be integers.
        Values must support literalization with to_literal and from_literal.
    '''
    
    def __init__(self, file, to_literal, from_literal):
        self.conn = sqlite3.connect(file)
        try:
            it = self.conn.execute('SELECT * FROM storage')
        except sqlite3.OperationalError:
            self.conn.execute('''
            CREATE TABLE storage(
                key INTEGER NOT NULL PRIMARY KEY,
                value TEXT
            )
            ''')
            it = iter(tuple())
        super().update({k: from_literal(literal_eval(v)) for k,v in it})
        self.to_literal = to_literal
        self.from_literal = from_literal
    
    def __setitem__(self, key, value):
        assert isinstance(key, int), f'key {key} must be an integer'
        lit = self.to_literal(value)
        lit_str = str(lit)
        assert literal_eval(lit_str) == lit, f'literal_eval("{lit_str}") != {lit}'
        assert self.from_literal(lit) == value, f'from_literal({lit}) != {value}'
        self.conn.execute('REPLACE INTO storage VALUES (?,?)', (key, lit_str))
        self.conn.commit()
        super().__setitem__(key, value)


class PosetStorage():
    
    '''
    Set of posets with replace (add or replace) operation and sampling
    The cls parameter is the Poset class that can parse literals
    '''
    
    def __init__(self, cls, file=':memory:'):
        to_literal = lambda posets: [V.to_literal() for V in posets]
        from_literal = lambda lits: [cls.from_literal(lit) for lit in lits]
        self.h = KVStorage(file, to_literal, from_literal)
        self.file = file
        self.n = sum(len(v) for v in self.h.values())
        self.pop = Population(*self.h.keys())
    
    def __repr__(self):
        n = self.n
        n_elems = '1 element' if n==1 else f'{n} elements'
        return f'PosetStorage(file="{self.file}") containing {n_elems}'
    
    def add(self, V):
        'adds the lattice V to the collection. Combines cached properties if already present'
        posets = self.h.get(hash(V), [])
        U = next((U for U in posets if U==V), None)
        if U is not None:
            U.__dict__.update(V.__dict__)
        else:
            V.assert_lattice()
            posets.append(V)
            self.n += 1
        self.h[hash(V)] = posets
        self.pop.add(hash(V))
        return
    
    def remove(self, V):
        'remove the lattice V from the collection'
        self.h.get((hash(V)), []).remove(V)
        self.pop.remove(hash(V))
    
    def __contains__(self, V):
        return next((U for U in self.h.get(hash(V), tuple()) if U==V), None)
    
    def __iter__(self):
        yield from (U for posets in self.h.values() for U in posets)
    
    def __len__(self):
        return self.n
    
    def sample(self, k=None, filter=None):
        it = (V for V in self.iter_random(filter))
        out = [V for V,_ in zip(it, range(1 if k is None else k))]
        return out[0] if k is None else out
    
    def iter_random(self, filter=None):
        it = (random.choice(self.h[key]) for key in self.pop.iter_random())
        if filter is None:
            yield from it
        else:
            yield from (V for V in it if filter(V))
        
import random

class Population():
    '''
    Population that can grow, shrink and be sampled quickly.
    Elements must be unique and hashable
    '''
    
    def __init__(self, *elems):
        self.l = list(elems)
        self.idx = {e:i for i,e in enumerate(elems)}
    
    def __repr__(self):
        s = repr(self.l)
        return f'Pop({s[1:-1]})'
    
    def __len__(self):
        return len(self.l)
    
    def __iter__(self):
        return iter(self.l)
    
    def __contains__(self, elem):
        return elem in self.idx
    
    def append(self, elem):
        assert elem not in self.idx, f'Can not add multiple instances of {elem}'
        self.idx[elem] = len(self.l)
        self.l.append(elem)
    
    def add(self, elem):
        if elem not in self:
            self.append(elem)
    
    def pop(self, i=None):
        'pop a random element, or that in the given position'
        if i is None:
            i = random.choice(range(len(self)))
        self._swap(i, len(l)-1)
        elem = self.l.pop()
        del self.idx[elem]
        return elem
    
    def _swap(self, i, j):
        'swap the ith element with the last element'
        self.l[j], self.l[i] = ei, ej = self.l[i], self.l[j]
        self.idx[ei] = j
        self.idx[ej] = i
    
    def remove(self, elem):
        return self.pop(self.idx[elem])
    
    def sample(self, k=None):
        if k is None:
            return self.sample(1)[0]
        idx = random.sample(range(len(self)), k)
        return [self.l[i] for i in idx]
    
    def iter_random(self):
        for n in range(len(self), 0, -1):
            i = random.choice(range(n))
            self._swap(i, n-1)
            yield self.l[n-1]


def as_html(file, posets):

    def as_table(V):
        svg1 = V.graphviz().create_svg().decode()
        svg2 = V.subgraph(V.irreducibles).graphviz(labels=V.irreducibles).create_svg().decode()
        return f'''<table><tr>
        <td>{svg1}</td>
        <td>{svg2}
        <pre># f: {V.num_f_lub_bottom}</pre>
        <pre>series decomp: {len(V.decompose_series())}</pre>
        </td>
        </tr></table>'''

    tables = '\n'.join(as_table(V) for V in posets)

    content=f'''
    <html><head><style>
    table tr td{{ padding:0.3em; page-break-inside: avoid; border: 1px solid}}
    table{{border-collapse: collapse}}
    table{{zoom: 0.5; display:inline-block; margin: 0.3em; text-align: justify;}}
    body {{
      display: flex;
      flex-wrap: wrap;
      justify-content: space-around;
      align-items: center;
    }}
    @media print {{
        table tbody tr td:before,
        table tbody tr td:after {{
            content : "" ;
            height : 4px ;
            display : block ;
        }}
    }}
    </style></head>
    <body>{tables}</body></html>
    '''
    with open(file, 'w') as f:
        f.write(content)
    return

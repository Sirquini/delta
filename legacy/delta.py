import random

from sys import stdin
from itertools import product
from itertools import combinations
#from graphviz import Graph
from time import time


# Matrix of LUBs
def lubs():
	global ml
	ml = [[0 for i in range(C)] for i in range(C)]
	for i in range(C):
		for j in range(C):
			ml[i][j] = pairlub(i,j)

# Matrix of GLBs
def glbs():
	global mg
	mg = [[0 for i in range(C)] for i in range(C)]
	for i in range(C):
		for j in range(C):
			mg[i][j] = pairglb(i,j)

# Functions ub(a,b), lb(a,b), minl(v) and maxl(v) are used to 
# calculate least upper bounds and greatest lower bounds for pairs (a,b)
# such that neither a < b nor b > a.

# Upper bounds of {a, b}
def ub(a,b):
	ub = []
	for i in range(C):
		if l[i][a] == 1 and l[i][b] == 1: ub.append(i)
	return ub

# Lower bounds of {a, b}
def lb(a,b):
	lb = []
	for i in range(C):
		if l[a][i] == 1 and l[b][i] == 1: lb.append(i)
	return lb

# Minimum of list v
def minl(v):
	A = set()
	for r in v:
		t = True
		for s in v:
			if t and l[s][r] == 0: t = False
		if t: A = A | {r}
	assert (len(A) == 1)
	return A.pop()

# Maximum of list v
def maxl(v):
	A = set()
	for r in v:
		t = True
		for s in v:
			if t and l[r][s] == 0: t = False
		if t: A = A | {r}
	assert (len(A) == 1)
	return A.pop()

# Least upper bound of {a, b}
def pairlub(a,b):
	c = 0
	if l[a][b] == 1: c = a
	elif l[b][a] == 1: c = b
	else: c = minl(ub(a,b))
	return c

# Greatest lower bound of {a, b}
def pairglb(a,b):
	c = 0
	if l[a][b] == 1: c = b
	elif l[b][a] == 1: c = a
	else: c = maxl(lb(a,b))
	return c

# Least upper bound of list t
def lub(t):
	r = 0
	for i in t:
		r = ml[r][i]
	return r

# Greatest lower bound of list t
def glb(t):
	r = 0
	for i in t:
		r = mg[r][i]
	return r

# Order in functions
# Returns true if f1 < f2. Otherwise false
def leq(f1,f2):
	i = 0
	t = True
	while i < len(f1) and t:
		if l[f2[i]][f1[i]] == 1: i+=1
		else: t = False
	return t

# Max function of list F
def maxf(F):
	A = set()
	for r in F:
		t = True
		for s in F:
			if t and not leq(s,r): t = False
		if t: A = A | {tuple(r)}
	assert (len(A) == 1)
	return A.pop()

# Returns true if l is a lattice.
def IsLattice(l):
	# List of pairs
	lt = list(combinations(range(C),2))
	test = True
	i = 0
	while test and i < len(lt):
		if lub([lt[i][0],lt[i][1]]) == -1: test = False
		else: i+=1
	return test

# Generating functions that preserve bottom
def gen_func():
	F = []
	A = []
	F.append([0])
	j = C-1
	while j > 0:
		while len(F) > 0:
			a = F.pop()
			b = a.copy()
			for i in range(C):
				a.append(i)
				A.append(a)
				a = b.copy()
		j-=1
		F = A.copy()
		A.clear()
	return F

# Filtering distributive functions
def dist_f(sf):
	# List of pairs
	t = list(combinations(range(C),2))
	d = []
	for f in sf:
		test = True
		i = 0
		while i < len(t) and test:
			t0 = t[i][0]
			t1 = t[i][1]
			if not( f[lub([t0,t1])] == lub([f[t0],f[t1]]) ): test = False
			i+=1
		if test: d.append(f)
	return d

# Calculating lubs: all combinations of n elements
# L stores tuples and c the respective least upper bound, i.e.,
# L[i] is such Lub(L[i]) = c[i]
def c_lubs(n):
	tuples = list(product(range(C),repeat=n))
	L = []
	c = []
	for t in tuples:
		L.append(t)
		c.append(lub(list(t)))
	return L,c

# Delta* for a set of functions F
def delta_ast(F):
	n = len(F)
	# All possibles tuples such that for each t in tuples and c belonging to the lattice, lub(t) > c
	tuples = list(product(range(C),repeat=n))
	delta = []
	for c in range(C):
		s = []
		for t in tuples:
			# lub(t) > c
			if l[lub(list(t))][c] == 1:
				sl = 0
				i = 0
				while i < n:
					# e represents the value of the spatial function F[i] at t[i]
					e = F[i][t[i]]
					if l[e][sl] == 1: sl = e
					elif l[sl][e] == 1: sl = sl
					else: sl = lub([sl,e])
					i+=1
				s.append(sl)
		if 0 in s: delta.append(0)
		else: delta.append(glb(s))
	return delta

# Delta* for a set of functions F (precalculating lubs)
def delta_ast2(F,L,w):
	# List of all pairs of elements in the lattice
	n = len(F)
	m = len(w)
	delta = []
	for c in range(C):
		s = C-1
		forms_derive_c = 0
		top_in_tuple = 0
		for u in range(m):
			if l[w[u]][c] == 1:
				forms_derive_c +=1
				e = [F[i][L[u][i]] for i in range(n)]
				if (C-1) in e:
					top_in_tuple+=1
					j = C-1
				else: j = lub(e)
				if l[s][j]==1: s = j
				elif l[j][s]==1: s = s
				else: s = glb([s,j])
		delta.append(s)
		print('#---- c:',c, 'Forms to derive c: ', forms_derive_c)
		print('#---- Top in tuple: ', top_in_tuple)
	return delta

# Delta for a set of functions S
def delta_n(F, S):
	cand = []
	num_fun_under = 0
	for f in F:
		t = True
		j = 0
		while j < len(S) and t:
			if not leq(f,S[j]): t = False
			else: j+=1
		if t:
			cand.append(f)
			num_fun_under+=1
	print('#---- delta_n: num_fun_under ',num_fun_under)
	return maxf(cand)


#####################################################################################
########################## Complementary functions ##################################
#####################################################################################

# Print space functions
def print_func(sf):
	for j in range(len(sf)):
		print ('Function number: ', j)
		print(sf[j])
		print('----------')
	return 0	

# Transitive reduction of lattice l
def reduction(l):
	nl = list(combinations(range(C),2))
	j = 0
	while j < len(nl):
		s = nl[j][0]
		t = nl[j][1]
		test = True
		if l[t][s] == 0:
			test = False
			nl.remove((s,t))
		else:
			v = s + 1
			while test and v < t: 
				if l[v][s] == 1 and l[t][v] == 1:
					nl.remove((s,t))
					test = False
				v += 1
		if test: j+=1
	return nl

# Generates a pdf with a graphic representation of the transitive closure of the lattice
def print_lattice(l):
	g = Graph()
	for i in range(C): g.node(str(i))
	for i in range(C):
		for j in range(C):
			if i != j and l[i][j] == 1:
				g.edge(str(i),str(j))
	#print(g.source)
	g.render('lattice', view=True)
	return 0
# Generates a pdf with a graphic representation of the lattice
def print_lattice_red(l):
	g = Graph()
	for i in range(C): g.node(str(i))
	for i in l:
		g.edge(str(i[1]),str(i[0]))
	#print(g.source)
	g.render('lattice_red', view=True)
	return 0

#####################################################################################
#####################################################################################

def main():
	s_time = time()

	# Elements in the lattice
	global C
	C = 6

	# Lattice
	global l

	l = [[0 for i in range(C)] for i in range(C)]
	for i in range(C):
		l[i][0] = 1
		l[i][i] = 1
		l[C-1][i] = 1

	# Matrices of least upper bounds and greatest lower bounds
	lubs()
	glbs()

	# Calculating space functions
	h = gen_func()
	h1 = time()
	print('Bot preserving funtions: ', len(h))
	sf = dist_f(h)
	
	aux_t = time()-s_time
	print('Space functions preprocessing time: ',aux_t)
	print('Space funtions: ', len(sf))
	print('________________________________________________________________\n')


	## Code to calculate delta and delta* for a given list n of space funtions (agents)
	## 'print_func' generates space functions and their identifiers 

	# n = [ int(x) for x in stdin.readline().split()]
	# m = len(n)
	# while m > 0:
	# 	h1 = time()
	# 	L,c = c_lubs(m)
	# 	at = time()-h1
	# 	print('Preprocessing LUBs: ', at)
	# 	I = [sf[i] for i in n]
	# 	print('\n__________SPACE FUNCTIONS__________\n')
	# 	for i in range(m):
	# 		print('Id: ', n[i], '|', I[i])
	# 	print('___________________________________\n')
	# 	h1 = time()
	# 	print('delta*: ' + repr(delta_ast(I)))
	# 	h1 = time()-h1
	# 	print('-- Delta* time: ',h1, '\n')
	# 	h1 = time()
	# 	print('delta*(PRE): ' + repr(delta_ast2(I,L,c)))
	# 	h1 = time()-h1
	# 	print('-- With preprocessing: ', h1)
	# 	print('-- Without preprocessing: ', at + h1, '\n')
	# 	h1 = time()
	# 	print('delta: ' + repr(delta_n(sf,I)))
	# 	h1 = time()-h1
	# 	print('-- With preprocessing: ',h1)
	# 	print('-- Without preprocessing: ', aux_t + h1)
	# 	print('------------------------------------------------------------------\n')
	# 	n = [ int(x) for x in stdin.readline().split()]
	# 	m = len(n)
	
	## Code to calculate delta and delta* for a random list n of space funtions (agents) len(n) = m

	m = 1
	while m < 8:
		n = random.sample(range(len(sf)), m)
		h1 = time()
		L,c = c_lubs(m)
		at = time()-h1
		print('Preprocessing LUBs time: ', at)
		I = [sf[i] for i in n]
		print('\n__________SPACE FUNCTIONS__________\n')
		for i in range(m):
			print('Id: ', n[i], '|', I[i])
		print('___________________________________\n')
		h1 = time()
		print('\ndelta*: ' + repr(delta_ast(I)))
		h1 = time()-h1
		print('-- Delta* time: ',h1, '\n')
		h1 = time()
		print('\ndelta*(PRE): ' + repr(delta_ast2(I,L,c)))
		h1 = time()-h1
		print('-- With preprocessing: ', h1)
		print('-- Without preprocessing: ', at + h1, '\n')
		h1 = time()
		print('\ndelta: ' + repr(delta_n(sf,I)))
		h1 = time()-h1
		print('-- With preprocessing: ',h1)
		print('-- Without preprocessing: ', aux_t + h1)
		print('------------------------------------------------------------------\n')
		m +=1

	# print(IsLattice(l))
	# print_lattice_red(reduction(l))
	# print_lattice(l)
	# print_func(sf)

	e_time = time() - s_time
	print('Total time: ', e_time)
main()

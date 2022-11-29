# -*- coding: utf-8 -*-
"""
Created on Sun Aug 21 12:50:51 2022

@author: jogib
"""

import networkx as nx
import itertools
import numpy as np
import matplotlib.pyplot as plt
import copy
import functools

#%%
pauli = {
    "I": np.array(np.mat("1+0j,0;0,1")),
    "Z": np.array(np.mat("1+0j,0;0,-1")),
    "X": np.array(np.mat("0+0j,1;1,0")),
    "Y": np.array(np.mat("0,-1j;1j,0")),
}
#%%
def theoretical_size(m,n):
    """
    Parameters
    ----------
    m : int
        cardinality, (size of clique)
    n : int
        length of pauli string.
    """
    
    if not m%2:
        return 0
    else:
        s = lambda n,k: (4**n)/(2**k) - (1-k%2) 
        return int(1/np.math.factorial(m)*np.prod([s(n,k) for k in range(0,m-1)]))
#%%
def anticommute_check(A: str, B: str):
    """
    Checks if two strings anti commute
    
    TODO: refactor XY calcu two bit string quiskit does it
    
    Parameters
    ----------
    A : str
    B : str
    Returns
    -------
    bool

    """
    ch_arr = []
    for i in range(len(A)):
        a_t = pauli[A[i]]
        b_t = pauli[B[i]]
        check = np.dot(a_t, b_t) + np.dot(b_t, a_t)
        ch_arr.append(not np.any(check))
    return np.count_nonzero(ch_arr) % 2

def gen_graph(n:int =2):
    """
    Example
    ----------
    G = gen_graph(n)
    G.remove_node("".join("I" for i in range(n)))

    Parameters
    ----------
    n : int, optional
        String length. The default is 2.

    Returns
    -------
    G : networkx graph
    """
    lis = ["X", "Y", "Z", "I"]
    combos = itertools.product(lis, repeat=n)
    d = {i: "".join(j for j in x) for i, x in enumerate(combos)}
    G = nx.empty_graph(len(d.values()))
    G = nx.relabel_nodes(G, d)
    # add edges
    for i in itertools.combinations(list(d.values()), 2):

        A = i[0]
        B = i[1]
        if anticommute_check(A, B):
            G.add_edge(i[0], i[1], **{"color": "tab:blue", "width": 0.2})
    return G

def gen_paulis(lis):
    """
    Parameters
    ----------
    lis : list
        list of the form ['XYZ','XXI'...] so that each the first char in
        each element of the list is pauli corresponding to the char.

    Returns
    -------
    paul_lis : lis
        [[pauli(x),pauli(y) pauli(z)],].

    """
    pl = [[pauli[j] for j in i] for i in list(lis)]
    pauli_lis = [[row[i] for row in pl] for i in range(len(list(lis)[0]))]
    return pauli_lis

def multiply_lis(lis):
    list_prod = [functools.reduce(np.dot, i) for i in lis]
    return list_prod

def find_pauli(lis:np.array):
    """
    takes a list of numerical paulies and returns a list of pauli labels
    Parameters
    ----------
    lis : np.array
        Must be a list of 2x2 arrays each corresponding to a pauli matrix

    Returns
    -------
    st : TYPE
        list of pauli strings
    """

    st = ""
    for i in lis:
        if i[0][0] != 0:
            if i[0][0] + i[1][1] == 0:
                st = st + "Z"
            else:
                st = st + "I"
        elif i[0][1] + i[1][0] == 0:
            st = st + "Y"
        else:
            st = st + "X"
    return st

def order_of_repeats(cliques):
    """
    TODO: not complete need to change how it does the counting

    Parameters
    ----------
    cliques : TYPE
        DESCRIPTION.

    Returns
    -------
    total : TYPE
        DESCRIPTION.

    """
    total = {}
    for element in cliques:
        each = {}
        for st in element:
            count = {}
            for s in st:
              if s in count:
                count[s] += 1
                each[s] += 1
              else:
                count[s] = 1
                each[s] = 1
        order = max(list(each.values()))
        if order in total:
            total[order]+=1
        else:
            total[order]=1
    return total,each,count

def add_triang(cliques):
    """
    TODO missing the odd number correction
    
    Find the triangles for the group of triangles for one stringlength greater
    than the original
    
    
    Example
    ----------
        n = 3
        G = gen_graph(n)
        G.remove_node("".join("I" for i in range(n)))
        cliques = []
        for i in nx.find_cliques(G):
            if len(i) == 3:
                cliques.append(set(i))
        add_triang(cliques)
        
    Parameters
    ----------
    cliques : list of anticomuting pairs

    Returns
    -------
    new_triang : TYPE
        returns thee triangles for one string length greater.



    todo: need a better way to generate things you can add
    the weird odd additions
    """
    new_triang=set()
    n=len(list(cliques)[0])
    for tria in cliques:
        trian= list(tria)
        for i in range(n+1):
            st=frozenset([x[:i]+'I'+x[i:] for x in trian])
            new_triang.add(st)
            for p in {'X','Y','Z'}:
                one = trian[0][:i]+'I'+trian[0][i:]
                two = trian[1][:i]+p+trian[1][i:]
                three=trian[2][:i]+p+trian[2][i:]
                st=frozenset([one,two,three])
                new_triang.add(st)
                one = trian[0][:i]+p+trian[0][i:]
                two = trian[1][:i]+'I'+trian[1][i:]
                three=trian[2][:i]+p+trian[2][i:]
                st=frozenset([one,two,three])
                new_triang.add(st)
                one = trian[0][:i]+p+trian[0][i:]
                two = trian[1][:i]+p+trian[1][i:]
                three=trian[2][:i]+'I'+trian[2][i:]
                st=frozenset([one,two,three])
                new_triang.add(st)
    return new_triang
#%%
n = 3
m=7
G = gen_graph(n)
G.remove_node("".join("I" for i in range(n)))
print(nx.algorithms.approximation.large_clique_size(G))
# print(clique_counter(G, 3))
cliques = []
for i in nx.find_cliques(G):
    if len(i) == m:  # max_xs:
        cliques.append(set(i))
print(theoretical_size(m,n),len(cliques))
#%%
n = 2
G = gen_graph(n)
G.remove_node("".join("I" for i in range(n)))
print(nx.algorithms.approximation.large_clique_size(G))
# print(clique_counter(G, 3))
#%%
cliques = []
for i in nx.find_cliques(G):
    if len(i) == 3:  # max_xs:
        cliques.append(set(i))
#%%
n = 4
m= 3
G = gen_graph(n)
G.remove_node("".join("I" for i in range(n)))
cliques = []
for i in nx.find_cliques(G):
    if len(i) == m:
        cliques.append(set(i))
print(len(cliques),theoretical_size(m,n))
trigs = add_triang(cliques)
print(len(trigs),theoretical_size(m,n+1))
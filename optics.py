#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""A clustering algorithm 
from Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Jörg Sander (1999). 
"OPTICS: Ordering Points To  Identify the Clustering Structure".
ACM SIGMOD international conference on Management of data. 
ACM Press. pp. 49–60."""

import heapq
import math
import numpy as np
import itertools 

def optics(points,eps,M=15):
    """ Optics algorithm
    :param points:
    :param eps:
    :param M:
    """
    nb_points = len(points)
    o = 0
    order = [0] * nb_points
    reach = [0] *  nb_points
    processed = [0] * nb_points
    

    for p in range(nb_points):
        if not processed[p]:
            N = get_neighbors(p, eps,points)
            processed[p] = 1
            order[o]=p
            o += 1
            seeds = PriorityQueue()
            if core_dist(p,N,M,points):
                update(N,p,seeds,eps,M,points,processed,reach)
                empty = False
                while not empty:
                    try :
                        q = seeds.pop()
                    except KeyError:
                        empty = True
                    else:
                        N2 = get_neighbors(p,eps,points)
                        processed[q] = 1
                        order[o]=q
                        o += 1   
                        if core_dist(q,N2,M,points):
                            update(N2,q,seeds,eps,M,points,processed,reach)
                    
    reach = [reach[o] for o in order]
    return order,reach

def core_dist(p,N,M,points):
    """Return the list of the distance to the M nearest points"""
    if len(N) < M:
        return None
    else:
        dist = [d(points[p],points[n]) for n in N]
        dist.sort()
        return dist[M-1]

def d(a,b):
    return math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

def update(N,p,seeds,eps,M,points,processed,reach):
    cd = core_dist(p,N,M,points)
    for n in N:
        if not processed[n]:
            new_rd = max(cd,d(points[p],points[n]))
            if not reach[n]: #n is not in the seed
                reach[n] = new_rd
                seeds.add(n,new_rd)
            else:
                if new_rd < reach[n]:
                    reach[n] = new_rd
                    seeds.add(n,new_rd)
                    

def get_neighbors(p,eps,points):
    ## O(n²) !!!!
    neighbors = []
    for k,q in enumerate(points):
        if d(points[p],q) < eps:
            neighbors.append(k)
    return neighbors

class PriorityQueue(object):
    """A priority queue inspired by http://docs.python.org/2/library/heapq.html"""
    def __init__(self):
        self.pq = []                         # list of entries arranged in a heap
        self.entry_finder = {}               # mapping of items to entries
        self.REMOVED = '<removed-item>'      # placeholder for a removed item
        self.counter = itertools.count()     # unique sequence count

    def add(self,item, priority=0):
        'Add a new item or update the priority of an existing item'
        if item in self.entry_finder:
            self.remove(item)
        count = next(self.counter)
        entry = [priority, count, item]
        self.entry_finder[item] = entry
        heapq.heappush(self.pq, entry)

    def remove(self,item):
        'Mark an existing item as REMOVED.  Raise KeyError if not found.'
        entry = self.entry_finder.pop(item)
        entry[-1] = self.REMOVED

    def pop(self):
        'Remove and return the lowest priority item. Raise KeyError if empty.'
        while self.pq:
            priority, count, item = heapq.heappop(self.pq)
            if item is not self.REMOVED:
                del self.entry_finder[item]
                return item
        raise KeyError('pop from an empty priority queue')

    def get(self):
        """Get last entry or None"""
        try:
            return self.pop()
        except KeyError:
            return None

## Test
points = np.array([[ 15.,  70.],
                  [ 31.,  87.],
                  [ 45.,  32.],
                  [  5.,   8.],
                  [ 73.,   9.],
                  [ 32.,  83.],
                  [ 26.,  50.],
                  [  7.,  31.],
                  [ 43.,  97.],
                  [ 97.,   9.]])
result = [0, 1, 5, 6, 2, 7, 8, 3, 4, 9]
    
def find_cluster(reach,ksi=0.005,M=15):
    steep_down_list = []
    clusters = []
    mib = 0
    i = 0
    nc = 0
    while i < len(reach)-2:
        #print("i={}, reach={}".format(i,reach[i]))
        mib = max(mib,reach[i])
        if reach[i]*(1-ksi) <= reach[i+1]:
            print("STEEP DOWN starting in {}".format(i))
            steep_down_list = update_and_filter(steep_down_list,reach,i,ksi,mib)
            sigma = 0
            while reach[i]<=reach[i-1] and sigma < M:
                i += 1
                if reach[i]*(1-ksi) <= reach[i+1]:
                    sigma = 0
                else:
                    sigma += 1
            print("Ending in {}\n".format(i))
            i +=1
            steep_down_list.append([i,0])
            mib = reach[i]
        elif reach[i]*(1-ksi) >= reach[i+1]:
            print("STEEP UP starting in {}".format(i))
            steep_down_list = update_and_filter(steep_down_list,reach,i,ksi,mib)
            sigma = 0
            while reach[i]>=reach[i-1] and sigma < M:
                i += 1
                if reach[i]*(1-ksi) >= reach[i+1]:
                    sigma = 0
                else:
                    sigma += 1
            print("Ending in {}\n".format(i))
            i +=1
            mib = reach[i]
            for D in steep_down_list:
                if reach[i]*(1-ksi) >= D[1]:
                    s,e = compute_cluster(i,reach,D,ksi)
                    print("Cluster possible in [{},{}] ({})\n".format(s,e,nc))
                    nc += 1
                    if cluster_cond(s,e,reach,M,ksi):
                        clusters.append((s,e))
                        print("Cluster confirmed\n")
        else:
            i += 1
    return clusters

def update_and_filter(steep_down,r,i,ksi,global_mib):
    sd = []
    for (s,mib) in steep_down:
        mib = max(mib,global_mib)
        if r[s] * (1-ksi) < global_mib:
            sd.append([s,mib])
    return sd

def compute_cluster(i,r,D,ksi):
    if r[D[0]] * (1-ksi) >= r[i+1]:
        s = 0 #faux
        return(s,i)
    elif r[i+1] * (1-ksi) >= r[D[0]]:
        e = 0 #faux
        return(D[0],e)
    else:
        return(D[0],i)

def cluster_cond(s,e,r,M,ksi):
    if e-s < M:
        return False
    elif not r[e] <= r[e+1] * (1-ksi):
        return False
    elif not r[s] >= r[s+1] * (1-ksi):
        return False
    else:
        return True


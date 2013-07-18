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
 
import geometry as geo


def optics(points,eps,M=15):
    """ Optics algorithm
    :param points:
    :param eps:
    :param M:
    """
    nb_points = len(points)
    o = 0
    order = [0] * nb_points
    reach = [None] *  nb_points
    processed = [0] * nb_points
    

    for p in range(nb_points):
       
        # for each unprocessed point p of DB
        if not processed[p]:
            # Get p neighbors (eps-radius)
            N = geo.get_neighbors(p, eps,points)
        
            # mark p as processed and output it to the ordered list
            processed[p] = 1
            order[o]=p
            o += 1

            # Seeds = empty priority queue
            seeds = PriorityQueue()

            # if (core-distance(p, eps, Minpts) != UNDEFINED)            
            if core_dist(p,N,M,points) != None:
                update(N,p,seeds,eps,M,points,processed,reach)
                empty = False
                while not empty:
                    try :
                        q = seeds.pop()
                    except KeyError:
                        empty = True
                    else:                     
                        # Get q neighbors (eps-radius)
                        N2 = geo.get_neighbors(q,eps,points)

                        # mark q as processed and output it to the ordered list
                        processed[q] = 1
                        order[o]=q
                        o += 1  

                        #if (core-distance(q, eps, Minpts) != UNDEFINED)
                        if core_dist(q,N2,M,points) != None:
                            update(N2,q,seeds,eps,M,points,processed,reach)

    reach = [reach[o] if reach[o] != None else 0 for o in order]
    return order,reach

def core_dist(p,N,M,points):
    """Return the list of the distance to the M nearest points"""
    if len(N) < M:
        return None
    else:
        dist = [geo.d(points[p],points[n]) for n in N]
        dist.sort()
        return dist[M-1]


def update(N,p,seeds,eps,M,points,processed,reach):
    cd = core_dist(p,N,M,points)
    for n in N:
        if not processed[n]:
            new_rd = max(cd,geo.d(points[p],points[n]))
            # If n is not in the priority queue
            if reach[n] == None:
                reach[n] = new_rd
                seeds.add(n,new_rd)
            else:
                #If we can improve its reachability distance
                if new_rd < reach[n]:
                    reach[n] = new_rd
                    seeds.add(n,new_rd) #update
                    


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

def find_cluster(reach,ksi=0.001,M=5):
    color = ["k"] * len(reach)
    steep_down_list = []
    clusters = []
    mib = 0
    i = 0
    nc = 0
    while i < len(reach)-2:
        # Global mib-value containing the maximum between the end of the last
        # steep and region and the current index.
        mib = max(mib,reach[i])

        ## If i is the start of a steep down area.
        if reach[i]*(1-ksi) > reach[i+1]:
            #print("STEEP DOWN starting in {}".format(i) +
            #      " ({} are opened)".format(len(steep_down_list)))
            steep_down_list = update_and_filter(steep_down_list,reach,i,ksi,mib)

            sD = i   # start of current downward region.
            color[i] = "r"
            
            # Go to the end of the steep area.
            sigma = 0
            while (reach[i+1]<=reach[i] 
                   and sigma < M
                   and i+1<len(reach)-2
                   and not reach[i+1] < reach[i+2]*(1-ksi)):
                # Test if there is no more than M consecutives
                # points that are not ksi-steep.
                if reach[i+1]*(1-ksi) >= reach[i+2]: # i+1 is ksi-steep downward.
                    sigma = 0
                else:
                    sigma += 1
                i += 1
                color[i] = "r"
            #print("Ending in {}\n".format(i))

            eD = i   # end of current downward region.
            #Save the current downward area : [start,end,mib]
            steep_down_list.append([sD,eD,0]) 

            i +=1
            mib = reach[i]

        ##If i is the start of a steep up area.
        elif reach[i] < reach[i+1]*(1-ksi):
            #print("STEEP UP starting in {}".format(i))
            steep_down_list = update_and_filter(steep_down_list,reach,i,ksi,mib)

            sU = i  # start of current upward region.
            color[i] = "g"

            # Go to the end of the steep area.
            sigma = 0
            while (reach[i+1]>=reach[i] 
                   and sigma < M 
                   and i+1<len(reach)-2
                   and not reach[i+1]*(1-ksi) > reach[i+2]):
                # Test if there is no more than M consecutives
                # points that are not ksi-steep.              
                if reach[i+1] <= reach[i+2]*(1-ksi): #i is ksi steep-upward.
                    sigma = 0
                else:
                    sigma += 1
                i += 1
                color[i] = "g"
            #print("Ending in {}\n".format(i))
            

            eU = i #end of the current upward region.

            i +=1
            mib = reach[i]

            for D in steep_down_list:
                ## Compare the end pf the steep-up U area multiplied by (1-ksi)
                ## with the mib value of the steep down area D thus satisfying 
                ## the condition (sc2*).
                #print("Compare {} and local mib {}".format(reach[eU]*(1-ksi),D[2]))
                if reach[eU]*(1-ksi) !=  D[2]:
                    s,e = compute_cluster((sU,eU),D,reach,ksi)
                    #print("Cluster possible in [{},{}] ({})\n".format(s,e,nc))
                    nc += 1
                    if cluster_cond(s,e,reach,M,ksi):
                        clusters.append((s,e))
                        #print("Cluster confirmed\n")
        else:
            i += 1
    return clusters, color

def update_and_filter(steep_down,r,i,ksi,global_mib):
    sd = []
    for (s,e,mib) in steep_down:
        mib = max(mib,global_mib)
        if r[s] * (1-ksi) >= global_mib:
            sd.append([s,e,mib])
    return sd

def compute_cluster(U,D,r,ksi):
    if r[D[0]] * (1-ksi) >= r[U[1]+1]:
        s = max([x for x in range(D[0],D[1]+1) if (r[x] > r[U[1]+1])])
        return(s,U[1])
    elif r[U[1]+1] * (1-ksi) >= r[D[0]]:
        e = min([x  for x in range(U[0],U[1]+1) if (r[x] < r[D[0]])])
        return(D[0],e)
    else:
        return(D[0],U[1])

def cluster_cond(s,e,r,M,ksi):
    if e-s < M:
        return False
    elif not r[s] * (1-ksi) >= r[s+1]:
        return False
    elif not r[e] <= r[e+1] * (1-ksi):
        return False
    else:
        return True


import visual
def test(pts):
    try:
        o,r = optics(np.transpose(pts),100000,15)
    except:
        o = 0
        r = 0
    
    try:
        c,co = find_cluster(r,0.001,15)
    except:
        c = 0
        co = 0 
    try:
        visual.plot_clust(r,co,c,False)
    except:
        print("Error in plot_clust")
    try:
        visual.plot_pclust(c,pts,r,o,True)
    except:
        print("Error in pclust")
    return o,r,c


def find_cluster_threshold(reach,threshold=None,M=5):
    color = ["k"] * len(reach)
    if threshold == None:
        threshold = max(reach)/3. #max(reach)-min(reach))/5. + min(reach)
    clusters = []
    s = None
    e = None
    for k,(r,rm1) in enumerate(zip(reach[1:],reach[:-1])):
        if s and not e:
            color[k] = "b"
        if r > threshold and rm1 < threshold and s != None and k-1-s > M:
            e = k-1
            clusters.append((s,e))
            s = None
            e = None
        elif r < threshold and rm1 > threshold:
            s = k-1
    return clusters, color

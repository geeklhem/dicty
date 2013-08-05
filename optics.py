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

    Output a dict containing :
    - order : List of n int. order[k] gives the position
    in the original array of the k-th point visited during the 
    walk.
    - reach: List of n float. Smallest reachability distance   
    with respect to a point on the left of k in the walk's order.
    - ofs: List of n float. Outliner factor. The OF is the  
    average of the ratios of the local reachability distance
    (lrds) of the M nearest neighbors and the one of k. 
    """

    nb_points = len(points)
    o = 0
    order = [0] * nb_points
    reach = [None] *  nb_points
    processed = [0] * nb_points
    lrd = [None] * nb_points
       
    for p in range(nb_points):
       
        # for each unprocessed point p of DB
        if not processed[p]:
            # Get p neighbors (eps-radius)
            N = geo.get_neighbors(p,eps,points)
        
            # mark p as processed and output it to the ordered list
            processed[p] = 1
            order[o]=p
            lrd[o] = compute_lrd(p,N,M,points)
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
                        lrd[o] = compute_lrd(q,N2,M,points)
                        o += 1  

                        #if (core-distance(q, eps, Minpts) != UNDEFINED)
                        if core_dist(q,N2,M,points) != None:
                            update(N2,q,seeds,eps,M,points,processed,reach)


    ofs = compute_ofs(lrd,order)
    reach = [reach[o] if reach[o] != None else 0 for o in order]
    return {"order":order,"reach":reach,"ofs":ofs}

def core_dist(p,N,M,points):
    """Return the list of the distance to the M nearest points"""
    if len(N) < M:
        return None
    else:
        dist = [geo.d(points[p],points[n]) for n in N]
        dist.sort()
        return dist[M-1]


def update(N,p,seeds,eps,M,points,processed,reach):
    """Update the priority queue """
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
                    #update, if the reachability distance is lowered, then
                    #this point should change position in the queue.
                    seeds.add(n,new_rd) 


#~~~~~~~~~~
# OPTICS-OF
#~~~~~~~~~~

def compute_lrd(p,N,M,points):
    """Compute the local reachability distance of p
    Which is the inverse of the average reachability distance
    between p and its M nearest neighbors"""
    rd = [0] * M
    cd = core_dist(p,N,M,points)
    neighM = N[:M]

    for k,n in enumerate(neighM):
        rd[k] = max(cd,
                    geo.d(points[p],points[n]))
    lrd =  1/(sum(rd)/float(len(rd)))
    return [lrd,neighM]

def compute_ofs(lrd,order):
    """ Compute the outlier factor of all points.
    The OF is the average of the ratios of the local
    reachability distance (lrds) of the M nearest neighbors
    and the one of p.
    It captures the degree of "outlierness" of p."""
    ofs = [None] * len(order)
    
    for j,l in enumerate(lrd):
        
        lrd_neighbors = [lrdn[0]/l[0] 
                         for k,lrdn
                         in enumerate(lrd) 
                         if order[k] in l[1]]
        ofs[j] = (sum(lrd_neighbors)
                  /float(len(lrd_neighbors)))
    return ofs


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Priority queue
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class PriorityQueue(object):
    """A priority queue inspired by 
    http://docs.python.org/2/library/heapq.html"""
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Cluster extraction
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def find_clusters_sander(reach,ratio=0.80,M=5):
    """ Extract cluster from the reachability distance plot 
    using the method described in Sander et al. 2003 : 
    Automatic Extraction of Clusters from Hierarchical Clustering Representations"""
    
    color = ["k"] * len(reach)

    ## STEP ONE : Get an ordered list of locals maxima ##
    P = []
    for i,r in enumerate(reach):
        if i>0:
            if (r >= max(reach[max(0,i-M):i]) 
                and r >= max(reach[i:min(i+M,len(reach))])):
                P.append((i,r))
                color[i] = "r"
    P.sort(key=lambda x:-x[1]) 


    ## STEP TWO : Process the list and build the cluster tree ##
    T = {"childs":[],
         "s":0,"e":len(reach),
         "N":len(reach),
         "parent":None,"i":1,"depth":0,
         "leaf":False}
    clusters = [T]
    
    cluster_tree(T,P,ratio,reach,M,clusters)

    return clusters, color

def cluster_tree(N,P,ratio,reach,M,clusters):
    """
    Recursively construct the cluster tree.  
    :param N: current node; the root of the tree in the first call.
    :type N: dict
    :param P: Local maxima points, sorted in descending order of reachability
    :type P: list
    """
    
    #    print("Current node is now {} and span from {} to {}".format(N["i"],
    #                                                                   N["s"],
    #                                                                   N["e"]))
            

    # Terminal recursion
    if len(P) == 0:
        N["leaf"] = True
        return None 
    
    # take the next largest local maximum point as a possible separation between clusters
    N["split_point"] = P[0]
    s = P[0] 
    P.pop(0)

    # STEP ONE : Two new nodes are created.
    
    # N1 = p in N["points"] | p is left of s in the reachability plot 
    N1 = {"s":N["s"],
          "e":s[0],
          "N":s[0]-N["s"],
          "childs":[],
          "leaf":False}

    # N2 = p in N["points"] | p is right of s in the reachability plot  
    N2 = {"s":s[0]+1,
          "e":N["e"],
          "N":N["e"]-s[0]-1,
          "childs":[],      
          "leaf":False}
    L1 =  [p for p in P if p[0]<s[0]] #p in L | p lies to the left of s in the reachability plot
    L2 =  [p for p in P if p[0]>s[0]] #p in L | p lies to the right of s in the reachability plot
    NL = [(N1,L1),(N2,L2)]

    points_N1 = reach[N1["s"]:N1["e"]]
    points_N2 = reach[N2["s"]:N2["e"]]

    if N1["N"] > M:
        N1["mean_r"] = sum(points_N1)/float(len(points_N1))
        N1["sd"] = math.sqrt((float(1)/(N1["N"]-1)) 
                             * sum([(p-N1["mean_r"])**2 for p in points_N1]))
    else :
        N1["mean_r"] = 0

    if N2["N"] > M:
        N2["mean_r"] = sum(points_N2)/float(len(points_N2))
        N2["sd"] = math.sqrt((float(1)/(N2["N"]-1)) 
                             * sum([(p-N2["mean_r"])**2 for p in points_N2]))
    else: 
        N2["mean_r"] = 0

    # STEP TWO : Is the new separation significant 
    if (N1["mean_r"]/s[1] > ratio) and (N2["mean_r"]/s[1] > ratio):
    #if split point s is not significant, ignore s and continue
        return cluster_tree(N,P,ratio,reach,M,clusters)

    # STEP THREE : Does the nodes are large enough ?
    else:
        if N1["N"] < M:
            NL.pop(0)
        if N2["N"] < M:
            NL.pop(-1)
        if len(NL) == 0:
            N["leaf"] = True
            return None #Parent_of_N is a leaf.

        # STEP FOUR : Add the nodes to the tree
        # check if the nodes can be moved up one level = if the current node is within one
        # standard deviation of the parent's value. 
        if (N["parent"] != None 
            and N["parent"]["parent"] != None
            and abs((s[1]-N["parent"]["mean_r"])/N["parent"]["sd"]) < 1) :
            # Adds the NLs to the parent node.
            for Ni,Li in NL:
                Ni["parent"] = N["parent"]
                Ni["i"] = len(clusters)+1
                Ni["depth"] = N["depth"]

                N["parent"]["childs"].append(Ni)
                clusters.append(Ni)


        else:
            # Adds the NLs to the current node
            for Ni,Li in NL:
                Ni["parent"] = N
                Ni["i"] = len(clusters)+1            
                Ni["depth"] = N["depth"] +1

                N["childs"].append(Ni)
                clusters.append(Ni)

        # RECURSIVITY !
        for (Ni,Li) in NL:
            cluster_tree(Ni,Li,ratio,reach,M,clusters)



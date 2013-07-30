#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""Partition frames using the group informations and the cells positions.
Always take a list of groups and a list of points and 
return a len(points) list with the index of the group they belong to """

import math
import itertools
import numpy as np
import optics as op
import copy
import geometry as geo

def voronoi(center,points):
    def d(a,b):
        return math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)

    attribution = [0]*len(points[0,:])

    for k,x,y in itertools.izip(itertools.count(),points[0,:],points[1,:]):
        nearest = 0
        mdist = d((center[0,0],center[1,0]),(x,y))+1

        for i,cx,cy in itertools.izip(itertools.count(),
                                      center[0,:],
                                      center[1,:]):
            dist = d((cx,cy),(x,y))
            if dist < mdist:
                mdist = dist
                attribution[k] = i
    return attribution


def optics_th(points,threshold,eps=9000,M=15):
    pts = np.transpose(points)
    order,reach = op.optics(pts,eps,M)
    c = 1
    attr = [0]*(len(order)+1)
    already = 0
    for r,o in zip(reach,order):
        if r > threshold and not already:
            c += 1 
            already = 1
            attr[o] = 0
        elif r>threshold and already: 
            already = 0
            attr[o] = 0
        else:
            attr[o] = c
    return attr,reach,order

def optics_clust(points,eps=9000,M=15,ksi=0.001,method="threshold"):
    pts = np.transpose(points)
    optics_results = op.optics(pts,eps,M)
    order = optics_results["order"]
    reach = optics_results["reach"]
    
    
    if method == "ksistep":
        clust_tuples,color_histo = op.find_cluster(reach,ksi,M)    
        clusters = [cluster_dict(c,order, points)
                    for c 
                    in clust_tuples]



    elif method == "threshold":
        clust_tuples,color_histo = op.find_cluster_threshold(reach,M=M)
        clusters = [cluster_dict(c,order, points)
                    for c 
                    in clust_tuples]


    elif method == "sander":
        clusters,color_histo = op.find_clusters_sander(reach,
                                                    ratio=0.75,
                                                    M=M)
        clusters = [exclude_loners(c,
                                   optics_results["ofs"],
                                   ofmax=1.1,
                                   M=M)
                    for c in clusters]

        clusters = [c for c in clusters if c != None]

        clusters = [cluster_properties(c,order,points)
                    for c in clusters]

        

    else:
        raise NameError


    attribution = [None]*len(pts)
    for k,c in enumerate([cl for cl in clusters if cl["leaf"]]):
        for o in order[c["s"]:c["e"]]:
            attribution[o] = k
    # attribution = [a 
    #                if optics_results["ofs"][order.index(k)] < 1.
    #                else None
    #                for k,a in enumerate(attribution)]
    loners = sum([1 for a in attribution if a==None])

    return {"attribution":attribution,
            "reach": reach,
            "order":order,
            "clusters":clusters,
            "loners":loners,
            "color_histo":color_histo,
            "ofs":optics_results["ofs"]}

def exclude_loners(cluster,ofs,ofmax,M):
    """Exlude particles with a too high outliner factor."""
    while ofs[cluster["s"]] > ofmax:
        cluster["s"] += 1
    while ofs[cluster["e"]-1] > ofmax:
        cluster["e"] -= 1
    cluster["N"] = cluster["e"] - cluster["s"]
    if cluster["N"] <= M:
        return None
    else:
        return cluster

        
def cluster_properties(cluster,order,points):
    if "leaf" in cluster.keys() and cluster["leaf"]:
        cluster["points_index"] = order[cluster["s"]:cluster["e"]]
        x = points[0,cluster["points_index"]]
        y = points[1,cluster["points_index"]]
        cluster["x_mean"] = sum(x)/cluster["N"]
        cluster["y_mean"] = sum(y)/cluster["N"]
        cluster["centroid"] = (cluster["x_mean"],cluster["y_mean"])
        cluster["convex_hull"] = geo.convex_hull(np.array((x,y)))
    return cluster

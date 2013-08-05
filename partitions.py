#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""Partition frames using the group informations and the cells positions. """

import math
import itertools
import numpy as np
import optics as op
import copy
import geometry as geo

def voronoi(center,points):
    """ Partitions the points using voronoi cells.
    :param center: Voronoi cells center's position.
    :type center: list
    :param points: Points position.
    :points type: np.array
    :return: a dict. with a len(points) list with the index of the group they belong to"""

    d = geo.d

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
    return {"attribution": attribution}


def optics_clust(points,eps=9000,M=15):
    """ Make the partition using OPTICS clustering algorithm. 
    :return: A dictionary.
    
    - "attribution": a len(points) list with the index of the group
    they belong to.
    - "reach": Minimal reachability distance of a point to the previous
    in the walk.
    - "order": Order of the points in the walk,
    - "clusters": List of dicts. containing clusters information.,
    - "loners": Number of loners,
    - "color_histo": Color of each bar of the reachability plot,
    - "ofs": Local outlier factor of each point.
    """
    pts = np.transpose(points)
    optics_results = op.optics(pts,eps,M)
    order = optics_results["order"]
    reach = optics_results["reach"]
    
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

    attribution = [None]*len(pts)
    for k,c in enumerate([cl for cl in clusters if cl["leaf"]]):
        for o in order[c["s"]:c["e"]]:
            attribution[o] = k
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
    """Complete a leaf cluster's dictionnary 
    by computing additional properties."""
    if cluster["leaf"]:
        cluster["points_index"] = order[cluster["s"]:cluster["e"]]
        x = points[0,cluster["points_index"]]
        y = points[1,cluster["points_index"]]
        cluster["x_mean"] = sum(x)/cluster["N"]
        cluster["y_mean"] = sum(y)/cluster["N"]
        cluster["centroid"] = (cluster["x_mean"],cluster["y_mean"])
        cluster["convex_hull"] = geo.convex_hull(np.array((x,y)))
    return cluster

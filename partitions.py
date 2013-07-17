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
    order,reach = op.optics(pts,eps,M)
    
    if methods == "ksistep":
        clusters,color_histo = [cluster_dict(c,order, points[f])
                                for c 
                                in op.find_cluster(reach,ksi,M)]

    elif methods == "threshold":
        clusters,color_histo = [cluster_dict(c,order, points[f])
                                for c 
                                in op.find_cluster_threshold(reach,M=M)]
    else:
        raise NameError

    attribution = [None]*len(points)
    for k,c in enumerate(clusters):
        for o in order[c["s"]:c["e"]]:
            attribution[o] = k
    loners = sum([1 for a in attribution if a!=None])

    return {"attribution":attribution,
            "reach": reach,
            "order":order,
            "clusters":clusters,
            "loners":loners,
            "color_histo":color_histo}

def cluster_dict(clust,order,points):
    f = {}
    f["s"],f["e"] = clust
    f["N"] = f["e"] - f["s"]
    f["points_index"] = order[f["s"]:f["e"]]

    x = points[0,f["points_index"]]
    y = points[1,f["points_index"]]
    f["x_mean"] = sum(x)/f["N"]
    f["y_mean"] = sum(y)/f["N"]
    f["centroid"] = (f["x_mean"],f["y_mean"])
    f["convex_hull"] = geo.convex_hull(np.array((x,y)))
    return f

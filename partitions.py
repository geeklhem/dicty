#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""Partition frames using the group informations and the cells positions.
Always take a list of groups and a list of points and 
return a len(points) list with the index of the group they belong to """

import math
import itertools
import numpy as np
import optics as op

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


def optics(points,threshold):
    pts = np.transpose(points)
    order,reach = op.optics(pts,500,5)
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
    return attr,reach

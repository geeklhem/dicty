#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
""" Geometrical functions"""

import numpy as np
import math



def get_neighbors(p,eps,points,mx=0):
    ## O(n²) !!!!
    if eps == None:
        eps = float("inf")

    neighbors = []
    for k,q in enumerate(points):
        dist =  d(points[p],q)
        if dist < eps:
            neighbors.append(k)
            if mx:
                distances.append(dist)
    if mx:
        zipped = zip(neighbors,dist)
        zipped.sort(key=lambda x:x[1])
        return zipped[0:mx]
    else:
        return neighbors

def convex_hull(points):
    """ Return a list of points forming the convex hull of
    a set of points using Graham Scan algorithm."""

    # Find the point with the lowest y-coordinate. 
    ymin = min(points[1,:])
    iymin = [i for i,j in enumerate(points[1,:]) if j==ymin]
    
    # If the lowest y-coordinate exist in more than one points in the set...
    if len(iymin) > 1:
        # ...take the point with the lowest x-coordinate out of the candidates.
        ip = points[0,:].index(min(points[0,iymin]))
    else:
        ip = iymin[0]

    p = (points[0,ip],points[1,ip])
    # Order the set of points in increasing order of the angle they and the point P
    # Make with the x-axis.
    angles = [(math.acos((x-p[0])/d((x-p[0],y-p[1]),(0,0))),(x,y)) 
              for x,y
              in zip(points[0,:],points[1,:])]
    angles.sort(key=lambda x:x[0])
    print(p)
    hull = [angles[0][1],angles[1][1]]
    for a,xy in angles[2:]:
        while len(hull)>=2 and cross_product(hull[-2],hull[-1],xy)<=0:
            hull.pop()
        hull.append(xy)
    return np.array(hull),angles


def cross_product(p1,p2,p3):
    """ Three points p1,p2 and p3 are a :
 
     - counter clockwise turn if cp > 0 
     - clockwise turn if cp < 0
     - colinear if cp = 0"""
    cp1 = (p2[0]-p1[0])*(p3[1]-p1[1]) 
    cp2 = (p2[1]-p1[1])*(p3[0]-p1[0]) 
    return cp1-cp2


def d(a,b):
    """Return the euclidian distance between a and b"""
    return math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)


def intersect_polygons(clip,subject):
    """Sutherland-Hodgman algorithm"""
    # Put the vertices in the right order
    def right_order(polygon):
        if cross_product(polygon[0],polygon[1],polygon[2]) > 0:
            polygon.reverse()
            return polygon
        else :
            return polygon
    clip = right_order(clip)
    subject = right_order(subject)

    # Use float for non int intersections points.
    clip = [(float(x),float(y)) for x,y in clip]
    subject = [(float(x),float(y)) for x,y in subject]

    output = subject
    clipEdges = [(v1,v2) for v1,v2 in zip(clip[:-1],clip[1:])] + [(clip[-1],clip[0])]
    for edge in clipEdges:
        inpt = [(v1,v2) for v1,v2 in zip(output[:-1],output[1:])] + [(output[-1],output[0])] 
        output = []
        for s,e in inpt:
            if cross_product(edge[0],edge[1],e)<0:
                if not cross_product(edge[0],edge[1],s)<0:
                    output.append(intersection((s,e),edge))
                output.append(e)
            elif cross_product(edge[0],edge[1],s)<0:
                output.append(intersection((s,e),edge))
    return output

def intersection(a,b):
    """return the intersection of two lines (a[0],a[1] and (b[0],b[1])"""
    X = (a[0][0],a[1][0],b[0][0],b[1][0])
    Y = (a[0][1],a[1][1],b[0][1],b[1][1])
    
    d = (X[0]-X[1]) * (Y[2]-Y[3]) - (Y[0]-Y[1]) * (X[2]-X[3])

    x = (X[0]*Y[1] - Y[0]*X[1])*(X[2]-X[3]) - (X[0]-X[1])*(X[2]*Y[3]-Y[2]*X[3])
    x /= d 

    y = (X[0]*Y[1] - Y[0]*X[1])*(Y[2]-Y[3]) - (Y[0]-Y[1])*(X[2]*Y[3]-Y[2]*X[3])
    y /= d 

    return x,y

def area(polygon):
    """Return the area of a given polygon
    :param polygon: Polygon, list of vertices
    :type polygon: list
    :return: (float) - Area of the polygon.
    """
    area = 0
    # j is the previous vertex to i
    for i,j in zip(polygon,polygon[-1:]+polygon[0:-1]):
        area += (j[0]+i[0]) * (j[1]-i[1])
    return abs(float(area/2))
pts = np.loadtxt("pts")


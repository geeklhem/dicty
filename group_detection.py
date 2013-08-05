#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""Detecting groups of cells using various technics. Always retruns """

import numpy as np
import math

def from_file(centerfile, csl=1000):
    """ Import groups from a csv file written with imageJ.
    Format : Id, Area, X, Y
    :param centerfile: csv file.
    :type centerfile: string
    :param csl: Minimal particle size to be considererd as a group.
    :type csl: int
    :return: Dict (np.array, np.array, int) - pos:Group position,
    area:Group area, N:Number of groups.
    """
    #Load data
    final = np.genfromtxt(centerfile,
                          usecols=(1,2,3),
                          skip_header=1,
                          delimiter=",",
                          names=("area","x","y"))
    #Get groups position
    centers = np.transpose(np.array([(x,y) for ar,x,y in final if ar>csl]))
    center_size = np.array(([ar for ar,x,y in final if ar>csl]))
    C = len(centers[0,:])

    return {"pos":centers, "area":center_size, "N":C}

def particle_size(data,frame,csl):
    """ Import groups from a frame of the data object using particle size

    .. warning ::
       Not implemented yet

    :param centerfile: data object.
    :type centerfile: data.Data
    :param frame: Frame number
    :type frame: int
    :param csl: Minimal particle size to be considererd as a group.
    :type csl: int
    :return: Dict (np.array, np.array, int) - pos:Group position,
    area:Group area, N:Number of groups.
    """
    return {"pos":[], "area":[], "N":0}

def from_attribution(points,attr,loners_index=0):
    """ Detect groups from a clustering analysis

    :return: Dict (np.array, np.array, int) - pos:Group position,
    area:Group area, N:Number of groups.
    """
    N = max(attr)+1
    # max_x = [0]*N
    # min_x = [float("inf")]*N
    sum_x = [0]*N
    # max_y = [0]*N
    # min_y = [float("inf")]*N
    sum_y = [0]*N
    number = [0]*N
    loners= 0

    for i,p in enumerate(zip(points[0,:],points[1,:])):
        g = attr[i] #get the group it belongs to.
        if g != loners_index:
            # max_x[g] = max(max_x[g], p[0])
            # min_x[g] = min(min_x[g], p[0])
            sum_x[g] += p[0]
            # min_y[g] = min(min_y[g], p[1])
            # max_y[g] = max(max_y[g], p[1])
            sum_y[g] += p[1]
            number[g] += 1
            # print("group {} +1 ({})".format(g,number[g]))
        else:
            loners += 1
            # print("One more loner ({})".format(loners))
    pos = np.transpose(np.array([(sx/n,sy/n) if n!=0
                                 else (0,0)
                                 for sx,sy,n in zip(sum_x,sum_y,number)]))
    area = [g**2*math.pi for g in number] # The radius is the number of elements

    return {"pos":pos, "area":area, "N":N,"cbs":number,"loners":loners}


if __name__ == "__main__":
    pass

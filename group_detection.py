#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""Detecting groups of cells using various technics. Always retruns """

import numpy as np


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

def density(data,frame,threshold):
    """ Import groups from a frame of the data object using particle density
    :param centerfile: data object.
    :type centerfile: data.Data
    :param frame: Frame number
    :type frame: int
    :param threshold: Density threhold between groups
    :type threshold: float
    :return: Dict (np.array, np.array, int) - pos:Group position,
    area:Group area, N:Number of groups.
    """
    return {"pos":[], "area":[], "N":0}

def algo(data,frame):
    """ Import groups from a frame of the data object using an algorithmic method.
    :param centerfile: data object.
    :type centerfile: data.Data
    :param frame: Frame number
    :type frame: int
    :param threshold: Density threhold between groups
    :type threshold: float
    :return: Dict (np.array, np.array, int) - pos:Group position,
    area:Group area, N:Number of groups.
    """
    return {"pos":[], "area":[], "N":0}

if __name__ == "__main__":
    pass

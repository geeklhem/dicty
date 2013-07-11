#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""Visualisation functions"""

import itertools
import matplotlib.pyplot as plt
import collections
from matplotlib.patches import Circle
import math
import numpy as np


def get_color():
       return itertools.cycle(['r', 'g', 'b', 'c', 'm', 'y', 'k'])

def plot_particle(points,attribution,xlim,ylim,groups=None,cbs=None,show=True):
    #Colors
    color_iter = get_color()
    if groups:
        colors = [next(color_iter) for c in range(groups["N"])]
    else:
        colors = [next(color_iter) for c in range(max(attribution)+1)]
    colorlist = [colors[i] for i in attribution]

    # Points
    plt.scatter(points[0,:],points[1,:],color=colorlist)
       
    if groups:
       # Groups
       ax = plt.gca()
       for x,y,c,s in zip(groups["pos"][0,:],groups["pos"][1,:],colors,groups["area"]):
              ax.add_artist(Circle(xy=(x,y),
                                 radius=math.sqrt(s/math.pi),
                                 facecolor=c,
                                 alpha=0.3))
       if cbs:
              for x,y,s in zip(groups["pos"][0,:],groups["pos"][1,:],cbs):
                     plt.text(x, y,"[{}]".format(s),alpha=0.5)
    # Axis details
    plt.xlim((0,xlim))
    plt.ylim((0,ylim))

    if show:
        plt.show()

def correlation(x,y,regline,xlab,ylab,show=True):
    plt.plot(x,regline)
    plt.scatter(x,y)
    ax = plt.gca()
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if show:
        plt.show()

def time(x,name="value",xlim=None,show=True):
    plt.plot(x)
    if xlim:
           plt.ylim((0,xlim))
    ax = plt.gca()
    ax.set_ylabel(name)
    ax.set_xlabel("frame")
    if show:
       plt.show()

def distrib(x,show=True):
    plt.plot(x)
    ax = plt.gca()
    ax.set_ylabel("Distribution")
    ax.set_xlabel("Group size")
    if show:
       plt.show()
       

def areaplot(cbc,show=True):
       #Colors
       color_iter = get_color()
       
       sections = zip(*cbc)
       x = range(len(cbc))
       y = [0]*len(cbc)  
       colorlist = [next(color_iter) for c in range(len(cbc[0]))]

       for ynow,c in zip(sections,colorlist):
              yc = [i+j for i,j in zip(ynow,y)]
              plt.plot(x,yc,color=c)
              plt.fill_between(x,y,yc,color=c,alpha=0.5)
              y = [i for i in yc]
     
       if show:
              plt.show()

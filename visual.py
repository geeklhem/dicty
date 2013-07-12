#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""Visualisation functions"""

import itertools
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import math
import random
import numpy as np


def get_color():
       index = [i/10. for i in range(11)]
       colors = ([(1,0,i) for i in index] +
                 [(1-i,0,1) for i in index]+
                 [(0,i,1) for i in index]+
                 [(0,1,1-i) for i in index]+
                 [(0,1-i,0) for i in index]+
                 [(i,0,0) for i in index]) 
       random.shuffle(colors)
       colors = ["k"] + colors #non optimal...
       return itertools.cycle(colors)
#       return itertools.cycle(['r', 'g', 'b', 'c', 'm', 'y', 'k'])

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

       
def optics_reachability(reach, show=True):
       """Plots optics algorithm results"""     
       plt.bar(range(len(reach)), reach, facecolor='k')
       ax = plt.gca()
       ax.set_ylabel("Reachability distance")
       ax.set_xlabel("Points")
       if show:
              plt.show()


def plot_clust(r,co,cl,show=True):
    plt.bar(range(len(r)),r,color=co,edgecolor=co); 
    clust = [(e-s,s,e) for s,e in cl]
    clust = sorted(clust, key=lambda x:x[0])
    ystep = 700./len(clust)
    y = 0
    crossed = [0] * len(r)
    line_breack=False
    for l,s,e in clust:
       for x in range(max(0,s-3),min(e+3,len(r))):
              if not crossed[x]:
                     crossed[x] = 1
              else:
                     line_breack = True
                     crossed = [0] * len(r)
       if line_breack:
              y = y-ystep
              line_breack=False
       plt.hlines(xmin=s,xmax=e,y=y, linewidth=1)
    if show:
           plt.show()

from matplotlib.collections import PolyCollection
import matplotlib as mpl

def plot_pclust(clust,pts,r,order,show=True):
    verts = []
    for c in clust:
        pts_order = order[c[0]:(c[1])]
        pts_list = [pts[:,o] for o in pts_order]
        verts.append(pts_list) 


    z = np.random.random(len(verts)) * 500

    fig, ax = plt.subplots()

    # Make the collection and add it to the plot.
    coll = PolyCollection(verts, array=z, cmap=mpl.cm.jet,alpha=0.5,facecolor="None")
    ax.add_collection(coll)
    
    ax.autoscale_view()
    plot_particle(pts,[0]*len(pts),6000,5000,show=show)

from matplotlib.path import Path
import matplotlib.patches as patches

def plot_path(pts,o,show=True):
    pts_list = [pts[:,o] for o in order]
    codes = [Path.MOVETO] + [Path.LINETO]*(len(pts_list)-1)
    path = Path(pts_list, codes)
    fig, ax = plt.subplots()
    patch = patches.PathPatch(path, facecolor="None", lw=2)
    ax.add_patch(patch)
    ax.autoscale_view()
    if show:
        plt.show()



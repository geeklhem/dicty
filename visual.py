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

def plot_particle(points,
                  attribution,
                  xlim,ylim,
                  clusters,
                  show=True):

    colorlist = [clusters[a]["color"] if a != None else "k" for a in attribution]

    # Points
    plt.scatter(points[0,:],points[1,:],color=colorlist)
    
    for cluster in clusters:
              plot_polygon(cluster["convex_hull"],cluster["color"],show=False)
              plt.text(cluster["x_mean"], cluster["y_mean"],
                       "[{}]".format(cluster["N"]),
                       alpha=0.5,
                       horizontalalignment='center', verticalalignment='center')
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
    clust = [(c["N"],c["s"],c["e"]) for c in cl]
    clust = sorted(clust, key=lambda x:x[0])
    ystep = 700./len(clust)
    y = -10
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


def plot_polygon(polygon,color,alpha=0.3,show=True):
    codes = [Path.MOVETO] + [Path.LINETO]*(len(polygon)-1)
    path = Path(polygon, codes)
    ax = plt.gca()
    patch = patches.PathPatch(path, facecolor=color,alpha=alpha, lw=2)
    ax.add_patch(patch)
    if show:
       ax.autoscale_view()
       plt.show()


def tree_trace(tclusters,clusters,show=True):
    previous_y = 0
    ax = plt.gca()
    lim_x,lim_y = 0,0
    #print clusters
    maxn = float(max([max([c["N"] for c in f]) for f in clusters]))
    minn = float(min([min([c["N"] for c in f]) for f in clusters]))
    for it,tclust in enumerate(tclusters):
        for f,frm in enumerate(tclust["cluster_indicies"]):
            current_y = f * 20
            current_x = 20*it
            previous_x = [current_x,]
            lim_x = max(current_x,lim_x)
            lim_y = max(current_y,lim_y)
            x_step = 20/(len(frm)+1.)
            max_r = 0
            for c in frm:
                if c:
                    n = clusters[f][c]["N"]
                    r = 0.3*x_step+0.6*x_step*(n-minn)/maxn
                    max_r = max(r,max_r)
                    ax.add_artist(Circle(xy=(current_x,current_y),
                                         radius=r,
                                         facecolor=clusters[f][c]["color"],
                                         alpha=0.5))
                    plt.text(current_x,
                             current_y,
                             "{}".format(n),
                             horizontalalignment='center',
                             verticalalignment='center')
                    #    if f > 0:
                #         for px in previous_x:
                #             plt.arrow(px,previous_y,
                #                       current_x-px,current_y-max_r-previous_y,
                #                       facecolor=clusters[f][c]["color"],
                #                       width=x_step/10.,
                #                       length_includes_head=True,
                #                       alpha = 0.5)
                # 
                current_x += x_step
                previous_y = current_y - max_r 
    plt.xlim((-10,lim_x+10))
    plt.ylim((-10,lim_y+10))
    if show:
           plt.show()

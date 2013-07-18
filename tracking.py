#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
""" Tracking functions"""

import numpy as np
import math

import geometry as geo


def track_cluster(clusters_frame):
    """ Return a list of list, giving the indices of the clusters
    tracking[frame][cluster] = cluster index or None"""
    tracking = []
    nb_clust = len(clusters_frame[0])
    # Each 
    tracking.append([(None,i) for i in range(len(clusters_frame[0]))])
    centroids = [c["centroid"] for c in clusters_frame[0]]

    #loop through frames starting at the second one
    for f,clusters in enumerate(clusters_frame[1:]):
        tracking.append([])
        for index,clust in enumerate(clusters):

            # Get the nearest neighbors ordered by distance
            neighbors = geo.get_neighbors(-1,
                                          None,
                                          centroids+[clust["centroid"]],
                                          mx = 3)
            neighbors = neighbors[1:] #supress the first element (same point, dist =0)
            clust["neighbors"] = []

            #print("clust {}".format(clust))
            for ni,dist in neighbors:
                #print("index {} (dist {}) in {}".format(ni,dist,len(clusters_frame[f])))
                #print("-----\n{}".format(clust))
                clip = geo.intersect_polygons(
                    clusters_frame[f][ni]["convex_hull"],
                    clust["convex_hull"] )
                area = geo.area(clip)
                clust["neighbors"].append({"index":ni,
                                           "distance":dist,
                                           "inter_area":area})
            
            #Order neihgbors by increasing overlap area
            clust["neighbors"].sort(key=lambda x:-x["inter_area"])

            # Take the most overlaping one (if there's an overlap).
            if len(clust["neighbors"]) and clust["neighbors"][0]["inter_area"] != 0:
                #print "Tracking {} : {}".format(f,tracking[f])
                prev_index_list = [x[1] for x in tracking[f]]
                #print "Neigh : {}".format(clust["neighbors"])
                #print "Index of {} in {} ".format(clust["neighbors"][0]["index"], prev_index_list)
                
                prev_index = prev_index_list.index(clust["neighbors"][0]["index"])
                tracking[f+1].append((prev_index,index))
            # Or just start a new trace
            else:
                tracking[f+1].append((None,index))
                
        # get the centroid of this frame for the treatment of the following frame
        centroids = [c["centroid"] for c in clusters]
    return tracking

def traces(tracking):
    nb_start = 0
    traces = []
    delete = []
    for f,trk in enumerate(tracking):
        for t in traces:
            t.append([])
        for old,new in trk:
            if old == None:
                nb_start += 1 
                traces.append([])
                for i in range(f):
                    traces[-1].append([None])
                traces[-1].append([new])
            else:
                for j,l in enumerate(last):
                    if old in l:
                        p_index = j
                traces[p_index][f] += [new]

        for j,fused in enumerate(traces):
            if fused[-1] == []:
                for r,fusion in enumerate(traces):
                    if fusion[f-1] == fused[f-1] and fusion[f-1] != [] and j != r:
                        # traces fusion
                        #print("Fusion : {} with {}".format(fusion,fused))
                        for k,v in enumerate(fused[:f-1]):
                            fusion[k].extend(v)
                            #print("Add {} to {} (trace {})".format(v,fusion[k],k))
                        fused[f-1] = [None]
                        delete.append(j)
        #print "Traces: {}".format(traces)
        last = [t[-1] for t in traces] 
        #print "Last: {}".format(last)
    traces = [t for i,t in enumerate(traces) if i not in delete]
    #print(delete)
    return traces


fake_tracking = [[(None,1),(None,2),(None,3),(None,4),(None,5),(None,6),(None,7)],
                 [(1,4),(2,2),(2,1),(4,3),(5,3),(None,5),(6,6),(7,6)],
                 [(4,1),(2,2),(1,4),(3,3),(6,5),(6,6)],
                 [(1,3),(4,2),(3,1)]]

fake_clusters = [[{"N":90},{"N":10},{"N":10},{"N":10},{"N":10},{"N":10},{"N":10},{"N":10}], 
                 [{"N":10},{"N":10},{"N":100},{"N":10},{"N":10},{"N":10},{"N":10},{"N":10}],
                 [{"N":120},{"N":10},{"N":10},{"N":10},{"N":10},{"N":10},{"N":10},{"N":10}],
                 [{"N":10},{"N":10},{"N":110},{"N":10},{"N":10},{"N":10},{"N":10},{"N":10}]] 

# fake_traces = traces(fake_tracking) 
# fake_tclust=traced_clusters_list(fake_traces,fake_clusters)
# tree_trace(fake_tclust,fake_clusters)
# plt.show()

import colorsys
import copy

def traced_clusters_list(traces,clusters):
    tclusters = []
    current_color = [0,.4,1]
    hue_step = .9/len(traces)
    for it,trc in enumerate(traces):
        dct = {
            "cluster_by_frame":[],
            "N_agg":[],
            "primary_color":colorsys.hls_to_rgb(*current_color),
            "cluster_indicies":trc
        }
        Nagg = 0
        for f,frm in enumerate(trc):
            color = copy.copy(current_color)
            cbf = sum([1 for k in frm if k != None])
            for c in frm:
                if c != None:
                    Nagg += clusters[f][c]["N"] 
                    clusters[f][c]["color"] = colorsys.hls_to_rgb(*color)
                    if len(frm) > 1:
                        hue_small_step = 0.5*hue_step/(len(frm)-1)
                        lum_step = 0.4/(len(frm)-1)
                        color = [c+mod for c,mod in zip(color,[hue_small_step,lum_step,0])]

            dct["cluster_by_frame"].append(cbf)
            dct["N_agg"].append(Nagg)
        current_color[0] += hue_step
        tclusters.append(dct)
    return tclusters





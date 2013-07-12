#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""Image analysis"""

import export
import visual
import data 
import group_detection
import partitions
import analysis
import optics


def mock():
    # Load data and create export
    output = export.HtmlExport("mock_export")
    experiment = data.Data("data/stack.csv")
    experiment.groups = group_detection.from_file("data/centers35.csv")
    experiment.attribution = []
    experiment.cbs = []
    experiment.distrib = []
    experiment.corr_cbs_ls = []
    experiment.corr_cl_ls = []

    experiment.attr_loners = (partitions.voronoi(experiment.groups["pos"],
                                                 experiment.points[-1]))
    experiment.loners = analysis.cell_by_section(experiment.attr_loners,
                                                 experiment.groups["N"])

    mf = 32

    # ANALYSIS
    for f in range(mf):
        experiment.attribution.append(partitions.voronoi(experiment.groups["pos"],
                                                         experiment.points[f]))
        experiment.cbs.append(analysis.cell_by_section(experiment.attribution[f],
                                                       experiment.groups["N"]))
        experiment.distrib.append(analysis.group_size_distrib(experiment.cbs[f],
                                                              experiment.loners,
                                                              experiment.groups["N"]))
        experiment.corr_cbs_ls.append(analysis.lin_correlate(experiment.cbs[f],
                                                             experiment.groups["area"]))
        experiment.corr_cl_ls.append(analysis.lin_correlate(experiment.cbs[f],
                                                            experiment.loners))




    experiment.cc_cbs_ls = [c["corrcoef"] for c in experiment.corr_cbs_ls]
    experiment.cc_cl_ls = [c["corrcoef"] for c in experiment.corr_cl_ls]


    # EXPORT 
    output.add_text("Analysis of a film")


    output.add_title("Global")
    output.add_text("{} frames".format(experiment.frame_nb))

    output.add_fig("cbsls",
                   visual.time,(experiment.cc_cbs_ls,"Correlation",1))

    output.add_fig("clls",
                   visual.time,(experiment.cc_cl_ls,"Correlation",1))

    output.add_fig("cbs",
                   visual.areaplot,(experiment.cbs,),
                   proportions=(3,1))


    output.add_title("Frame by frame".format(f))
    for f in range(mf):
        output.add_title("Frame {0}".format(f),3)
        output.add_fig("particle_{0}".format(f),
                       visual.plot_particle,(experiment.points[f],
                                             experiment.attribution[f],
                                             experiment.X,
                                             experiment.Y,
                                             experiment.groups,
                                             experiment.cbs[f],
                                         ),
                       proportions=(2*float(experiment.X)/float(experiment.Y),2))

        output.add_fig("distrib_{0}".format(f),
                       visual.distrib,(experiment.distrib[f]["crowding"],))

        output.add_fig("cbs_corr_{0}".format(f),
                       visual.correlation,(experiment.cbs[f],
                                           experiment.groups["area"],                 
                                           experiment.corr_cbs_ls[f]["regline"],
                                           "Cell by section",
                                           "Group Area"))


        output.add_fig("cl_corr_{0}".format(f),
                       visual.correlation,(experiment.cbs[f],
                                           experiment.loners,                 
                                           experiment.corr_cl_ls[f]["regline"],
                                           "Cell by section",
                                           "Loners by sections"))

        output.add_text("""
        <strong>Mean crowding</strong> :   {0[mean_crowding]},
        <strong>Mean group size</strong> : {0[mean_gs]} 
        """.format(experiment.distrib[f]))

        output.add_text("""
        <strong>Linear Correlation : </strong> 
        Group Area =
        {0[a]} Cell by section + {0[b]} 
        <strong> Correlation coefficient :</strong> {0[corrcoef]} 
        <strong> Mean square error :</strong>  {0[mse]}) 
        """.format(experiment.corr_cbs_ls[f]))

        output.add_text("""
        <strong>Linear Correlation : </strong> 
        Loners by section =
        {0[a]} Cell by section + {0[b]} 
        <strong> Correlation coefficient :</strong> {0[corrcoef]}
        <strong> Mean square error</strong>  : {0[mse]}) 
        """.format(experiment.corr_cl_ls[f]))



    ##  export 
    output.export()

import random
def opt_exp(mif=0,maf=None,M=15):
    # Load data and create export
    output = export.HtmlExport("02_optics")
    print("Loading...")
    experiment = data.Data("data/stack.csv")
    experiment.attribution = []
    experiment.reach = []
    experiment.groups = []
    experiment.order = []
    experiment.clust = []
    experiment.distrib = []

    if not maf:
        maf = experiment.frame_nb
                        
    print("Analysis...")
    # ANALYSIS
    for f in range(mif,maf):
        print("{}/{}".format(f+1,maf))
        #random.shuffle(experiment.points[f])
        attr,reach,order,cluster,color = partitions.optics_clust(experiment.points[f],
                                                                 eps=9000,
                                                                 ksi=0.001,
                                                                 M=M)
        experiment.reach.append(reach)
        experiment.order.append(order)
        experiment.attribution.append(attr)
        experiment.groups.append(group_detection.from_attribution(experiment.points[f],
                                                                  attr,
                                                                  loners_index=0))
        experiment.clust.append((cluster,color))
        experiment.distrib.append(analysis.group_size_distrib(experiment.groups[-1]["cbs"],
                                                              experiment.groups[-1]["loners"],
                                                              experiment.groups[-1]["N"]))

    # EXPORT 
    print("Export")
    output.add_text("Analysis of a film")
    output.add_text("<strong>Clustering algorithm</strong> : Optics")
    output.add_text("<strong>Parameters</strong> : M = {}".format(M))


    output.add_title("Global")
    output.add_text("{} frames".format(experiment.frame_nb))

    output.add_title("Frame by frame".format(f))
    for f in range(maf-mif):
        print("{}/{}".format(f+mif+1,maf))
        output.add_title("Frame {0}".format(f+mif+1),3)
        output.add_text("{} clusters, {} Loners detected".format(
            experiment.groups[f]["N"],
            experiment.groups[f]["loners"]))
        output.add_fig("particle_{0}".format(f+mif+1),
                       visual.plot_particle,
                       (experiment.points[f+mif],
                        experiment.attribution[f],
                        experiment.X,
                        experiment.Y,
                        experiment.groups[f],
                        experiment.groups[f]["cbs"]),
                       proportions=(2*float(experiment.X)/float(experiment.Y),2))

        output.add_fig("creach_{0}".format(f+mif+1),
                       visual.plot_clust,
                       (experiment.reach[f],
                        experiment.clust[f][1],
                        experiment.clust[f][0]),
                       proportions=(2*float(experiment.X)/float(experiment.Y),2))

        output.add_fig("distrib_{0}".format(f),
                       visual.distrib,(experiment.distrib[f]["crowding"],))


        output.add_text("""
        <strong>Mean crowding</strong> :   {0[mean_crowding]},
        <strong>Mean group size</strong> : {0[mean_gs]} 
        """.format(experiment.distrib[f]))


    # ##  export 
    output.export()
    return output,experiment

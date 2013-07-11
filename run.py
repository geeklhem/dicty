#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""Image analysis"""

import export
import visual
import data 
import group_detection
import partitions
import analysis


def mock():
    # Load data and create export
    output = export.HtmlExport("mock_export")
    experiment = data.Data("stack.csv")
    experiment.groups = group_detection.from_file("centers35.csv")
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


def opt_exp(mif=0,maf=20,thres=400):
    # Load data and create export
    output = export.HtmlExport("optics")
    print("Loading...")
    experiment = data.Data("stack.csv")
    experiment.attribution = []
    experiment.reach = []
                        
    print("Analysis...")
    # ANALYSIS
    for f in range(mif,maf):
        print("{}/{}".format(f+1,maf))
        attr,reach = partitions.optics(experiment.points[f], thres)
        experiment.reach.append(reach)
        experiment.attribution.append(attr)

    # EXPORT 
    print("Export")
    output.add_text("Analysis of a film")


    output.add_title("Global")
    output.add_text("{} frames".format(experiment.frame_nb))

    output.add_title("Frame by frame".format(f))
    for f in range(maf-mif):
        print("{}/{}".format(f+mif+1,maf))
        output.add_title("Frame {0}".format(f+mif+1),3)
        output.add_fig("particle_{0}".format(f+mif+1),
                       visual.plot_particle,
                       (experiment.points[f+mif],
                        experiment.attribution[f],
                        experiment.X,
                        experiment.Y),
                       proportions=(2*float(experiment.X)/float(experiment.Y),2))

    ##  export 
    output.export()
    return output,experiment

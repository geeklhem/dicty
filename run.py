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
import os.path

class Experiment(object):
    """A run of this software"""
    def __init__(self,name,datafile="data/stack.csv",force=False,options={}):
        """Experiment constructor"""

        self.path = os.path.join("exports/",name+"/")
        self.savefile = os.path.join(self.path,name+".dicty")

        # Write if file doesn't exist, overwrite if force option enabled.
        if not os.path.isfile(self.savefile) or force:
            self.output = export.HtmlExport(name)
            self.data = data.Data(datafile)
            self.run(**options)
        else:
            self.output,self.data = self.load()
        

        # Run analysis, export result and save the object.
        self.output.export()
        self.save((self.output,self.data))

        def load(self):
            with open(self.savefile, 'rb') as fichier:
                unpickler = pickle.Unpickler(fichier)
                return unpickler.load()
                
        def save(self,s):
            with open(self.savefile, 'wb') as fichier:
                pickler = pickle.Pickler(fichier)
                pickler.dump(s)
    
        def run(self,**options):
            pass

    
class VoronoiAnalysis(Self.Data):
    """A voronoi analysis of the data"""
    def run(self,mf=None):
    # Load data and create export
        self.data.groups = group_detection.from_file("data/centers35.csv")
        self.data.attribution = []
        self.data.cbs = []
        self.data.distrib = []
        self.data.corr_cbs_ls = []
        self.data.corr_cl_ls = []

        self.data.attr_loners = (partitions.voronoi(self.data.groups["pos"],
                                                     self.data.points[-1]))
        self.data.loners = analysis.cell_by_section(self.data.attr_loners,
                                                     self.data.groups["N"])



        if not maf:
            maf = self.data.frame_nb


        # ANALYSIS
        for f in range(mf):
            self.data.attribution.append(partitions.voronoi(self.data.groups["pos"],
                                                             self.data.points[f]))
            self.data.cbs.append(analysis.cell_by_section(self.data.attribution[f],
                                                           self.data.groups["N"]))
            self.data.distrib.append(analysis.group_size_distrib(self.data.cbs[f],
                                                                  self.data.loners,
                                                                  self.data.groups["N"]))
            self.data.corr_cbs_ls.append(analysis.lin_correlate(self.data.cbs[f],
                                                                 self.data.groups["area"]))
            self.data.corr_cl_ls.append(analysis.lin_correlate(self.data.cbs[f],
                                                                self.data.loners))




        self.data.cc_cbs_ls = [c["corrcoef"] for c in self.data.corr_cbs_ls]
        self.data.cc_cl_ls = [c["corrcoef"] for c in self.data.corr_cl_ls]


        # EXPORT 
        self.output.add_text("Analysis of a film")


        self.output.add_title("Global")
        self.output.add_text("{} frames".format(self.data.frame_nb))

        self.output.add_fig("cbsls",
                       visual.time,(self.data.cc_cbs_ls,"Correlation",1))

        self.output.add_fig("clls",
                       visual.time,(self.data.cc_cl_ls,"Correlation",1))

        self.output.add_fig("cbs",
                       visual.areaplot,(self.data.cbs,),
                       proportions=(3,1))


        self.output.add_title("Frame by frame".format(f))
        for f in range(mf):
            self.output.add_title("Frame {0}".format(f),3)
            self.output.add_fig("particle_{0}".format(f),
                           visual.plot_particle,(self.data.points[f],
                                                 self.data.attribution[f],
                                                 self.data.X,
                                                 self.data.Y,
                                                 self.data.groups,
                                                 self.data.cbs[f],
                                             ),
                           proportions=(2*float(self.data.X)/float(self.data.Y),2))

            self.output.add_fig("distrib_{0}".format(f),
                           visual.distrib,(self.data.distrib[f]["crowding"],))

            self.output.add_fig("cbs_corr_{0}".format(f),
                           visual.correlation,(self.data.cbs[f],
                                               self.data.groups["area"],                 
                                               self.data.corr_cbs_ls[f]["regline"],
                                               "Cell by section",
                                               "Group Area"))


            self.output.add_fig("cl_corr_{0}".format(f),
                           visual.correlation,(self.data.cbs[f],
                                               self.data.loners,                 
                                               self.data.corr_cl_ls[f]["regline"],
                                               "Cell by section",
                                               "Loners by sections"))

            self.output.add_text("""
            <strong>Mean crowding</strong> :   {0[mean_crowding]},
            <strong>Mean group size</strong> : {0[mean_gs]} 
            """.format(self.data.distrib[f]))

            self.output.add_text("""
            <strong>Linear Correlation : </strong> 
            Group Area =
            {0[a]} Cell by section + {0[b]} 
            <strong> Correlation coefficient :</strong> {0[corrcoef]} 
            <strong> Mean square error :</strong>  {0[mse]}) 
            """.format(self.data.corr_cbs_ls[f]))

            self.output.add_text("""
            <strong>Linear Correlation : </strong> 
            Loners by section =
            {0[a]} Cell by section + {0[b]} 
            <strong> Correlation coefficient :</strong> {0[corrcoef]}
            <strong> Mean square error</strong>  : {0[mse]}) 
            """.format(self.data.corr_cl_ls[f]))

  
class OpticsAnalysis(Self.Data):
    def run(self,mif=0,maf=None,M=15):
        # Load data and create export
        self.data.attribution = []
        self.data.reach = []
        self.data.groups = []
        self.data.order = []
        self.data.clust = []
        self.data.distrib = []

        if not maf:
            maf = self.data.frame_nb

        print("Analysis...")
        # ANALYSIS
        for f in range(mif,maf):
            print("{}/{}".format(f+1,maf))
            attr,reach,order,cluster,color = partitions.optics_clust(self.data.points[f],
                                                                     eps=9000,
                                                                     ksi=0.001,
                                                                     M=M)
            self.data.reach.append(reach)
            self.data.order.append(order)
            self.data.attribution.append(attr)
            self.data.groups.append(group_detection.from_attribution(
                self.data.points[f],
                attr,
                loners_index=0))
            self.data.clust.append((cluster,color))
            self.data.distrib.append(analysis.group_size_distrib(
                self.data.groups[-1]["cbs"],
                self.data.groups[-1]["loners"],
                self.data.groups[-1]["N"]))

        # EXPORT 
        print("Export")
        self.output.add_text("Analysis of {}".format(self.name))
        self.output.add_text("<strong>Clustering algorithm</strong> : Optics")
        self.output.add_text("<strong>Parameters</strong> : M = {}".format(M))

        self.output.add_title("Global")
        self.output.add_text("{} frames".format(self.data.frame_nb))

        self.output.add_title("Frame by frame".format(f))
        for f in range(maf-mif):
            print("{}/{}".format(f+mif+1,maf))
            self.output.add_title("Frame {0}".format(f+mif+1),3)
            self.output.add_text("{} clusters, {} Loners detected".format(
                self.data.groups[f]["N"],
                self.data.groups[f]["loners"]))
            self.output.add_fig("particle_{0}".format(f+mif+1),
                           visual.plot_particle,
                           (self.data.points[f+mif],
                            self.data.attribution[f],
                            self.data.X,
                            self.data.Y,
                            self.data.groups[f],
                            self.data.groups[f]["cbs"]),
                           proportions=(2*float(self.data.X)/float(self.data.Y),2))

            self.output.add_fig("creach_{0}".format(f+mif+1),
                           visual.plot_clust,
                           (self.data.reach[f],
                            self.data.clust[f][1],
                            self.data.clust[f][0]),
                           proportions=(2*float(self.data.X)/float(self.data.Y),2))

            self.output.add_fig("distrib_{0}".format(f),
                           visual.distrib,(self.data.distrib[f]["crowding"],))


            self.output.add_text("""
            <strong>Mean crowding</strong> :   {0[mean_crowding]},
            <strong>Mean group size</strong> : {0[mean_gs]} 
            """.format(self.data.distrib[f]))

#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
""" Experiment.py
Contain the experiment class, which create object interfacing the data and the analyses
algorithms.
New kind of analysis should be described in subclasses of the general Experiment one.
"""

import os.path
import pickle
import copy

import export
import visual
import data 
import group_detection
import partitions
import analysis
import tracking
import traceback

class Experiment(object):
    """ The experiment object handle the data and its analysis.
 
    - It comes with saving and loading functions (using pickle) for 
      easy recovering 
    - The kind of analysis performed is precised in the run function.
    In this generic class nothing is done. You have to create inherited
    classes of it to precise it. """

    def __init__(self,name,datafile="data/stack.csv",force=False,options={}):
        """Experiment constructor
        Load data, initialize the attributes and start the analysis.

        :param name: The name of the experiment. Export is under exports/name.
        :type name: str
        :param datafile: Path to the csv file containing particles positions.
        :type datafile: str
        :param force: If true overwrite previously created files.
        :type force: bool
        :param options: Parameters to the run function.
        :type options: dict
        """

        self.options = options

        if name == "Same":
            # If no name is precised, the saved file will have the same name
            # as the datafile.
            self.name = os.path.splitext(os.path.basename(datafile))[0]
        else:
            self.name = name

        self.path = os.path.join("exports/",name+"/")
        self.savefile = os.path.join(self.path,name+".dicty")

        # Write if file doesn't exist, overwrite if force option enabled.
        if not os.path.isfile(self.savefile) or force:
            print("Creating new file...")
            self.output = export.HtmlExport(name)
            #self.data = data.Data(datafile,fileformat="sim",interval=50)
            self.data = data.Data(datafile)
            self.save("Loading")
        else:
            print("Loading")
            self.output,self.data,self.saved_step = self.load()
            print("Data loaded, the file was recovered at step {}".format(self.saved_step))
        
        # Run analysis, export results and save the object.
        self.run(**options)
        
        
    def load(self):
        """ Loader for this object, using pickle."""
        with open(self.savefile, 'rb') as fichier:
            unpickler = pickle.Unpickler(fichier)
            return unpickler.load()

    def save(self,saved_step):
        """ Save function for this object, using pickle.
        :param saved_step: The step at which the save is made,
        for future loadings.
        :type saved_step: str
        """
        self.saved_step = saved_step
        with open(self.savefile, 'wb') as fichier:
            pickler = pickle.Pickler(fichier)
            pickler.dump((self.output,self.data,self.saved_step))

    def run(self,**options):
        """ Run analysis, export results and save the object. 
        Here nothing is done. Please subclass to precise the behavior
        :param options: Options for the analysis. Usually come from the -p flag.
        :type options: dict"""
        pass

  
class OpticsAnalysis(Experiment):
    """ Analyse the data using the optics clustering algorithm."""
    def run(self,mif=0,maf=None,M=15,step=None,eps=None):
        """ Run analysis, export results and save the object. 
            Save are done at each step. """
        
        if not maf:
            maf = self.data.frame_nb
        
        if self.saved_step == "Loading" or step == "Analysis":
            try:
                print("Analysis...")
                self.analysis(mif,maf,M,eps)
            except Exception:
                print("ERROR: Error in analysis")
                traceback.print_exc()
            else:
                self.save("Analysis")

        if self.saved_step == "Analysis" or step == "Tracking":
            try:
                print("Tracking...")
                self.tracking()
            except Exception:
                print("ERROR: Error in tracking")
                traceback.print_exc()
            else:
                self.save("Tracking")

        if self.saved_step == "Tracking" or step == "Export":
            try:
                print("Export...")
                self.export(mif,maf)
            except Exception:
                print("ERROR: Error in export")
                traceback.print_exc()
            else:
                self.save("Export")

    def analysis(self,mif,maf,M,eps):
        """ For each frame, execute the clustering algorithm, 
        extract clusters and compute their characteristics."""
        self.data.attribution = []
        self.data.reach = []
        self.data.order = []
        self.data.clusters = []
        self.data.distrib = []
        self.data.frame_info = []
        self.data.reach_color = []
        self.data.all_clusters= []
        self.data.ofs = []

        # ANALYSIS
        for f in range(mif,maf):
            print("Frame : {}/{}".format(f+1,maf))
            partition = partitions.optics_clust(self.data.points[f],
                                                eps=eps,
                                                M=M)
            self.data.frame_info.append({
                "loners":partition["loners"],
                "N":len(self.data.points[f][0:][0]),})

            self.data.reach.append(partition["reach"])
            self.data.order.append(partition["order"])
            self.data.ofs.append(partition["ofs"])
            self.data.attribution.append(partition["attribution"])
            self.data.clusters.append([c for c in partition["clusters"] if c["leaf"]])
            self.data.reach_color.append(partition["color_histo"])
            self.data.all_clusters.append(partition["clusters"])

            self.data.distrib.append(analysis.group_size_distrib(
                [c["N"] for c in self.data.clusters[-1]],
                self.data.frame_info[-1]["loners"],
                self.data.frame_info[-1]["N"]))


    def tracking(self):
        """Once each frame has been processed, reconstruct the genealogy of groups"""

        # Find the ancestor of each cluster.
        self.data.tracking = tracking.track_cluster(self.data.clusters)
        
        # Reconstruct the traces
        self.data.traces = tracking.traces(self.data.tracking) 
        
        # Get the index of the tracked clusters
        self.data.tclust = tracking.traced_clusters_list(self.data.traces,
                                                         self.data.clusters)

    def export(self,mif,maf):
        """ Construct the html export""" 
        self.output.elements = []
        
        # General figures
        self.output.add_text("Analysis of {}".format(self.name))
        self.output.add_text("<strong>Clustering algorithm</strong> : Optics")
        self.output.add_text("<strong>Parameters</strong> : ".format(self.options))

        self.output.add_title("Global")
        self.output.add_text("{} frames".format(self.data.frame_nb))

        self.output.add_fig("tracking",
                            visual.tree_trace,
                            (self.data.tclust,
                             self.data.clusters),
                            proportions=(10,5))

        self.output.elements.append("<slideshow>")


        self.output.add_fig("Particles",
                            visual.time,
                            ([f["N"] for f in self.data.frame_info],
                             "Particles"),
                            proportions=(2,1))

        self.output.add_fig("Loners",
                            visual.time,
                            ([f["loners"] for f in self.data.frame_info],
                             "Loners"),
                            proportions=(2,1))

        self.output.add_fig("Clusters",
                            visual.time,
                            ([len(f) for f in self.data.clusters],
                             "Clusters"),
                            proportions=(4,2))

        self.output.add_fig("Mean group size",
                            visual.time,
                            ([f["mean_gs"] for f in self.data.distrib],
                             "Mean group size"),
                            proportions=(4,2))

        self.output.add_fig("Mean crowding",
                            visual.time,
                            ([f["mean_crowding"] for f in self.data.distrib],
                             "Mean crowding"),
                            proportions=(4,2))


        # Frame by frame figures
        self.output.add_title("Frame by frame".format(f))
        for f in range(maf-mif):
            print("Frame {}/{}".format(f+mif+1,maf))
            self.output.add_title("Frame {0}".format(f+mif+1),3)
            self.output.add_text("{} clusters, {} Loners detected".format(
                self.data.frame_info[f]["N"],
                self.data.frame_info[f]["loners"]))
            self.output.add_fig("particle_{0:02}".format(f+mif+1),
                           visual.plot_particle,
                           (self.data.points[f+mif],
                            self.data.attribution[f],
                            self.data.X,
                            self.data.Y,
                            self.data.clusters[f]),
                           proportions=(2*float(self.data.X)/float(self.data.Y),2))

            self.output.add_fig("creach_{0:02}".format(f+mif+1),
                           visual.plot_clust,
                           (self.data.reach[f],
                            self.data.reach_color[f],
                            self.data.all_clusters[f]),
                           proportions=(2*float(self.data.X)/float(self.data.Y),2))

            self.output.add_fig("distrib_{0:02}".format(f+mif+1),
                           visual.distrib,(self.data.distrib[f]["crowding"],))


            self.output.add_text("""
            <strong>Mean crowding</strong> :   {0[mean_crowding]},
            <strong>Mean group size</strong> : {0[mean_gs]} 
            """.format(self.data.distrib[f]))
       
        self.output.add_slideshow("particle_")

        # Make html
        self.output.export()

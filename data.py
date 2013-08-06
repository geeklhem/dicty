#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""Loading data"""

import numpy as np

class Data(object):
    """Hold data about a film """
    def __init__(self,pointsfile,csl=1000,fileformat=None,interval=None):
        if fileformat == "sim":
            self.load_data_sim(pointsfile,interval)
        else:
            self.load_data(pointsfile,csl)
        self.X = max([max(p[0,:]) for p in self.points]) 
        self.Y = max([max(p[1,:]) for p in self.points]) 
        self.X *= 1.1
        self.Y *= 1.1

    def load_data(self,pointsfile,csl):
        """ Load a csv file containing for each particle : "area","x","y","frame"."""
        images = np.genfromtxt(pointsfile,
                               skip_header=1,
                               delimiter=",",
                               usecols=(1,2,3,4),
                               names=("area","x","y","img"))
        

        # Points for successives images
        self.points = []
        self.frame_names = []
        self.N = []

        for img in range(1,int(images["img"][-1])+1):
            self.points.append((
                np.transpose(np.array([(x,y) 
                                       for ar,x,y,sl 
                                       in images
                                       if (ar < csl and
                                           sl == img)]))))
            self.frame_names.append("slice_{}".format(img))
            self.N.append(len(self.points[-1]))
        
        self.frame_nb = len(self.frame_names)

    def load_data_sim(self,pointsfile,interval):
        print("Load text file")
        images = np.genfromtxt(pointsfile,
                               delimiter=",",
                               usecols=(0,1,3),
                               names=("x","y","img"))

        # Points for successives images
        self.points = []
        self.frame_names = []
        self.N = []

        for img in range(1,int(images["img"][-1])/interval):
            print("Load {}/{}".format(img,int(images["img"][-1])/interval))
            img = img*interval
            self.points.append((
                np.transpose(np.array([(x,y) 
                                       for x,y,sl 
                                       in images
                                       if  sl == img]))))
            self.frame_names.append("slice_{}".format(img))
            self.N.append(len(self.points[-1]))
        
        self.frame_nb = len(self.frame_names)

        

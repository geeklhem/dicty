#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""Loading data"""

import numpy as np

class Data(object):
    """Hold data about a film """
    def __init__(self,pointsfile,csl=1000,
                    x=(9744*500/789),y=(7280*500/789)):
              self.X = x
              self.Y = y
              self.load_data(pointsfile,csl)
              
    def load_data(self,pointsfile,csl):
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

    def density(self,precision=100):
        """Compute the density of the particle and return a
        frame_nb(precision*precsion) list of array"""
        pass

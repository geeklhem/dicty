#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""Data Analysis"""

import collections
import numpy as np
import scipy as sp


def cell_by_section(attribution,nb_sections):
    cc = collections.Counter(attribution)

    cell_by_center = [0] * nb_sections
    for k,v in cc.iteritems():
        cell_by_center[k] = v
    return cell_by_center


def group_size_distrib(cellBySection,lonersBySection,nb_sections):
    g = [0] * (max(cellBySection)+1)
    g[1] = sum(lonersBySection)

    for c in cellBySection:
        g[c] += c
    gs = [float(i)/float(n) if n
          else 0
          for n,i in enumerate(g)]

    return {"group_size_distrib":gs,
            "crowding":g,
            "mean_gs":np.mean(gs),
            'mean_crowding':np.mean(g)}


def lin_correlate(x,y):
    (a,b)=sp.polyfit(x,y,1)
    line=sp.polyval([a,b],x)
    err=np.sqrt(sum((line-y)**2)/len(x))
    corr = np.corrcoef(x,y)[0,1]

    return {"a":a,"b":b,"regline":line,"mse":err,"corrcoef":corr}

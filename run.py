#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""
DICTy : Dictyostelium Images Clustering and Tracking for pYthon 

Usage:
  run.py <file> [-a=<analysis>] [-o=<outputName>] [-p=<parameters>] [-f]
  run.py [--help|-h|--version|--license]

Where <file> is a csv file containing for each particle : "area","x","y","frame".

Options:
  -a=<analysis>            Analysis to run [default: OpticsAnalysis]
  -o=<outputName>          Export name [default: Same] (Same as data basename)
  -p=<parameters>          Analysis parmeters in the format "p1=value;p2=value;p3=v"
  -f                       Force overwrite
  -h --help                Show this screen.
  --version                Show version.
  --license                Show license information.
"""

import sys
import ast

import docopt
import experiment 

__author__ = "Guilhem Doulcier"
__copyright__ = "Copyright 2013, Guilhem Doulcier"
__license__ = "GPLv3"
__version__ = "alpha"

def main(args):
    if  args["--license"]:
        print('\n  DICTy : Dictyostelium Images Clustering and Tracking for pYthon  - v'+__version__)
        print("    Copyright (C) 2013 Guilhem DOULCIER")
        print("    This program comes with ABSOLUTELY NO WARRANTY.")
        print("    This program is free software: you can redistribute it and/or modify")
        print("    it under the terms of the GNU General Public License as published by") 
        print("    the Free Software Foundation, either version 3 of the License, or") 
        print("    (at your option) any later version.\n")
        sys.exit(2)
    
    param = {}
    if args["-p"]:
        custom_p = args["-p"].split(",")
        for p in custom_p:
            p = p.split("=")
            try:
                p[1] = ast.literal_eval(p[1])
            except:
                print("Error in parsing {0} argument".format(p[0]))
            else:
                param[p[0]]=p[1]
    
    ExperimentClass = getattr(experiment,args["-a"])
    
    exp = ExperimentClass(args["-o"],
                          datafile=args["<file>"],
                          force=args["-f"],
                          options=param)

    return exp

if __name__ == "__main__":
    args = docopt.docopt(__doc__, version=__version__)
    main(args)
    

fake_args = {'--help': False,
             '--license': False,
             '--version': False,
             '-a': 'OpticsAnalysis',
             '-f': True,
             '-o': 'Same',
             '-p': "mif=15,maf=18",
             '<file>': "data/stack.csv"}
    
argtot = {'--help': False,
             '--license': False,
             '--version': False,
             '-a': 'OpticsAnalysis',
             '-f': True,
             '-o': 'Same',
             '-p': "",
             '<file>': "data/stack.csv"}
    


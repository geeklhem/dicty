#!/usr/bin/env/ python
# -*- coding: utf-8 -*-
"""HTML export"""

import os
import errno
import time 
import socket
import matplotlib.pyplot as plt

class HtmlExport(object):
    def __init__(self,name):
        self.date = time.strftime("%Y-%m-%d %H:%M:%S",time.gmtime())
        self.host = socket.gethostname()
       
        self.header = """
        <html>
        <head>
        <title>Experimental picture analysis</title>
        <style type="text/css">
        img {max-width:100%;}
        footer {color:grey;}
        </style>
        </head><body>
        <h1>Experimental picture analysis</h1>
        """
        self.footer = """
        <footer>
        Created {date} on {host} 
        </footer>
        </body>""".format(date=self.date,host=self.host)
        self.elements = []
        self.name = name
        self.path = os.path.join("exports",self.name)
        try:
            os.mkdir(self.path)
        except OSError as ex:
            if ex.errno != errno.EEXIST:
                raise
    def add_all_img(self,pattern):
        pics = ""
        for i in sorted(glob.glob(os.path.join(name,pattern))):
               pics += """
               <img src ="{path}" alt="figure matching '{pattern}'"/>
               """.format(path=os.path.basename(i))
        self.elements.append(pics)

    def add_text(self,text):
        self.elements.append("<p>{}</p>".format(text))

    def add_title(self,text,lvl=2):
        self.elements.append("<h{l}>{title}</h{l}>".format(title=text,l=lvl))

    def add_fig(self,name,graphical_function,args=(),kargs={},proportions=(1,1)):
       
        p = os.path.join(self.path,name+".png")

        graphical_function(*args,show=False,**kargs)
        f = plt.gcf()
        f.set_dpi(150)
        f.set_size_inches((5*proportions[0],5*proportions[1]))
        plt.savefig(p,bbox_inches="tight")
        plt.clf()

        self.elements.append("""
        <img src ="{path}" alt="plot generated by {name}"/>
        """.format(path=name+'.png',name=graphical_function.__name__))
        

    def export(self):
        page = self.header + "".join(self.elements) + self.footer
        with open(os.path.join(self.path,"index.html"), 'w') as f:
            f.write(page)

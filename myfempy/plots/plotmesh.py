# -*- coding: utf-8 -*-
"""
========================================================================
~~~ MODULO DE SIMULACAO ESTRUTURAL PELO METODO DOS ELEMENTOS FINITOS ~~~
       	                    __                                
       	 _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
       	| '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
       	| | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
       	|_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
       	            |___/                       |_|     |___/ 

~~~      Mechanical studY with Finite Element Method in PYthon       ~~~
~~~                PROGRAMA DE ANÁLISE COMPUTACIONAL                 ~~~
~~~              copyright @ 2022, all rights reserved               ~~~
========================================================================
"""

# import numpy as np
import vedo as vd 
from myfempy.tools.tools import get_version, get_logo



# @profile
def post_show_mesh(file2plot, plotset):
    
        win = vd.Plotter(title='POST-PROCESS', sharecam=False,  screensize=(1280, 720))
        
        mesh = vd.Mesh(file2plot+'.vtk').lineWidth(1).flat()

        if plotset['edge'] == False:
           mesh = vd.Mesh(file2plot+'.vtk')
        else:
            pass
        
        cname = vd.colorMap(range(21), 'jet')
        mesh.cmap(cname, on=plotset['apply']).addScalarBar(title=plotset['text_plot'],c='w')
        text = vd.Text2D('MYFEMPY v'+get_version()+' < '+plotset['text_plot']+' >\nPress "q" to continue...',  s = 1, font = 'Arial', c= 'white')
        win.show(text, mesh, viewup='y', bg='black', axes = 4)
        # win.close()
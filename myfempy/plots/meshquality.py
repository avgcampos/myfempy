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

import vedo as vd 
import numpy as np
from myfempy.tools.tools import get_version, get_logo

class MeshProp:
    
       
    def __init__(self, plotset):
        
        self.plotset = plotset
    
    def mesh_numbering(self):

        win = vd.Plotter(title='PRE-PROCESS', sharecam=False,  screensize=(1280, 720))
    
        mesh = vd.Mesh(self.plotset['RENDER']['filename']+'.vtk').lineWidth(0.1).flat()
        if self.plotset['LABELS']['lines'] == False:
           mesh = vd.Mesh(self.plotset['RENDER']['filename']+'.vtk')
        else:
            pass
        
        mesh.cmap("RdYlBu", on='cells', n=4)
        text = vd.Text2D('MYFEMPY v'+get_version()+' < mesh numb. > ',  s = 1, font = 'Arial', c= 'white')
        
        
        nodes = self.plotset['inci'][:,4:4+self.plotset["nodecon"]].reshape((self.plotset["nnode"]*self.plotset["nodecon"],))
        noduni, idx = np.unique(nodes, return_index=True)
        
        
        labs0 = mesh.labels(content=nodes[np.sort(idx)].astype(int), # 'id'
                        cells=False,
                        scale=self.plotset['LABELS']['scale'],
                        font='Arial',
                        c='white')
        
            
        labs1 = mesh.labels('cellid',#content=self.plotset['inci'][:,0].astype(int), # 'cellid'
                        cells=True,
                        scale=self.plotset['LABELS']['scale'],
                        font='Arial',
                        c='green')
              
        win.show(text, mesh, labs0, viewup='y', bg='black', axes=4)
        # win.render()
        # win.clear()
        # show(text, mesh, labs, viewup='y', bg='white', axes=4, screensize=(1280, 720))
        # vd.close()
            
                    
    def  mesh_quality(self):
        
        win = vd.Plotter(title='PRE-PROCESS', sharecam=False,  screensize=(1280, 720))
        
        mesh = vd.Mesh(self.plotset['RENDER']['filename']+'.vtk').lineWidth(0.1).flat()
            
        if self.plotset['QUALITY']['lines'] == False:
           mesh = vd.Mesh(self.plotset['RENDER']['filename']+'.vtk')
        else:
            pass
        
        mesh.cmap("RdYlBu", on='cells', n=16).addScalarBar()
        text = vd.Text2D('MYFEMPY v'+get_version()+' < mesh quality > ',  s = 1, font = 'Arial', c = 'white')
        
        mesh.addQuality(measure=self.plotset['QUALITY']['method']).cmap('RdYlBu', on='cells').print()
        # hist = vd.pyplot.histogram(mesh.celldata["Quality"], xtitle='mesh quality', bc='black')
        # hist.scale(self.plotset['quality']["scale"]).pos(self.plotset['quality']['posx'],self.plotset['quality']['posy']).useBounds(True)
        mesh.addScalarBar3D(c='white', title=('Meth. '+str(self.plotset['QUALITY']['method'])))
        
        # labs = mesh.labels(cells=True,
        #                    precision=3,
        #                    scale=self.plotset['quality']['scale'],
        #                    font='Arial',
        #                    c='black',
        #                    )
                
        win.show(text, mesh, viewup='y', bg='black', axes = 4)
        # win.render()
        # win.clear()
        # win.close()
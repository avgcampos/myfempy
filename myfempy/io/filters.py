# -*- coding: utf-8 -*-
"""
_______________________________________________________________________________
~~~~~~                        MYFEMPY                                   ~~~~~~~
>> innputfilter

===============================================================================
@author: ANTONIO VINICIUS GARCIA CAMPOS
@copyright: 3D EASYCAE SERVIÇOS DE ANÁLISE COMPUTACIONAL, 2022
@licence: GPL-3.0 License
"""

#-----------------------------------------------------------------------------#
import sys
import os
import imp
import numpy as np
from colorama import Fore, Back, Style

#-----------------------------------------------------------------------------#
def input_filter(path_user,file_name):
    pointlist = np.zeros((1,4))
    linelist = np.zeros((1,3))
    propgeonewlist = np.zeros((1,6))
    propgeobiblist = np.zeros((1,6))
    compmaterial = np.zeros((1,4))
    propmatlist = np.zeros((1,8))
    solutconfig = np.zeros((1,3))
    meshconfig = np.zeros((1,5))
    forcelist = np.zeros((1,9))
    boncdlist = np.zeros((1,7))
    solvercfg = np.zeros((1,7))
    outputcfg = np.zeros((1,4))
        
    file_input = str(path_user+'/'+file_name)
    
    with open(file_input,'r') as file_object:
        line = file_object.readline()
        line = file_object.readline()
        line = file_object.readline()
        line = file_object.readline()
        line = file_object.readline()
        line = file_object.readline()
        line = file_object.readline()              
        if line == '#SOLVER_CFG\n':
            line = file_object.readline() 
            while line != "\n":
                lineaux = line.split()
                linearray = np.array([lineaux[0],lineaux[2],lineaux[3],lineaux[4],lineaux[5],lineaux[6],lineaux[7]])
                solvercfg = np.append(solvercfg,[linearray],axis=0)
                line = file_object.readline()
                

        line = file_object.readline()
        if line == '#IMPORT_GEO\n':
            line = file_object.readline()
        
        line = file_object.readline()
        if  line == '#NEW_GEO\n':
            line = file_object.readline()
            if  line == '##POINT_COORD\n':
                line = file_object.readline()
                while line != "\n":
                    lineaux = line.split()
                    linearray = np.array([lineaux[0],lineaux[2],lineaux[3],lineaux[4]])
                    pointlist = np.append(pointlist,[linearray],axis=0)
                    line = file_object.readline()
            
            line = file_object.readline()
            if  line == '##LINE_CONEC\n':
                line = file_object.readline()
                while line != "\n":
                    lineaux = line.split()
                    linearray = np.array([lineaux[0],lineaux[2],lineaux[3]])
                    linelist = np.append(linelist,[linearray],axis=0)
                    line = file_object.readline()
                
            line = file_object.readline()
            planelist = []
            if line == '##PLANE_CONEC\n':
                line = file_object.readline()
                while line != "\n":
                    lineaux = line.split()
                    planelist.append(lineaux[2:])
                    line = file_object.readline()

                
        line = file_object.readline()              
        if  line == '#PROP_GEOMETRIA\n':
            line = file_object.readline()
            if  line == '##GEO_INP\n':
                line = file_object.readline()
                while line != "\n":
                    lineaux = line.split()
                    linearray = np.array([lineaux[0],lineaux[2],lineaux[3],lineaux[4],lineaux[5],lineaux[6]])
                    propgeonewlist = np.append(propgeonewlist,[linearray],axis=0)
                    line = file_object.readline()
            
            line = file_object.readline()
            if  line == '##GEO_BIB\n':
                line = file_object.readline()
                while line != "\n":
                    lineaux = line.split()
                    linearray = np.array([lineaux[0],lineaux[2],lineaux[3],lineaux[4],lineaux[5],lineaux[6]])
                    propgeobiblist = np.append(propgeobiblist,[linearray],axis=0)
                    line = file_object.readline()
        
        line = file_object.readline()              
        if  line == '#PROP_MATERIAL\n':
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                linearray = np.array([lineaux[0],lineaux[2],lineaux[3],lineaux[4]])
                compmaterial = np.append(compmaterial,[linearray],axis=0)
                line = file_object.readline()
        
        line = file_object.readline()              
        if  line == '##MAT_CFG\n':
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                linearray = np.array([lineaux[0],lineaux[2],lineaux[3],lineaux[4],lineaux[5],lineaux[6],lineaux[7],lineaux[8]])
                propmatlist = np.append(propmatlist,[linearray],axis=0)
                line = file_object.readline()
                
        line = file_object.readline()
        if  line == '#MODEL_CFG\n':
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                linearray = np.array([lineaux[0],lineaux[2],lineaux[3]])
                solutconfig = np.append(solutconfig,[linearray],axis=0)
                line = file_object.readline()

        line = file_object.readline()              
        if  line == '#MESH_CFG\n':
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                linearray = np.array([lineaux[0],lineaux[2],lineaux[3],lineaux[4],lineaux[5]])
                meshconfig = np.append(meshconfig,[linearray],axis=0)
                line = file_object.readline()


        line = file_object.readline()              
        if  line == '#SOLUTION_CFG\n':
            line = file_object.readline()
            if  line == '##LOADS\n':
                line = file_object.readline()
                while line != "\n":
                    lineaux = line.split()
                    linearray = np.array([lineaux[0],lineaux[2],lineaux[3],lineaux[4],lineaux[5],lineaux[6],lineaux[7],lineaux[8],lineaux[9]])
                    forcelist = np.append(forcelist,[linearray],axis=0)
                    line = file_object.readline()
            
            line = file_object.readline()
            if  line == '##BOND_COND\n':
                line = file_object.readline()
                while line != "\n":
                    lineaux = line.split()
                    linearray = np.array([lineaux[0],lineaux[2],lineaux[3],lineaux[4],lineaux[5],lineaux[6],lineaux[7]])
                    boncdlist = np.append(boncdlist,[linearray],axis=0)
                    line = file_object.readline()
            
        line = file_object.readline()              
        if  line == '#VIEWSOLUTION_CFG\n':
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                linearray = np.array([lineaux[0],lineaux[2],lineaux[3],lineaux[4]])
                outputcfg = np.append(outputcfg,[linearray],axis=0)
                line = file_object.readline()
                
        line = file_object.readline()
        if  line == '#GRAPH_OUT\n':
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                graphout = np.zeros((1,1))
                graphout = np.append(graphout,np.array([[lineaux[0]]]),axis=1)
                graphout = np.append(graphout,np.array([lineaux[2:]]),axis=1)
                line = file_object.readline()
        
        line = file_object.readline()
        if  line == '#FILE_OUT\n':
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                fileout = np.zeros((1,1))
                fileout = np.append(fileout,np.array([[lineaux[0]]]),axis=1)
                fileout = np.append(fileout,np.array([lineaux[2:]]),axis=1)
                line = file_object.readline()
                    
                
    pointlist = pointlist[1::][::]
    linelist = linelist[1::][::]
    # planelist = planelist[1::][::]
    propgeonewlist = propgeonewlist[1::][::]  
    propgeobiblist = propgeobiblist[1::][::] 
    propmatlist = propmatlist[1::][::]
    compmaterial = compmaterial[1::][::]
    solutconfig = solutconfig[1::][::]
    meshconfig = meshconfig[1::][::]
    forcelist = forcelist[1::][::]
    boncdlist = boncdlist[1::][::]
    solvercfg = solvercfg[1::][::]
    outputcfg = outputcfg[1::][::]
    graphout = graphout[0,1:]
    fileout = fileout[0,1:]
    
    return solutconfig,solvercfg,pointlist,linelist,planelist,propgeonewlist,propgeobiblist,compmaterial,propmatlist,meshconfig,forcelist,boncdlist,solvercfg,outputcfg,graphout,fileout

#-----------------------------------------------------------------------------#





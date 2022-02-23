# -*- coding: utf-8 -*-
"""
_______________________________________________________________________________
~~~~~~                        MYFEMPY                                   ~~~~~~~
>> ioctrl

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
def export_vtkmesh(postdata,datamesh,coord,inci):
    if  postdata['vtkcelltype'] == 3: # VTK_LINE
        with open(postdata['filename']+'.vtk','w') as file_object:
            file_object.write('# vtk DataFile Version 2.0\n')
            file_object.write('mesh vtk view\n')
            file_object.write('ASCII\n')
            file_object.write('DATASET UNSTRUCTURED_GRID\n')
            file_object.write('POINTS '+str(datamesh['lencoord'])+' double\n')
            for ii in range(0,datamesh['lencoord']):
                file_object.write(str(coord[ii,1])+' '+str(coord[ii,2])+' '+str(coord[ii,3])+'\n')
            
            file_object.write('\n')
            file_object.write('CELLS'+' '+str(datamesh['leninci'])+' '+str(3*datamesh['leninci'])+'\n')
            for ii in range(0,datamesh['leninci']):
                file_object.write('2'+' '+str(int(inci[ii,4])-1)+' '+str(int(inci[ii,5])-1)+'\n')
            
            file_object.write('\n')
            file_object.write('CELL_TYPES'+' '+str(datamesh['leninci'])+'\n')
            for ii in range(0,datamesh['leninci']): 
                file_object.write('3'+'\n')
        
        
        
            if postdata['datatype'] == 'elm':
                file_object.write('\n')
                file_object.write('CELL_DATA'+' '+str(datamesh['leninci'])+'\n')
                for jj in range(postdata['datalist'].shape[1]):
                    file_object.write('SCALARS '+str(postdata['datatitle'][jj])+' float 1\n')
                    file_object.write('LOOKUP_TABLE default\n')
                    for ii in range(0,datamesh['leninci']): 
                        file_object.write(str(postdata['datalist'][ii,jj])+'\n')
                    file_object.write('\n')
            
            elif postdata['datatype'] == 'avr':
                file_object.write('\n')
                file_object.write('POINT_DATA'+' '+str(datamesh['lencoord'])+'\n')
                for jj in range(postdata['datalist'].shape[1]):
                    file_object.write('SCALARS '+str(postdata['datatitle'][jj])+' float 1\n')
                    file_object.write('LOOKUP_TABLE default\n')
                    for ii in range(0,datamesh['lencoord']): 
                        file_object.write(str(postdata['datalist'][ii,jj])+'\n')
                    file_object.write('\n')
    
    elif postdata['vtkcelltype'] == 5: # VTK_TRIANGLE
        with open(postdata['filename']+'.vtk','w') as file_object:
            file_object.write('# vtk DataFile Version 2.0\n')
            file_object.write('mesh vtk view\n')
            file_object.write('ASCII\n')
            file_object.write('DATASET UNSTRUCTURED_GRID\n')
            file_object.write('POINTS '+str(datamesh['lencoord'])+' double\n')
            for ii in range(0,datamesh['lencoord']):
                file_object.write(str(coord[ii,1])+' '+str(coord[ii,2])+' '+str(coord[ii,3])+'\n')
            
            file_object.write('\n')
            file_object.write('CELLS'+' '+str(datamesh['leninci'])+' '+str(4*datamesh['leninci'])+'\n')
            for ii in range(0,datamesh['leninci']):
                file_object.write('3'+' '+str(int(inci[ii,4])-1)+' '+str(int(inci[ii,5])-1)+' '+str(int(inci[ii,6])-1)+'\n')
            
            file_object.write('\n')
            file_object.write('CELL_TYPES'+' '+str(datamesh['leninci'])+'\n')
            for ii in range(0,datamesh['leninci']): 
                file_object.write('5'+'\n')
        
            if postdata['datatype'] == 'elm':
                file_object.write('\n')
                file_object.write('CELL_DATA'+' '+str(datamesh['leninci'])+'\n')
                for jj in range(postdata['datalist'].shape[1]):
                    file_object.write('SCALARS '+str(postdata['datatitle'][jj])+' float 1\n')
                    file_object.write('LOOKUP_TABLE default\n')
                    for ii in range(0,datamesh['leninci']): 
                        file_object.write(str(postdata['datalist'][ii,jj])+'\n')
                    file_object.write('\n')
            
            elif postdata['datatype'] == 'avr':
                file_object.write('\n')
                file_object.write('POINT_DATA'+' '+str(datamesh['lencoord'])+'\n')
                for jj in range(postdata['datalist'].shape[1]):
                    file_object.write('SCALARS '+str(postdata['datatitle'][jj])+' float 1\n')
                    file_object.write('LOOKUP_TABLE default\n')
                    for ii in range(0,datamesh['lencoord']): 
                        file_object.write(str(postdata['datalist'][ii,jj])+'\n')  
                    file_object.write('\n')
    
    elif postdata['vtkcelltype'] == 9: # VTK_QUAD
        with open(postdata['filename']+'.vtk','w') as file_object:
            file_object.write('# vtk DataFile Version 2.0\n')
            file_object.write('mesh vtk view\n')
            file_object.write('ASCII\n')
            file_object.write('DATASET UNSTRUCTURED_GRID\n')
            file_object.write('POINTS '+str(datamesh['lencoord'])+' double\n')
            for ii in range(0,datamesh['lencoord']):
                file_object.write(str(coord[ii,1])+' '+str(coord[ii,2])+' '+str(coord[ii,3])+'\n')
            
            file_object.write('\n')
            file_object.write('CELLS'+' '+str(datamesh['leninci'])+' '+str(5*datamesh['leninci'])+'\n')
            for ii in range(0,datamesh['leninci']):
                file_object.write('4'+' '+str(int(inci[ii,4])-1)+' '+str(int(inci[ii,5])-1)+' '+str(int(inci[ii,6])-1)+' '+str(int(inci[ii,7])-1)+'\n')
            
            file_object.write('\n')
            file_object.write('CELL_TYPES'+' '+str(datamesh['leninci'])+'\n')
            for ii in range(0,datamesh['leninci']): 
                file_object.write('9'+'\n')
        
            if postdata['datatype'] == 'elm':
                file_object.write('\n')
                file_object.write('CELL_DATA'+' '+str(datamesh['leninci'])+'\n')
                for jj in range(postdata['datalist'].shape[1]):
                    file_object.write('SCALARS '+str(postdata['datatitle'][jj])+' float 1\n')
                    file_object.write('LOOKUP_TABLE default\n')
                    for ii in range(0,datamesh['leninci']): 
                        file_object.write(str(postdata['datalist'][ii,jj])+'\n')
                    file_object.write('\n')
            
            elif postdata['datatype'] == 'avr':
                file_object.write('\n')
                file_object.write('POINT_DATA'+' '+str(datamesh['lencoord'])+'\n')
                for jj in range(postdata['datalist'].shape[1]):
                    file_object.write('SCALARS '+str(postdata['datatitle'][jj])+' float 1\n')
                    file_object.write('LOOKUP_TABLE default\n')
                    for ii in range(0,datamesh['lencoord']): 
                        file_object.write(str(postdata['datalist'][ii,jj])+'\n')
                    file_object.write('\n')
                    

    elif postdata['vtkcelltype'] == 12: # VTK_HEXAHEDRON
        with open(postdata['filename']+'.vtk','w') as file_object:
            file_object.write('# vtk DataFile Version 2.0\n')
            file_object.write('mesh vtk view\n')
            file_object.write('ASCII\n')
            file_object.write('DATASET UNSTRUCTURED_GRID\n')
            file_object.write('POINTS '+str(datamesh['lencoord'])+' double\n')
            for ii in range(0,datamesh['lencoord']):
                file_object.write(str(coord[ii,1])+' '+str(coord[ii,2])+' '+str(coord[ii,3])+'\n')
            
            file_object.write('\n')
            file_object.write('CELLS'+' '+str(datamesh['leninci'])+' '+str(9*datamesh['leninci'])+'\n')
            for ii in range(0,datamesh['leninci']):
                file_object.write('8'+' '+str(int(inci[ii,4])-1)+' '+str(int(inci[ii,5])-1)+' '+str(int(inci[ii,6])-1)+' '+str(int(inci[ii,7])-1)+' '+str(int(inci[ii,8])-1)+' '+str(int(inci[ii,9])-1)+' '+str(int(inci[ii,10])-1)+' '+str(int(inci[ii,11])-1)+'\n')
            
            file_object.write('\n')
            file_object.write('CELL_TYPES'+' '+str(datamesh['leninci'])+'\n')
            for ii in range(0,datamesh['leninci']): 
                file_object.write('12'+'\n')
        
            if postdata['datatype'] == 'elm':
                file_object.write('\n')
                file_object.write('CELL_DATA'+' '+str(datamesh['leninci'])+'\n')
                for jj in range(postdata['datalist'].shape[1]):
                    file_object.write('SCALARS '+str(postdata['datatitle'][jj])+' float 1\n')
                    file_object.write('LOOKUP_TABLE default\n')
                    for ii in range(0,datamesh['leninci']): 
                        file_object.write(str(postdata['datalist'][ii,jj])+'\n')
                    file_object.write('\n')
            
            elif postdata['datatype'] == 'avr':
                file_object.write('\n')
                file_object.write('POINT_DATA'+' '+str(datamesh['lencoord'])+'\n')
                for jj in range(postdata['datalist'].shape[1]):
                    file_object.write('SCALARS '+str(postdata['datatitle'][jj])+' float 1\n')
                    file_object.write('LOOKUP_TABLE default\n')
                    for ii in range(0,datamesh['lencoord']): 
                        file_object.write(str(postdata['datalist'][ii,jj])+'\n')
                    file_object.write('\n')
                
                    
                
#-----------------------------------------------------------------------------#
def usrlog_out(usrlog_file,solutconfig,meshconfig,compmaterial,forcelist,boncdlist,solvercfg,outputcfg,graphout,fileout):
    with open(usrlog_file,'w') as file_object:
        
        for i in range(len(solvercfg)):
            file_object.write('solver_typ = '+"'"+solvercfg[i,1]+"'"+'\n')
            file_object.write('solver_def = '+"'"+solvercfg[i,2]+"'"+'\n')
            file_object.write('solver_opt = '+"'"+solvercfg[i,3]+"'"+'\n')
            file_object.write('solver_start = '+"'"+solvercfg[i,4]+"'"+'\n')
            file_object.write('solver_end = '+"'"+solvercfg[i,5]+"'"+'\n')
            file_object.write('solver_step = '+"'"+solvercfg[i,6]+"'"+'\n')
        
        for i in range(len(solutconfig)):
            file_object.write('mod_typ = '+"'"+solutconfig[i,1]+"'"+'\n')
            file_object.write('mod_opt = '+"'"+solutconfig[i,2]+"'"+'\n')
        
        for i in range( len(meshconfig)):
            file_object.write('mesh_typ = '+"'"+meshconfig[i,1]+"'"+'\n')
        
        for i in range(len(compmaterial)):
            file_object.write('mat_typ = '+"'"+compmaterial[i,1]+"'"+'\n')
            file_object.write('mat_opt = '+"'"+compmaterial[i,2]+"'"+'\n') 
            file_object.write('mat_def = '+"'"+compmaterial[i,3]+"'"+'\n') 
        
        for j in range(len(np.unique(forcelist[:,8]))):
            forceliststep = forcelist[np.where(forcelist[:,8]==str(j+1)),:][0]
            for i in range(len(forceliststep)):
                file_object.write('force_typ_'+str(i)+' ''= '+"'"+forceliststep[i,1]+"'"+'\n')
                file_object.write('force_opt_'+str(i)+' ''= '+"'"+forceliststep[i,4]+"'"+'\n')
                file_object.write('force_opt_dir'+str(i)+' ''= '+"'"+forceliststep[i,2]+"'"+'\n')
                file_object.write('force_step'+str(j)+' ''= '+"'"+forceliststep[i,8]+"'"+'\n')
        
        for i in range(len(boncdlist)):
            file_object.write('bc_typ_'+str(i)+' ''= '+"'"+boncdlist[i,1]+"'"+'\n')
            file_object.write('bc_opt_'+str(i)+' ''= '+"'"+boncdlist[i,3]+"'"+'\n')
            file_object.write('bc_opt_dir'+str(i)+' ''= '+"'"+boncdlist[i,2]+"'"+'\n')
        
        file_object.write('prev_sol = '+"'"+outputcfg[0,1]+"'"+'\n')
        file_object.write('graph_out = '+"'"+outputcfg[0,2]+"'"+'\n')
        file_object.write('save_fig = '+"'"+outputcfg[0,3]+"'"+'\n')
        
        for i in range(len(graphout)):
            file_object.write(f'graph_out_{i} = '+"'"+graphout[i]+"'"+'\n')
            
        for i in range(len(fileout)):
            file_object.write(f'file_out_{i} = '+"'"+fileout[i]+"'"+'\n')
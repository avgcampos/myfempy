# -*- coding: utf-8 -*-
"""
_______________________________________________________________________________
~~~~~~ MODULO DE SIMULACAO ESTRUTURAL PELO METODO DOS ELEMENTOS FINITOS ~~~~~~~
~~~~~~      Mechanical studY with Finite Element Method in PYthon       ~~~~~~~
                              _____   _____   __  __                 
          _ __ ___    _   _  |  ___| | ____| |  \/  |  _ __    _   _ 
         | '_ ` _ \  | | | | | |_    |  _|   | |\/| | | '_ \  | | | |
         | | | | | | | |_| | |  _|   | |___  | |  | | | |_) | | |_| |
         |_| |_| |_|  \__, | |_|     |_____| |_|  |_| | .__/   \__, |
                      |___/                           |_|      |___/ 

~~~~~~                   PROGRAMA DE PROPOSITO GERAL                    ~~~~~~~  
~~~~~~                  copyright all rights reserved                   ~~~~~~~
===============================================================================
@author: ANTONIO VINICIUS GARCIA CAMPOS
@copyright: 3D EASYCAE SERVIÇOS DE ANÁLISE COMPUTACIONAL
@licence: GPL-3.0 License
    
Copyright (c) 2021, 3D EASYCAE SERVIÇOS DE ANÁLISE COMPUTACIONAL
All rights reserved.
"""


#%% EXECUTANDO A ANALISE
import sys
import time
import imp
import numpy as np
from colorama import Fore, Back, Style

import warnings
warnings.filterwarnings("ignore")


#%% WELCOME AND INSTRUCTIONS
from myfempy.setup.myfempy_welcome import logo_myfempy, about_myfempy
logo_myfempy()
about_myfempy(welcome_true='false')


#%% IMPORT MYFEMPY PREPROCESS PACKs
from myfempy.setup.myfempy_preproc import create_user_dir,readInpdata,gen_GmshId2Geo_out,gen_UsrLog_out,mesh_type,loading_bar_v1,white_output_logofile


#%% CRIAR NOME E PASTA DE DADOS DA ANALISE
print(Style.RESET_ALL)
usr_analysi_name = input('>> DEFINA O NOME DA ANÁLISE: ')
file_dir = input('>> DEFINA UMA PASTA PARA SALVAR OS DADOS DA ANÁLISE: ')
path_user = create_user_dir('user/',file_dir)
# file_name = input('>> NOME DO ARQUIVO DE ENTRADA DA ANÁLISE: ')
# file_name = file_name+'.txt'
file_name = 'myfempy_input_cfg.txt'

print(Style.RESET_ALL)
print(Fore.GREEN + Style.BRIGHT+'\r************  P R E - P R O C E S S   L O A D I N G ************')

loading_bar_v1(0)
#%% CONFIGURACAO DA ANALISE
solutconfig,solvercfg,pointlist,linelist,planelist,propgeonewlist,propgeobiblist,compmaterial,propmatlist,meshconfig,forcelist,boncdlist,solvercfg,outputcfg,graphout,fileout  = readInpdata(path_user,file_name)
# USER LOG OUTPUT
gen_UsrLog_out('myfempy/setup/'+'usrlog.txt',solutconfig,meshconfig,compmaterial,forcelist,boncdlist,solvercfg,outputcfg,graphout,fileout)
usrlog = imp.load_source('usrlog','myfempy/setup/'+'usrlog.txt')

loading_bar_v1(10)
#%% GERAR ARQUIVO DE GEOMETRIA GMSH
gmshgeo_file = usr_analysi_name+'_usr_geo.geo'
gen_GmshId2Geo_out(path_user+'/'+gmshgeo_file,pointlist,linelist,planelist,meshconfig,propgeonewlist)

loading_bar_v1(20)
#%% GERAR ARQUIVO DE MALHA GMSH
from myfempy.solve.read_mesh import gen_GmshGeo2Msh_out

file_geo = path_user+'/'+gmshgeo_file
name_out_mesh = path_user+'/'+usr_analysi_name+'_usr_mesh.msh1'
geo_select = usrlog.mod_typ
gen_GmshGeo2Msh_out(file_geo,name_out_mesh,geo_select)

loading_bar_v1(40)
#%% DEFINICAO DA MALHA 
datamesh,coord,inci,type_elm = mesh_type(usrlog,name_out_mesh)

loading_bar_v1(60)
#%% CONFIG. MATERIAL E GEOMETRIA
if len(propgeonewlist) != 0:
    numgeolist = len(propgeonewlist)
    tabgeo = np.zeros((numgeolist,5))
    tabdimSection = []
    tabtypSection = []
    screen_beam_cs = 'false'
    for i in range(numgeolist):
        tabgeo[i,:] = np.array([[float(propgeonewlist[i,1]),float(propgeonewlist[i,2]),float(propgeonewlist[i,3]),float(propgeonewlist[i,4]),float(propgeonewlist[i,5])]])
    
    y_max = 1.0  
    y_min = -1.0
    z_max = 1.0 
    z_min = -1.0 
    r_max = 1.0
    CGcoord = np.array([y_max,y_min,z_max,z_min,r_max])

if len(propgeobiblist) != 0:
    from myfempy.lib.geometric_properties import geometryProp
    from myfempy.solve.read_mesh import import_mesh_tria
    
    numgeolist = len(propgeobiblist)
    tabgeo = np.zeros((numgeolist,5))
    tabdimSection = []
    tabtypSection = []
    screen_beam_cs = 'true'
    numPtsRef = 10
    for i in range(numgeolist):
        typSection = propgeobiblist[i,1]
        dimSection = np.array([float(propgeobiblist[i,2]),float(propgeobiblist[i,3]),float(propgeobiblist[i,4]),float(propgeobiblist[i,5])])
        A,Izz,Iyy,Jxx,CGcoord = geometryProp(typSection,dimSection)
        tabgeo[i,:] = np.array([[A,Izz,Iyy,Jxx,999]])
        tabdimSection.append(dimSection)
        tabtypSection.append(typSection)

nummatlist = len(propmatlist)
tabmat = np.zeros((nummatlist,7))
for i in range(nummatlist):
    tabmat[i,:] = np.array([[float(propmatlist[i,1]),float(propmatlist[i,2]),float(propmatlist[i,3]),\
                            float(propmatlist[i,4]),float(propmatlist[i,5]),float(propmatlist[i,6]),float(propmatlist[i,7])]])

loading_bar_v1(80)
#%% FORCAS 
from myfempy.lib.domain_physics import model_physics
forcenodeaply, boncdnodeaply = model_physics(usrlog,datamesh,coord,inci,forcelist,boncdlist,tabgeo)

loading_bar_v1(100)
#%% PREVIEW DA ANALISE
from myfempy.solve.read_mesh import export_mesh
from myfempy.solve.view_mesh import view_analysis, view_beam_crossSection
from myfempy.setup.myfempy_welcome import myfempy_ctrllVer
myfempy_version = myfempy_ctrllVer() 

color = (np.array(inci[:,3])[np.newaxis]).T
data_title = ['NODISPL']
file_dir = path_user+'/'+usr_analysi_name+'_mesh_view.vtk'
typeData = 'elm'
export_mesh(file_dir,type_elm,len(coord),len(inci),coord,inci,color,data_title,typeData)

file_savePNG=path_user+'/'+usr_analysi_name+'_frame.png'
file_saveSTL=path_user+'/'+usr_analysi_name+'_frame.stl'
project_name = usr_analysi_name
scala_view = 0.05*max(max(coord[:,1]),max(coord[:,2]),max(coord[:,3]))
title_win = 'STRUCTURE PLOTTER'


screen_on = 'false'
save_screen = 'false'
if usrlog.prev_sol == 'prev_sol_true':
    screen_on ='true'
if usrlog.save_fig == 'save_fig_true':
    save_screen = 'true'
    
view_analysis(file_dir,file_savePNG,file_saveSTL,screen_on,save_screen,inci,coord,forcenodeaply,boncdnodeaply,project_name,title_win,myfempy_version,scala_view,tabdimSection,tabtypSection,screen_beam_cs)

#%% DEFINICAO DA SOLUCAO
# VETORES DE FORCAS E CONDICOES DE CONTORNO
print(Style.RESET_ALL)
print(Fore.GREEN + Style.BRIGHT+'\r*********************  S O L U T I O N   C O M P U T I N G   *******************')
from myfempy.lib.domain_physics import bond_condi_phys,loads_phys

F = loads_phys(datamesh,forcenodeaply)
freedof, fixedof = bond_condi_phys(datamesh,boncdnodeaply)

print(Fore.GREEN + Style.BRIGHT+'\r*********************  S O L V I N G   E Q U A T I O N S   *******************')
from myfempy.solve.myfempy_solve import myfempySolve

KG, U, U2, U3, time_list = myfempySolve(usrlog,datamesh,coord,inci,tabmat,tabgeo,forcenodeaply,F,freedof)

#%% POS PROCESSAMENTO
print(Style.RESET_ALL)
print(Fore.GREEN + Style.BRIGHT+'\r********************  P O S T - P R O C E S S   L O A D I N G  ********************')
if usrlog.graph_out == 'plot_res_true':
    from myfempy.lib.myfempy_posproc import graphic_prosproc
    print(Style.RESET_ALL)
    print(Fore.YELLOW + Style.BRIGHT+'\r>>> SAVING SOLUTION ...')
    graphic_prosproc(path_user,usr_analysi_name,usrlog,graphout,fileout,time_list,datamesh,inci,coord,tabgeo,CGcoord,tabmat,type_elm,KG,U,U2,U3,screen_on,save_screen)
      

#%%
print(Style.RESET_ALL)
print('\n')
print('#===                  A N A L Y S I S   S U C C E S S F U L                  ===#')
print('#===                        M Y F E M P Y @ 2 0 2 1                          ===#')
print('#===                   T H A N K   Y O U   F O R   U S E !                   ===#')
print('=================================================================================')


#%%
# if __name__ == '__main__':
#     main()
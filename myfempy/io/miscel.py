# -*- coding: utf-8 -*-
"""
Created on Sat May 16 12:05:12 2020
@author: ANTONIO VINICIUS GARCIA CAMPOS
@version: v13
_______________________________________________________________________________
 ~~~~~~~~~~ MODULO DE SIMULACAO PELO METODO DOS ELEMENTOS FINITOS ~~~~~~~~~~

ESTE MODULO GERA E LER DADOS DE ARQUIVO DE MALHA DO GMSH NA EXT. msh1 [version 1 
legacy] E TAMBEM FAZ A EXPORTACAO DE MALHA NA EXTENCAO VTK

> FUNCIONALIDADES
--- ENTRADAS: dados de geometria
--- SAIDA: dados de malha

===============================================================================

> ATUALIZACOES DA VERSAO:
    --- [gen_geo2mesh]: definicao de geometria por meio de lista previa, geo_select --> geo_list:{'frame','shell','solid'}    
    --- [gen_geo2mesh]: Geracao de geometria via terminal interno, sem a necessidade de executar gmsh
    --- [import_mesh_]: datamesh determinado externamente
_______________________________________________________________________________
"""
import sys
import os
import imp
import numpy as np
from colorama import Fore, Back, Style
#-----------------------------------------------------------------------------#

#%%
def create_user_dir(path,file_dir):
    if not os.path.exists(path+file_dir):
        # print('nova pasta criada no diretório "myfempy/user/"')
        os.makedirs(path+file_dir)
    
    # print('esta pasta já existe no diretório  "myfempy/user/"')
    path_user  = str(path+file_dir)
    return path_user 

#%%
def loading_bar_v1(pct,name):
    # sys.stdout.write(Fore.CYAN + Style.BRIGHT+"\r|%-50s"%('#'*int(pct*0.5)+'|'+str(round(pct))+'%'))
    sys.stdout.write(Style.BRIGHT+"\r"+name+": "+"|%-50s"%('#'*int(pct*0.5)+'|'+str(round(pct))+'%'))
    sys.stdout.flush()

#-----------------------------------------------------------------------------#
def animation_mesh(path_user,usr_analysi_name,type_elm,inci,coord,datamesh,Udef,data_result,data_title,typeData,time_out,scale_mesh,ModeNumb):
    from myfempy.setup.myfempy_preproc import create_user_dir
    path_user = create_user_dir(path_user,'/anime_file')
    # scale_mesh_frac0 = np.linspace(-time_out,time_out,time_out)
    scale_mesh_frac0 = np.concatenate((np.linspace(-time_out,time_out,int(time_out/2)),np.linspace(time_out,-time_out,int(time_out/2))),axis=0)
    for itime in range(0,time_out):
        scale_mesh_frac1 = 0.01*scale_mesh*scale_mesh_frac0[itime]
        file_dir=path_user+'/'+usr_analysi_name+'_mode_shape_anime_'+str(ModeNumb)+'_'+(str(itime))+'.vtk'
        meshDefU_anime = np.concatenate((coord[:,0][:, np.newaxis],np.add(coord[:,1:],scale_mesh_frac1*Udef)),axis=1)
        export_mesh(file_dir,type_elm,datamesh[1],datamesh[2],meshDefU_anime,inci,data_result,data_title,typeData)
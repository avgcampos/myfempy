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
#------------------------------------------------------------------------------
import numpy as np
import sys


#------------------------------------------------------------------------------
#%% vetor cond. contorno
def bond_condi(datamesh,boncdnodeaply):
    # numero de tipos de restriçoes fixas
    ntbc = boncdnodeaply.shape[0] 
    
    # dofs fixos do sistema
    fixedof = np.zeros((1,datamesh['fulldof']))
    
    if datamesh['nodedof'] == 2:
        if (datamesh['keyelem'] == "spring20") or (datamesh['keyelem'] == "truss22") or (datamesh['keyelem'] =="plane32") or (datamesh['keyelem'] =="plane42"):
            for ii in range(ntbc):
                no = int(boncdnodeaply[ii,1])
                if int(boncdnodeaply[ii,0]) == 0:
                    fixedof[0,datamesh['nodedof']*no-2] = datamesh['nodedof']*no-1
                    fixedof[0,datamesh['nodedof']*no-1] = datamesh['nodedof']*no
                elif int(boncdnodeaply[ii,0]) == 1:
                    fixedof[0,datamesh['nodedof']*no-2] = datamesh['nodedof']*no-1
                elif int(boncdnodeaply[ii,0]) == 2:
                    fixedof[0,datamesh['nodedof']*no-1] = datamesh['nodedof']*no

        elif datamesh['keyelem'] == "beam21":
            for ii in range(ntbc):
                no = int(boncdnodeaply[ii,1])
                if int(boncdnodeaply[ii,0]) == 0:
                    fixedof[0,datamesh['nodedof']*no-2] = datamesh['nodedof']*no-1
                    fixedof[0,datamesh['nodedof']*no-1] = datamesh['nodedof']*no
                elif int(boncdnodeaply[ii,0]) == 2:
                    fixedof[0,datamesh['nodedof']*no-2] = datamesh['nodedof']*no-1
                elif int(boncdnodeaply[ii,0]) == 6:
                    fixedof[0,datamesh['nodedof']*no-1] = datamesh['nodedof']*no             
    

    elif datamesh['nodedof'] == 3:
        if datamesh['keyelem'] == "frame22":
            for ii in range(ntbc):
                no = int(boncdnodeaply[ii,1])
                if int(boncdnodeaply[ii,0]) == 0:
                    fixedof[0,datamesh['nodedof']*no-3] = datamesh['nodedof']*no-2
                    fixedof[0,datamesh['nodedof']*no-2] = datamesh['nodedof']*no-1
                    fixedof[0,datamesh['nodedof']*no-1] = datamesh['nodedof']*no
                elif int(boncdnodeaply[ii,0]) == 1:
                    fixedof[0,datamesh['nodedof']*no-3] = datamesh['nodedof']*no-2
                elif int(boncdnodeaply[ii,0]) == 2:
                    fixedof[0,datamesh['nodedof']*no-2] = datamesh['nodedof']*no-1
                elif int(boncdnodeaply[ii,0]) == 6:
                    fixedof[0,datamesh['nodedof']*no-1] = datamesh['nodedof']*no  
        
        elif datamesh['keyelem'] == "solid83":
            for ii in range(ntbc):
                no = int(boncdnodeaply[ii,1])
                if int(boncdnodeaply[ii,0]) == 0:
                    fixedof[0,datamesh['nodedof']*no-3] = datamesh['nodedof']*no-2
                    fixedof[0,datamesh['nodedof']*no-2] = datamesh['nodedof']*no-1
                    fixedof[0,datamesh['nodedof']*no-1] = datamesh['nodedof']*no
                elif int(boncdnodeaply[ii,0]) == 1:
                    fixedof[0,datamesh['nodedof']*no-3] = datamesh['nodedof']*no-2
                elif int(boncdnodeaply[ii,0]) == 2:
                    fixedof[0,datamesh['nodedof']*no-2] = datamesh['nodedof']*no-1
                elif int(boncdnodeaply[ii,0]) == 3:
                    fixedof[0,datamesh['nodedof']*no-1] = datamesh['nodedof']*no
    
    elif datamesh['nodedof'] == 6:
        if datamesh['keyelem'] == "frame23":
            for ii in range(ntbc):
                no = int(boncdnodeaply[ii,1])
                if int(boncdnodeaply[ii,0]) == 0:
                    fixedof[0,datamesh['nodedof']*no-6] = datamesh['nodedof']*no-5
                    fixedof[0,datamesh['nodedof']*no-5] = datamesh['nodedof']*no-4
                    fixedof[0,datamesh['nodedof']*no-4] = datamesh['nodedof']*no-3
                    fixedof[0,datamesh['nodedof']*no-3] = datamesh['nodedof']*no-2
                    fixedof[0,datamesh['nodedof']*no-2] = datamesh['nodedof']*no-1
                    fixedof[0,datamesh['nodedof']*no-1] = datamesh['nodedof']*no
                elif int(boncdnodeaply[ii,0]) == 1:
                    fixedof[0,datamesh['nodedof']*no-6] = datamesh['nodedof']*no-5
                elif int(boncdnodeaply[ii,0]) == 2:
                    fixedof[0,datamesh['nodedof']*no-5] = datamesh['nodedof']*no-4
                elif int(boncdnodeaply[ii,0]) == 3:
                    fixedof[0,datamesh['nodedof']*no-4] = datamesh['nodedof']*no-3
                elif int(boncdnodeaply[ii,0]) == 4:
                    fixedof[0,datamesh['nodedof']*no-3] = datamesh['nodedof']*no-2
                elif int(boncdnodeaply[ii,0]) == 5:
                    fixedof[0,datamesh['nodedof']*no-2] = datamesh['nodedof']*no-1
                elif int(boncdnodeaply[ii,0]) == 6:
                    fixedof[0,datamesh['nodedof']*no-1] = datamesh['nodedof']*no
            
    
    fixedof = fixedof[np.nonzero(fixedof)]
    fixedof = fixedof - np.ones_like(fixedof)
    # fixedofs = np.unique(fixedofs)
    # todos os dofs do sistema
    alldof = np.arange(0,datamesh['fulldof'],1,int)
    # dofs livres do sistema
    freedof = np.setdiff1d(alldof,fixedof)
    return freedof, fixedof


#------------------------------------------------------------------------------
#%% vetor forca
def load_apply(datamesh,forcenodeaply):
    forcelist = np.zeros((datamesh['fulldof'],len(np.unique(forcenodeaply[:,3]))))
    for fstep in range(len(np.unique(forcenodeaply[:,3]))):
        forceaply = forcenodeaply[np.where(forcenodeaply[:,3]==fstep+1),:][0]
        nload = forceaply.shape[0]  
        if datamesh['nodedof'] == 2:
            if (datamesh['keyelem'] == "spring20") or (datamesh['keyelem'] == "truss22") or (datamesh['keyelem'] =="plane32") or (datamesh['keyelem'] =="plane42"):
                for ii in range(nload):
                    if int(forceaply[ii,1])==1:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    elif int(forceaply[ii,1])==2:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']-1))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    
            elif datamesh['keyelem'] == "beam21":
                for ii in range(nload):
                    if int(forceaply[ii,1])==2:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    elif int(forceaply[ii,1])==6:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']-1))
                        forcelist[gdlload,fstep] += forceaply[ii,2]              
                    
        elif datamesh['nodedof'] == 3:
            if datamesh['keyelem'] == "frame22":
                for ii in range(nload):
                    if int(forceaply[ii,1])==1:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    elif int(forceaply[ii,1])==2:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']-1))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    elif int(forceaply[ii,1])==6:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']-2))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    
                    
                    
            elif datamesh['keyelem'] == "solid83":
                for ii in range(nload):
                    if int(forceaply[ii,1])==1:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    elif int(forceaply[ii,1])==2:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']-1))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    elif int(forceaply[ii,1])==3:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']-2))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    
    
        elif datamesh['nodedof'] == 6:
            if datamesh['keyelem'] == "frame23":
                for ii in range(nload):
                    if int(forceaply[ii,1])==1:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    elif int(forceaply[ii,1])==2:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']-1))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    elif int(forceaply[ii,1])==3:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']-2))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    elif int(forceaply[ii,1])==4:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']-3))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    elif int(forceaply[ii,1])==5:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']-4))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    elif int(forceaply[ii,1])==6:
                        gdlload = int(datamesh['nodedof']*forceaply[ii,0]-(datamesh['nodedof']-5))
                        forcelist[gdlload,fstep] += forceaply[ii,2]
                    
               
    return forcelist



#------------------------------------------------------------------------------
#%% step setting
def step_setting(usrlog,stepstart,stepend,stepstep):
    
    stepstart = int(usrlog.solver_start)
    stepend = int(usrlog.solver_end)
    stepstep = int(usrlog.solver_step)
    
    if (stepend-stepstart)==0:
        stepset = int(stepend)
    else:
        stepset = int((stepend-stepstart)/stepstep)
        
    return stepset




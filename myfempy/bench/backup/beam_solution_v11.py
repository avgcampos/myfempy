# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 19:34:35 2020
@author: ANTONIO VINICIUS GARCIA CAMPOS
@contact: elementosfinitos.querosaber@gmail.com
@version: BETA V10
@copyright: 3D EASYCAE SERVIÇOS DE ANÁLISE COMPUTACIONAL
@lisence:
    
BSD 3-Clause License

Copyright (c) 2021, Vinicius-Campos
All rights reserved.
_______________________________________________________________________________
~~~~~~ MODULO DE SIMULACAO ESTRUTURAL PELO METODO DOS ELEMENTOS FINITOS ~~~~~~~
~~~~~~                   PROGRAMA DE PROPOSITO GERAL                    ~~~~~~~  
~~~~~~                  copyright all rights reserved                   ~~~~~~~    
===============================================================================

> VERSAO DE LANCAMENTO: V10/2020

===============================================================================
> FUNCIONALIDADES DO PROGRAMA
--- Modulos de analise de estruturas esbeltas com secao constante

--- Integridade com GMSH, pode-se adicionar malha externa, em formato ".msh1"

--- Visualizacao por VTK, pode-se salvar arquivo para visualizar resultados no
paraview ou outros

--- Todo o codigo foi desenvolvido em python3, utilizou-se a biblioteca
Numpy/Scipy para calculo e manipulacao de matrizes, solver do sistema algebrico
utilizando o Numpy, com posibilidade de solver nativo proprio

--- Recomendado utilizar distribuicao Anaconda, com todas as bilbiotecas
atualizadas, necessario instalar via terminal a biblioteca vtk para python3
===============================================================================

> ATUALIZACOES DA VERSAO:
--- Melhoria no Plot da Estrutur, uso do VTK;

--- Geracao de malha internamento no programa, usando gmsh in bat
_______________________________________________________________________________
"""
#%% LOAD DE FUNCOES
import matplotlib.pylab as plt
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg
import sys
import time
from sys import getsizeof
from stf_fini_elm import frame22
from domain_physics import bond_condi_phys, vetor_force, search_nodeXY, search_nodeX, search_nodeY
from prev_process import gen_geo2mesh,import_mesh_line, export_mesh
from post_process import mesh_def_flexural_2d,mesh_def_axial_rigid_2d,intn_force_flexural_full, intn_force_bending_only, react_suport, stress_flexural_full_elm
from view_mesh import view_analysis, view_results, plotter_ShearForce_graphics,plotter_BendingMoment_graphics,plotter_AxialForce_graphics
from my_solver_linear import solver_linear

#%% GERACAO DA GEOMETRIA
geo_list = ['beam','shell','solid']
# Beam 1D Structure
# Shell 2D Structure
# Solid 3D Structure

geo_select = 'beam'

file_geo = 'frame_2d_truss.geo'
name_out_mesh = 'mesh_2d_truss.msh1'
gen_geo2mesh(file_geo,name_out_mesh,geo_select)

# PROPRIEDADES DA SECAO TRANSV. DA ESTRUTURA
secArea = 50*100
secIzz = 4.167E6
secJxx = 1
tabgeo = np.array([[secArea,secIzz]])

#%% GERACAO DA MALHA 
coord, conec, type_elm = import_mesh_line(name_out_mesh)

ngdl = 3        # num gdl por no
nnod = len(coord) # num de nos na malha
nelm = len(conec)  # num de elementos na malha
sgdl = ngdl*nnod  # num total de gdl da malha
datamesh = [ngdl,nnod,nelm,sgdl]

prop_elm = np.ones((datamesh[2],3))
inci = np.concatenate((conec[:,0][:, np.newaxis],prop_elm,conec[:,1:]),axis=1)

#%% PROPRIEDADE DE MATERIAL
modEls = 200E3
cofPoi = 0.3
denVol = 7800
tabmat = np.array([[modEls,cofPoi,denVol]])

#%% FORCAS E CARGAS 

# node_search_force1 = search_nodeXY(0,4000.0,coord,0.002)
# # node_search_force2 = search_nodeXY(8.0,4.0,coord,0.002)
# # node_search_force3 = search_nodeXY(12.0,2.0,coord,0.002)
# frcApy_vet = np.array([[node_search_force1,1,1000]],dtype=object)

node_search_force1 = search_nodeXY(4000.0,2000.0,coord,0.002)
node_search_force2 = search_nodeXY(8000.0,4000.0,coord,0.002)
node_search_force3 = search_nodeXY(12000.0,2000.0,coord,0.002)
frcApy_vet = np.array([[node_search_force1,2,-20000],\
                        [node_search_force2,2,-20000],\
                        [node_search_force3,2,-20000]],dtype=object)

#%% CONDICOES DE CONTORNO
# node_search_bc1 = search_nodeXY(0.0,0.0,coord,0.002)
# node_search_bc2 = search_nodeXY(4000.0,0.0,coord,0.002)
# # list_node_bc = np.array([node_search_bc1[0],node_search_bc2[0]])

# # type_dof : 0 - all dofs fixed; 1 - x diretion fixed; 2 - y diretion fixed
# # list_node_bc : list of nodes fixeds
# # type_dim : 1 - point; 2 - edge; 3 - surface
# # bondCond_vet [type_dof,list_node_bc,type_dim]

# # bondCond_vet = np.array([[0,node_search_bc1],[1,node_search_bc2]],dtype=object)
# bondCond_vet = np.array([[0,node_search_bc1,1],[2,node_search_bc2,1]],dtype=object)
# bondCond_vet = np.array([[0,node_search_bc1,1]],dtype=object)

node_search_bc1 = search_nodeXY(0.0,0.0,coord,0.002)
node_search_bc2 = search_nodeXY(16000.0,0.0,coord,0.002)
# list_node_bc = np.array([node_search_bc1[0],node_search_bc2[0]])

# # type_dof : 0 - all dofs fixed; 1 - x diretion fixed; 2 - y diretion fixed
# # list_node_bc : list of nodes fixeds
# # type_dim : 1 - point; 2 - edge; 3 - surface
# # bondCond_vet [type_dof,list_node_bc,type_dim]

# # bondCond_vet = np.array([[0,node_search_bc1],[1,node_search_bc2]],dtype=object)
bondCond_vet = np.array([[0,node_search_bc1,1],[0,node_search_bc2,1]],dtype=object)


#%% VISUALIZACAO DO PROBLEMA
project_name = 'FRAME ANALISYS'
file_dir='prev_viewmesh_vtk.vtk'

# color = np.random.rand(datamesh[2],1)
color = np.zeros([datamesh[2],1])
export_mesh(file_dir,type_elm,datamesh[1],datamesh[2],coord,conec,color)
# Visualizar a malha VTK
save_screen = 'yes'
file_save='frame_viewmesh_png.png'
scala_view = 800
title_win = 'ANALYSIS VISUALIZATION'
view_analysis(file_dir,file_save,save_screen,datamesh,coord,frcApy_vet,bondCond_vet,project_name,title_win,scala_view)

#%% SOLVER
vtF = vetor_force(datamesh,frcApy_vet)
freedof, fixedof = bond_condi_phys(datamesh,bondCond_vet)

KG = frame22(datamesh,inci,coord,tabmat,tabgeo)
# KG = truss20(datamesh,inci,coord,tabmat,tabgeo)
 
U = np.zeros((datamesh[3],1))
U[freedof] = np.linalg.solve(KG[:,freedof][freedof,:],vtF[freedof])


#%% POS-PROCESSAMENTO
# malha deformada
scale_mesh=600
meshDefUV,UVmeshRotZ = mesh_def_flexural_2d(datamesh,coord,U,scale_mesh)
# meshDefUV = mesh_def_axial_rigid_2d(datamesh,coord,U)

# Tensao normal nos elementos
y_max=+50
strs_elm = stress_flexural_full_elm(datamesh,U,inci,coord,tabmat,tabgeo,y_max)

data_analysis = np.array([[max(abs(U))],[max(strs_elm)],[sum(strs_elm)/datamesh[2]],[min(strs_elm)]])

# exporta malha para visualizacao no VTK
# malha de saida -> vtk ASCII Version 2.0
# visulização do deslocamento
file_dir='mesh_def_view_frame_beam_vtk.vtk'
# export_mesh(file_dir,type_elm,datamesh[1],datamesh[2],meshDefUV,conec,np.ones((datamesh[2],1)))
export_mesh(file_dir,type_elm,datamesh[1],datamesh[2],meshDefUV,conec,strs_elm)
# Visualizar a malha VTK
file_save='mesh_stress_view_frame_png.png'
save_screen = 'yes' 
scale_bar = 'yes'
title_win = 'NODAL DISPLACEMENTS'
title_scale_bar = 'STRESS NORMAL: AXIAL + BENDING' #STRESS: Eqv. Von Mises
view_results(file_dir,file_save,scale_bar,title_scale_bar,save_screen,project_name,title_win,data_analysis)

plotter_AxialForce_graphics(file_save,save_screen,datamesh,inci,coord,Fint[:,0])
# plotter_ShearForce_graphics(file_save,save_screen,datamesh,inci,coord,Fint[:,0])
# plotter_BendingMoment_graphics(file_save,save_screen,datamesh,inci,coord,Fint[:,0])

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
    --- Calculo da tensao média nodal, com base nos valores dos elementos;
    --- Calculo com elementos finitos Q4-iso;
    --- Melhoria no Plot da Estrutura, uso do VTK;
    --- Geracao de malha internamento no programa, usando gmsh in betch
_______________________________________________________________________________
"""
#%% LOAD DE PACOTES
import matplotlib.pylab as plt
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg
import sys
import time
from sys import getsizeof
from stf_fini_elm import plane42_sps
from domain_physics import bond_condi_phys, vetor_force, search_nodeXY, search_nodeX, search_nodeY
from prev_process import gen_geo2mesh, import_mesh_quad , export_mesh
from post_process import mesh_def_plane_2d, stress_vm_plane_q4, stress_avr_plane_q4
from view_mesh import view_analysis, view_results
from my_solver_linear import solver_linear

#%% GERACAO DA GEOMETRIA GMSH

geo_list = ['frame','shell','solid']
# Frame 1D Structure
# Shell 2D Structure
# Solid 3D Structure

geo_select = 'shell'

file_geo = 'plate_q4_mesh.geo'
name_out_mesh = 'file_q4_iso_2d.msh1'
gen_geo2mesh(file_geo,name_out_mesh,geo_select)

# PROPRIEDADES DA SECAO TRANSV. DA ESTRUTURA
secArea = 10 # area da secao, elementos esbeltos
secIzz = 200 # momento de inercia, elementos esbeltos
secJxx = 1   # momento polar de inercia, elementos esbeltos
tck = 10     # espessura, elemento de casca e placa
tabgeo = np.array([[secArea,secIzz,secJxx,tck]])


#%% GERACAO DA MALHA GMSH

coord, conec, type_elm = import_mesh_quad(name_out_mesh)

ngdl = 2           # num gdl por no
nnod = len(coord)  # num de nos na malha
nelm = len(conec)  # num de elementos na malha
sgdl = ngdl*nnod   # num total de gdl da malha
datamesh = [ngdl,nnod,nelm,sgdl]

prop_elm = np.ones((datamesh[2],3))
inci = np.concatenate((conec[:,0][:, np.newaxis],prop_elm,conec[:,1:]),axis=1)

#%% PROPRIEDADE DE MATERIAL
typeMechMat = "isotropic_planeStress_2d"
modEls = 200000 # modulo de elasticidade
cofPoi = 0.3 # coeficiente de poisson
denVol = 7800 # densidade do material
tabmat = np.array([[modEls,cofPoi,denVol]])

#%% FORCAS E CARGAS 
node_search_force1 = search_nodeXY(20,10,coord,0.002)
# node_search_force2 = search_nodeXY(100,120,coord,0.002)
# frcApy_vet = np.array([[node_search_force1,1,10000],[node_search_force2,3,5000]],dtype=object)
# frcApy_vet = np.array([[node_search_force1,2,-10000],[node_search_force2,2,-10000]],dtype=object)
frcApy_vet = np.array([[node_search_force1,2,-500]],dtype=object)

#%% CONDICOES DE CONTORNO

# coordx=(coord[np.where(coord[:,1]<=12.0),:])[0]

node_search_bc1 = search_nodeX(0.0,coord,0.002)
# node_search_bc2 = search_nodeXY(120.0,0.0,coord,0.002)

# type_dof : 0 - all dofs fixed; 1 - x diretion fixed; 2 - y diretion fixed
# list_node_bc : list of nodes fixeds
# type_dim : 1 - point; 2 - edge; 3 - surface
# bondCond_vet [type_dof,list_node_bc,type_dim]

# bondCond_vet = np.array([[0,node_search_bc1],[0,node_search_bc2]],dtype=object)
# bondCond_vet = np.array([[0,node_search_bc1,node_search_bc2]],dtype=object)
bondCond_vet = np.array([[0,node_search_bc1,2]],dtype=object)

#%% VISUALIZACAO DO PROBLEMA
project_name = 'Analise de Braco para porca'
file_dir='prev_viewmesh_vtk.vtk'
file_save='shell13_viewmesh_png.png'
# color = np.random.rand(datamesh[2],1)
color = np.zeros([datamesh[2],1])
typeData='elm'
export_mesh(file_dir,type_elm,datamesh[1],datamesh[2],coord,conec,color,typeData)
# Visualizar a malha VTK
save_screen = 'no'
scala_view = 2
title_win = 'ANALYSIS VISUALIZATION'
view_analysis(file_dir,file_save,save_screen,datamesh,coord,frcApy_vet,bondCond_vet,project_name,title_win,scala_view)

#%% SOLVER
vtF = vetor_force(datamesh,frcApy_vet)
freedof, fixedof = bond_condi_phys(datamesh,bondCond_vet)

KG,A = plane42_sps(datamesh,conec,coord,typeMechMat,tabmat,tck,2)

U = np.zeros((datamesh[3],1))
U[freedof,0] = sp.linalg.spsolve(KG[:,freedof][freedof,:],vtF[freedof])

#%% POS-PROCESSAMENTO
# malha deformada
scale_mesh=12
meshDef = mesh_def_plane_2d(datamesh,coord,U,scale_mesh)

# # # # tensao nos elementos von-mises
strs_elm_eqv,strs_elm_xx,strs_elm_yy,strs_elm_xy = stress_vm_plane_q4(U,datamesh,conec,coord,typeMechMat,tabmat,tck,2)

strs_avr_eqv, S = stress_avr_plane_q4(strs_elm_eqv,datamesh,inci,conec)

print('deslo max: ',max(abs(U)),' [mm]')
print('tensao EQV media: ',sum(strs_elm_eqv)/datamesh[2],' [MPa]')
print('tensao EQV max: ',max(strs_elm_eqv),' [MPa]')
print('tensao EQV min: ',min(strs_elm_eqv),' [MPa]')
print('----------------------------------------')
print("tensao XX",strs_elm_xx,'\n')
print("tensao YY",strs_elm_yy,'\n')
print("tensao XY",strs_elm_xy,'\n')

data_analysis1 = np.array([[max(abs(U))],[max(strs_avr_eqv)],[sum(strs_avr_eqv)/datamesh[1]],[min(strs_avr_eqv)]])

# exporta malha para visualizacao no VTK
# malha de saida -> vtk ASCII Version 2.0
# visulização do deslocamento
file_dir='meshview_shell13_vtk.vtk'
typeData='avr'
export_mesh(file_dir,type_elm,datamesh[1],datamesh[2],meshDef,conec,strs_avr_eqv,typeData)
# Visualizar a malha VTK
file_save='meshview_shell13_Sts_avr_png.png'
title_win = 'NODAL DISPLACEMENTS'
save_screen = 'yes'
scale_bar = 'yes'
title_scale_bar = 'STRESS EQV-VM NODAL AVR'
view_results(file_dir,file_save,scale_bar,title_scale_bar,save_screen,project_name,title_win,data_analysis1)

data_analysis2 = np.array([[max(abs(U))],[max(strs_elm_eqv)],[sum(strs_elm_eqv)/datamesh[2]],[min(strs_elm_eqv)]])

# exporta malha para visualizacao no VTK
# malha de saida -> vtk ASCII Version 2.0
# visulização do deslocamento
file_dir='meshview_shell132_vtk.vtk'
typeData='elm'
export_mesh(file_dir,type_elm,datamesh[1],datamesh[2],meshDef,conec,strs_elm_eqv,typeData)
# Visualizar a malha VTK
file_save='meshview_shell13_Sts_elm_png.png'
title_win = 'NODAL DISPLACEMENTS'
save_screen = 'yes'
scale_bar = 'yes'
title_scale_bar = 'STRESS EQV-VM ELM MEAN'
view_results(file_dir,file_save,scale_bar,title_scale_bar,save_screen,project_name,title_win,data_analysis2)

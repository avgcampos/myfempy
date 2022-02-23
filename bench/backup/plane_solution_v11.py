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
    --- Melhoria no Plot da Estrutura, uso do VTK;
    --- Geracao de malha internamento no programa, usando gmsh in bat
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
from stf_fini_elm import plane32_sps
from domain_physics import bond_condi_phys, vetor_force, search_nodeXY, search_nodeX, search_nodeY
from prev_process import gen_geo2mesh, import_mesh_tria , export_mesh
from post_process import mesh_def_plane_2d, stress_vm_plane
from view_mesh import view_analysis, view_results
from my_solver_linear import solver_linear

#%% GERACAO DA GEOMETRIA GMSH

geo_list = ['beam','plane','solid']
# Frame 1D Structure
# Shell 2D Structure
# Solid 3D Structure

geo_select = 'plane'

file_geo = 'shell_2d_plate.geo'
name_out_mesh = 'mesh_2d_plate.msh1'
gen_geo2mesh(file_geo,name_out_mesh,geo_select)

# PROPRIEDADES DA SECAO TRANSV. DA ESTRUTURA
secArea = 10 # area da secao, elementos esbeltos
secIzz = 200 # momento de inercia, elementos esbeltos
secJxx = 1   # momento polar de inercia, elementos esbeltos
tck = 20      # espessura, elemento de casca e placa
tabgeo = np.array([[secArea,secIzz,secJxx,tck]])


#%% GERACAO DA MALHA GMSH

coord, conec, type_elm = import_mesh_tria(name_out_mesh)

ngdl = 2           # num gdl por no
nnod = len(coord)  # num de nos na malha
nelm = len(conec)  # num de elementos na malha
sgdl = ngdl*nnod   # num total de gdl da malha
datamesh = [ngdl,nnod,nelm,sgdl]

prop_elm = np.ones((datamesh[2],3))
inci = np.concatenate((conec[:,0][:, np.newaxis],prop_elm,conec[:,1:]),axis=1)

#%% PROPRIEDADE DE MATERIAL
modEls = 200000 # modulo de elasticidade
cofPoi = 0.3 # coeficiente de poisson
denVol = 7800 # densidade do material
tabmat = np.array([[modEls,cofPoi,denVol]])

#%% FORCAS E CARGAS 
node_search_force1 = search_nodeXY(0,120.0,coord,0.002)
node_search_force2 = search_nodeXY(100,120,coord,0.002)
# frcApy_vet = np.array([[node_search_force1,1,10000],[node_search_force2,3,5000]],dtype=object)
# frcApy_vet = np.array([[node_search_force1,2,-10000],[node_search_force2,2,-10000]],dtype=object)
frcApy_vet = np.array([[node_search_force1,2,-10000]],dtype=object)

#%% CONDICOES DE CONTORNO
node_search_bc1 = search_nodeX(0.0,coord,0.002)
node_search_bc2 = search_nodeXY(200.0,0.0,coord,0.002)
# node_search_bc2 = search_nodeXY(120.0,0.0,coord,0.002)

# type_dof : 0 - all dofs fixed; 1 - x diretion fixed; 2 - y diretion fixed
# list_node_bc : list of nodes fixeds
# type_dim : 1 - point; 2 - edge; 3 - surface
# bondCond_vet [type_dof,list_node_bc,type_dim]

# bondCond_vet = np.array([[0,node_search_bc1],[0,node_search_bc2]],dtype=object)
# bondCond_vet = np.array([[0,node_search_bc1,node_search_bc2]],dtype=object)
bondCond_vet = np.array([[0,node_search_bc2*np.ones_like(node_search_bc1),1],[1,node_search_bc1,2]],dtype=object)

#%% VISUALIZACAO DO PROBLEMA
project_name = 'Shell Analysis 2021'
file_dir='prev_viewmesh_vtk.vtk'
file_save='shell11_viewmesh_png.png'
# color = np.random.rand(datamesh[2],1)
color = np.zeros([datamesh[2],1])
export_mesh(file_dir,type_elm,datamesh[1],datamesh[2],coord,conec,color)
# Visualizar a malha VTK
save_screen = 'yes'
scala_view = 20
title_win = 'ANALYSIS VISUALIZATION'
view_analysis(file_dir,file_save,save_screen,datamesh,coord,frcApy_vet,bondCond_vet,project_name,title_win,scala_view)


#%% SOLVER
vtF = vetor_force(datamesh,frcApy_vet)
freedof, fixedof = bond_condi_phys(datamesh,bondCond_vet)

KG,A = plane32_sps(datamesh,conec,coord,tabmat,tabgeo)

U = np.zeros((datamesh[3],1))
U[freedof,0] = sp.linalg.spsolve(KG[:,freedof][freedof,:],vtF[freedof])

#%% POS-PROCESSAMENTO
# malha deformada
scale_mesh=10000
meshDef = mesh_def_plane_2d(datamesh,coord,U,scale_mesh)

# # # # tensao nos elementos von-mises
stress_eqv_VM = stress_vm_plane(U,conec,coord,tabmat,datamesh)

print('deslo max: ',max(abs(U)),' [mm]')
print('tensao media: ',sum(stress_eqv_VM)/datamesh[2],' [MPa]')
print('tensao max: ',max(stress_eqv_VM),' [MPa]')
print('tensao min: ',min(stress_eqv_VM),' [MPa]')

data_analysis = np.array([[max(abs(U))],[max(stress_eqv_VM)],[sum(stress_eqv_VM)/datamesh[2]],[min(stress_eqv_VM)]])

# exporta malha para visualizacao no VTK
# malha de saida -> vtk ASCII Version 2.0
# visulização do deslocamento
file_dir='meshview_shell11_vtk.vtk'
export_mesh(file_dir,type_elm,datamesh[1],datamesh[2],meshDef,conec,np.sqrt(stress_eqv_VM))
# Visualizar a malha VTK
file_save='meshview_stress_shell11_1mm.png'
title_win = 'Stress_1mm'
save_screen = 'yes'
scale_bar = 'yes'
title_scale_bar = 'STRESS EQV-VM NORM SQRT'
view_results(file_dir,file_save,scale_bar,title_scale_bar,save_screen,project_name,title_win,data_analysis)

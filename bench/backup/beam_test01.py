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

from analysiconfig import *


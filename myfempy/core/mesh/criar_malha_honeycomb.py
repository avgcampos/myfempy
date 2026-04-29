# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 15:29:23 2026

@author: antvi
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def gerar_malha_honeycomb(nx, ny, a_geom, t_geom, mostrar_numeros):
    # 1. Definição das Dimensões Geométricas Fixas    
    # Proporções da célula para mapeamento
    # Largura total teórica = 3 * a
    # Altura total teórica = a * sqrt(3)
    L_total = 3 * a_geom
    H_total = a_geom * np.sqrt(3)
    
    # Matriz de densidade
    malha = np.zeros((ny, nx), dtype=int)
    
    mid_y = H_total / 2
    quarto_a = a_geom / 2

    # 2. Preenchimento da Estrutura baseada em coordenadas normalizadas
    for y_idx in range(ny):
        # Coordenada Y do centro do elemento
        y_coord = (y_idx + 0.5) * (H_total / ny)
        dy = abs(y_coord - mid_y)
        
        for x_idx in range(nx):
            # Coordenada X do centro do elemento
            x_coord = (x_idx + 0.5) * (L_total / nx)
            
            e_solido = False
            
            # A. Conectores laterais (Extremidades esquerda e direita)
            if dy <= t_geom: #* 0.8:
                if x_coord < quarto_a or x_coord >= (L_total - quarto_a):
                    e_solido = True
            
            # B. Faces Horizontais (Topo e Base central)
            if (a_geom <= x_coord < 2 * a_geom):
                if y_coord < t_geom or y_coord >= (H_total - t_geom):
                    e_solido = True
            
            # C. Paredes Inclinadas (Diagonais)
            # Relação linear para a inclinação
            x_ref_l = quarto_a + (dy * (a_geom - quarto_a) / mid_y)
            
            # Margem de tolerância baseada na espessura t
            # Usamos t_geom / cos(theta) aproximadamente para manter a espessura constante
            tolerancia = t_geom #* 0.8 
            
            # Metade esquerda
            if quarto_a <= x_coord < a_geom:
                if abs(x_coord - x_ref_l) < tolerancia:
                    e_solido = True
            
            # Metade direita (espelhada)
            x_ref_r = L_total - x_ref_l
            if (L_total - a_geom) <= x_coord < (L_total - quarto_a):
                if abs(x_coord - x_ref_r) < tolerancia:
                    e_solido = True

            if e_solido:
                malha[y_idx, x_idx] = 1

    # 3. Plotagem
    fig, ax = plt.subplots(figsize=(12, 8))
    cor_solido = '#2c3e50'
    cor_vazio = '#ffffff'
    
    for y in range(ny):
        for x in range(nx):
            indice = y * nx + x  
            cor = cor_solido if malha[y, x] == 1 else cor_vazio
            
            rect = Rectangle((x, y), 1, 1, facecolor=cor, edgecolor='black', linewidth=0.2)
            ax.add_patch(rect)
            
            if mostrar_numeros:
                txt_col = 'white' if malha[y, x] == 1 else 'black'
                ax.text(x + 0.5, y + 0.5, str(indice), va='center', ha='center', 
                        fontsize=7, color=txt_col, alpha=0.8)

    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    ax.set_aspect('equal')
    
    plt.title(f"Malha Honeycomb Discretizada: {nx}x{ny}\nDimensões: a={a_geom:.4f}, t={t_geom:.4f}", fontsize=12)
    plt.axis('off')
    plt.tight_layout()
    plt.show()

    indices_solidos = np.where(malha.flatten(order='C') == 1)[0]
    return indices_solidos

# # --- Execução ---
# # Defina aqui o número de elementos desejado para a malha
# # a = 3^(-3/4)
# a = 3**(-0.75) 
# # t/a = sqrt(3)/12 -> t = a * sqrt(3)/12
# t = a * (np.sqrt(3) / 12)
# ids_estruturais = gerar_malha_honeycomb_v2(nx=24, ny=12, a_geom = a, t_geom = t)

# print(f"Total de elementos sólidos: {len(ids_estruturais)}")
# print(f"Índices estruturais:")
# print(ids_estruturais)
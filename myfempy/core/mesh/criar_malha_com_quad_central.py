# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 16:11:00 2026

@author: Vinicius
"""

import matplotlib.pyplot as plt

def criar_malha_com_quadrado(nx, ny, tamanho_central=4):
    """
    Cria uma malha de elementos finitos quadrados, numera cada elemento
    a partir de 1 e retorna lista dos elementos dentro do quadrado central.
    
    nx, ny: dimensões da malha
    tamanho_central: tamanho do quadrado central (em número de elementos)
    """
    elementos_centro = []
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Calcular limites do quadrado central
    inicio_x = (nx - tamanho_central) // 2
    fim_x = inicio_x + tamanho_central - 1
    inicio_y = (ny - tamanho_central) // 2
    fim_y = inicio_y + tamanho_central - 1
    
    num = 0  # agora começa em 0
    for j in range(ny):
        for i in range(nx):
            # Cor padrão
            cor = "lightgray"
            # Se está dentro do quadrado central, muda cor
            if inicio_x <= i <= fim_x and inicio_y <= j <= fim_y:
                cor = "salmon"
                elementos_centro.append(num)
            
            # Desenha quadrado
            rect = plt.Rectangle((i, j), 1, 1, facecolor=cor, edgecolor="black")
            ax.add_patch(rect)
            
            # # Escreve número no centro
            # ax.text(i + 0.5, j + 0.5, str(num),
            #         ha="center", va="center", fontsize=6, color="black")
            
            num += 1
    
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    ax.set_aspect("equal")
    plt.title(f"Malha {nx}x{ny} com quadrado central {tamanho_central}x{tamanho_central}")
    plt.show()
    
    return elementos_centro

# # # Exemplo: malha 20x20 com quadrado central 4x4
# lista_centro = criar_malha_com_quadrado(10, 10, tamanho_central=4)

# print("Elementos dentro do quadrado central:", lista_centro)

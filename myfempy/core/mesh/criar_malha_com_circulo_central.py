# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 16:15:00 2026

@author: Vinicius
"""

import matplotlib.pyplot as plt
import math

def criar_malha_com_circulo(nx, ny, raio=4):
    """
    Cria uma malha de elementos finitos quadrados, numera cada elemento
    a partir de 1 e retorna lista dos elementos dentro de um círculo central.
    
    nx, ny: dimensões da malha
    raio: raio do círculo central (em número de elementos)
    """
    elementos_centro = []
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Centro da malha
    cx, cy = nx / 2, ny / 2
    
    num = 1  # numeração começa em 1
    for j in range(ny):
        for i in range(nx):
            # Coordenadas do centro do quadrado
            x_centro = i + 0.5
            y_centro = j + 0.5
            
            # Distância até o centro da malha
            dist = math.sqrt((x_centro - cx)**2 + (y_centro - cy)**2)
            
            # Cor padrão
            cor = "lightgray"
            # Se está dentro do círculo, muda cor
            if dist <= raio:
                cor = "salmon"
                elementos_centro.append(num)
            
            # Desenha círculo
            rect = plt.Rectangle((i, j), 1, 1, facecolor=cor, edgecolor="black")
            ax.add_patch(rect)
            
            # # Escreve número no centro
            ax.text(i + 0.5, j + 0.5, str(num),
                    ha="center", va="center", fontsize=6, color="black")
            
            num += 1
    
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    ax.set_aspect("equal")
    plt.title(f"Malha {nx}x{ny} com círculo central (raio={raio})")
    plt.show()
    
    return elementos_centro

# Exemplo: malha 20x20 com círculo central de raio 4
lista_centro = criar_malha_com_circulo(10, 10, raio=4)

print("Elementos dentro do círculo central:", lista_centro)

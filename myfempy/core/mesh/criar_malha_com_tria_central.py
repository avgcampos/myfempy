# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 16:20:00 2026

@author: Vinicius
"""

import matplotlib.pyplot as plt
import math

def ponto_no_triangulo(px, py, v1, v2, v3):
    """
    Verifica se o ponto (px, py) está dentro do triângulo definido por v1, v2, v3.
    Usa coordenadas baricêntricas.
    """
    denom = ((v2[1] - v3[1]) * (v1[0] - v3[0]) +
             (v3[0] - v2[0]) * (v1[1] - v3[1]))
    a = ((v2[1] - v3[1]) * (px - v3[0]) +
         (v3[0] - v2[0]) * (py - v3[1])) / denom
    b = ((v3[1] - v1[1]) * (px - v3[0]) +
         (v1[0] - v3[0]) * (py - v3[1])) / denom
    c = 1 - a - b
    return (0 <= a <= 1) and (0 <= b <= 1) and (0 <= c <= 1)

def criar_malha_com_triangulo(nx, ny, lado=8):
    """
    Cria uma malha de elementos finitos quadrados, numera cada elemento
    a partir de 1 e retorna lista dos elementos dentro de um triângulo equilátero central.
    
    nx, ny: dimensões da malha
    lado: tamanho do lado do triângulo (em número de elementos)
    """
    elementos_centro = []
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Centro da malha
    cx, cy = nx / 2, ny / 2
    
    # Altura do triângulo equilátero
    h = (math.sqrt(3) / 2) * lado
    
    # Vértices do triângulo central (apontando para cima)
    v1 = (cx, cy + h/2)          # vértice superior
    v2 = (cx - lado/2, cy - h/2) # vértice inferior esquerdo
    v3 = (cx + lado/2, cy - h/2) # vértice inferior direito
    
    num = 1
    for j in range(ny):
        for i in range(nx):
            # Centro do quadrado
            x_centro = i + 0.5
            y_centro = j + 0.5
            
            # Cor padrão
            cor = "lightgray"
            # Se está dentro do triângulo, muda cor
            if ponto_no_triangulo(x_centro, y_centro, v1, v2, v3):
                cor = "salmon"
                elementos_centro.append(num)
            
            # Desenha quadrado
            rect = plt.Rectangle((i, j), 1, 1, facecolor=cor, edgecolor="black")
            ax.add_patch(rect)
            
            # Escreve número no centro
            ax.text(i + 0.5, j + 0.5, str(num),
                    ha="center", va="center", fontsize=6, color="black")
            
            num += 1
    
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)
    ax.set_aspect("equal")
    plt.title(f"Malha {nx}x{ny} com triângulo equilátero central (lado={lado})")
    plt.show()
    
    return elementos_centro

# Exemplo: malha 20x20 com triângulo equilátero central de lado 8
lista_centro = criar_malha_com_triangulo(21, 21, lado=10)

print("Elementos dentro do triângulo central:", lista_centro)


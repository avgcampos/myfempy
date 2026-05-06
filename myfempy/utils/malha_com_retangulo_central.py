import matplotlib.pyplot as plt

def criar_malha_com_retangulo(nx, ny, a, b):
    """
    Cria uma malha de elementos finitos, numera cada elemento
    e destaca um retângulo central de dimensões a por b.
    
    nx, ny: dimensões totais da malha (nº de elementos)
    a: largura do retângulo central (dimensão x)
    b: altura do retângulo central (dimensão y)
    """
    elementos_centro = []
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Calcular limites do retângulo central para x (a)
    inicio_x = (nx - a) // 2
    fim_x = inicio_x + a - 1
    
    # Calcular limites do retângulo central para y (b)
    inicio_y = (ny - b) // 2
    fim_y = inicio_y + b - 1
    
    num = 0 
    for j in range(ny):
        for i in range(nx):
            # Cor padrão (fundo)
            cor = "lightgray"
            
            # Verifica se a célula atual (i, j) está dentro dos limites do retângulo
            if inicio_x <= i <= fim_x and inicio_y <= j <= fim_y:
                cor = "salmon"
                elementos_centro.append(num)
            
            # Desenha o elemento
            rect = plt.Rectangle((i, j), 1, 1, facecolor=cor, edgecolor="black")
            ax.add_patch(rect)
            
            num += 1
    
    # ax.set_xlim(0, nx)
    # ax.set_ylim(0, ny)
    # ax.set_aspect("equal")
    # plt.title(f"Malha {nx}x{ny} com retângulo central {a}x{b}")
    # plt.grid(False)
    # plt.show()
    
    return elementos_centro

# Exemplo de uso: Malha 20x15 com um retângulo central 8x4
# lista_centro = criar_malha_com_retangulo(nx=10, ny=10, a=4, b=6)

# print(f"Total de elementos no retângulo: {len(lista_centro)}")
# print("IDs dos elementos:", lista_centro)
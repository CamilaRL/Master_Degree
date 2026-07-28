import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
    

def Distribution(beta_R, dmn):

    n = 1/(np.exp(beta_R * dmn) - 1)
    
    f = n + 0.5
        
    return f
    

### MAIN

w = 1

beta_hot = 0.2085149450066301
beta_cold = 2.264392803150969

mu = 2.0

fh = Distribution(beta_hot, w)
fc = Distribution(beta_cold, w)

fi = [fc, fh]
cores = ['navy', 'red']
cores_meio = ['skyblue', 'pink']
nome = ['Cold', 'Hot']

fig, ax = plt.subplots(figsize=(5, 5))

for i in range(2):

    r = np.array([mu, mu])

    V = np.array([
        [fi[i], 0.0],
        [0.0, fi[i]]
    ])

    # Autovalores e autovetores de V
    values, vectors = np.linalg.eigh(V)

    # Ângulo de rotação da elipse em graus
    angle = np.degrees(np.arctan2(*vectors[:, 0][::-1]))

    width, height = 2 * np.sqrt(values)

    # Desenha a elipse e o ponto médio
    ellipse = Ellipse(xy=r, width=width, height=height, angle=angle, edgecolor=cores[i], fc=cores_meio[i], alpha=0.5, lw=2, label=nome[i])
    ax.add_patch(ellipse)
    ax.plot(r[0], r[1], 'o', color=cores[i])

# --- 4. Centraliza os Eixos no (0,0) ---
# Esconde as bordas superior e direita
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')

# Move os eixos esquerdo e inferior para a posição (0, 0)
ax.spines['left'].set_position(('data', 0))
ax.spines['bottom'].set_position(('data', 0))

# Adiciona setas nas pontas dos eixos (opcional, dá um visual mais limpo)
ax.plot(1, 0, ">k", transform=ax.get_yaxis_transform(), clip_on=False)
ax.plot(0, 1, "^k", transform=ax.get_xaxis_transform(), clip_on=False)

ax.set_xticks([])
ax.set_yticks([])

# --- 5. Ajustes Finais ---
# Define limites simétricos ao redor do zero para visualização centralizada
#max_val = max(abs(r[0]) + width, abs(r[1]) + height)
max_val = max(width, height)
ax.set_xlim(-max_val, max_val)
ax.set_ylim(-max_val, max_val)

# Rótulos dos eixos (posicionados nas pontas)
ax.set_xlabel('$q$', fontsize=20, loc='right')
ax.set_ylabel('$p$', fontsize=20, loc='top', rotation=0)

#ax.grid(True, linestyle='--', alpha=0.5)
ax.set_aspect('equal')  # Garante escala 1:1 para não distorcer a elipse
#plt.legend(loc='upper left')

plt.show()
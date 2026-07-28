import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
    

def Distribution(beta_R, dmn):

    n = 1/(np.exp(beta_R * dmn) - 1)
    
    f = n + 0.5
        
    return f
    

### MAIN

w = 1


rh, beta_hot, rc, beta_cold = np.loadtxt('./ThermalKinematics/initial_temperatures.txt', unpack=True)

cores = ['navy', 'red']
cores_meio = ['skyblue', 'pink']

for j in range(2):

    fh = Distribution(beta_hot[j], w)
    fc = Distribution(beta_cold[j], w)

    ri = [rc[j], rh[j]]
    fi = [fc, fh]

    
    fig, ax = plt.subplots(figsize=(5, 5))

    for i in range(2):

        nome = [r'$r_c$ = '+f'{ri[0]:.2f}', r'$r_h$ = '+f'{ri[1]:.2f}']
        
        r = np.array([0.0, 0.0])

        V = np.array([
            [fi[i]*np.cosh(2*ri[i]), -fi[i]*2*np.sinh(ri[i])*np.cosh(ri[i])],
            [-fi[i]*2*np.sinh(ri[i])*np.cosh(ri[i]), fi[i]*np.cosh(2*ri[i])]
        ])

        # Autovalores e autovetores de V
        values, vectors = np.linalg.eigh(V)

        # Ângulo de rotação da elipse em graus
        #angle = [0.0, 0.0] ## np.degrees(np.arctan2(*vectors[:, 0][::-1]))

        width, height = 2 * np.sqrt(values)

        # Desenha a elipse e o ponto médio
        ellipse = Ellipse(xy=r, width=width, height=height, edgecolor=cores[i], fc=cores_meio[i], alpha=0.5, lw=2)
        ax.add_patch(ellipse)
        ax.plot(r[0], r[1], 'o', color=cores[i], label=nome[i])

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
    plt.legend(loc='upper left', fontsize=15)
    
    plt.show()
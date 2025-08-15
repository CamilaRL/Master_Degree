import numpy as np
import matplotlib.pyplot as plt


tot = 7
cmap = plt.get_cmap('rainbow')
colors = iter(cmap(np.linspace(0.01, 1, tot)))

curvas = []
completion1 = []
completion2 = []

for i in range(1, tot, 1):

    tlist, completion = np.loadtxt(f'./ThermalKinematics_Aquecer/completion_{i}.txt', unpack=True)
    
    volta1 = np.where(0.71 == tlist)
    volta2 = np.where(2.2 == tlist)
    
    curvas.append(i)
    completion1.append(completion[volta1])
    completion2.append(completion[volta2])
    
    cor = next(colors)
    
    plt.plot(tlist, completion, color=cor, label=f'Curve {i}')
    plt.xlabel('Time')
    plt.ylabel('Degree of Completion')

plt.legend(loc='best', bbox_to_anchor=(1., 0.5, 0.5, 0.5))
plt.tight_layout()
plt.show()


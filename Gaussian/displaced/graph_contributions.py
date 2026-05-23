import numpy as np
import matplotlib.pyplot as plt


muList = [0.1, 0.5, 1.0, 2.0]
symbols = ['-', '--', ':']
labels = ['Total', 'Passive Contribution', 'Ergotropic Contribution']

## Leitura Arquivos

Iw_p_list = []
Iw_e_list = []
Sprod_p_list = []
Sprod_e_list = []

Iw_list = []
Sprod_list = []


for k, mu in enumerate(muList):
    
    tlist, Iw_p, Iw_e, Sprod_p, Sprod_e = np.loadtxt(f'./ThermalKinematics_Contributions/mu{mu}-heating-contributions.txt', unpack=True)
    
    
    Iw_p_list.append([Iw_p])
    Iw_e_list.append([Iw_e])
    Sprod_p_list.append([Sprod_p])
    Sprod_e_list.append([Sprod_e])

    tlist, Iw_p, Iw_e, Sprod_p, Sprod_e = np.loadtxt(f'./ThermalKinematics_Contributions/mu{mu}-cooling-contributions.txt', unpack=True)
    
    Iw_p_list[k].append(Iw_p)
    Iw_e_list[k].append(Iw_e)
    Sprod_p_list[k].append(Sprod_p)
    Sprod_e_list[k].append(Sprod_e)
    
    
    Iw, Sprod = np.loadtxt(f'./ThermalKinematics/mu{mu}-heating.txt', unpack=True, usecols=(1,6))
    
    Iw_list.append([Iw])
    Sprod_list.append([Sprod])
    
    Iw, Sprod = np.loadtxt(f'./ThermalKinematics/mu{mu}-cooling.txt', unpack=True, usecols=(1,6))
    
    Iw_list[k].append(Iw)
    Sprod_list[k].append(Sprod)
    

## Plots
'''
### Wigner Fisher Information

##### Heating

plt.figure(figsize=(15,5))

for i in range(len(muList)):
    
    num = 141 + i
    
    plt.subplot(num)
    plt.plot(tlist, Iw_list[i][0], color='red', linestyle=symbols[0], linewidth=2)
    plt.plot(tlist, Iw_p_list[i][0], color='red', linestyle=symbols[1], linewidth=2)
    plt.plot(tlist, Iw_e_list[i][0], color='red', linestyle=symbols[2], linewidth=2)
    
    plt.xscale('log')
    plt.title(r'$\mu$ = '+f'{muList[i]}', fontsize=12)
    plt.xlabel('Time', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    if i == 0:
        plt.ylabel('Wigner Fisher Information', fontsize=12)
    if i == len(muList)-1:
        plt.legend(fontsize=12, loc='center left', bbox_to_anchor=(1, 0.5))
        
    
    
plt.suptitle('Heating', fontsize=14)
plt.tight_layout()
plt.subplots_adjust(right=0.8)
plt.show()

##### Cooling

plt.figure(figsize=(15,5))

for i in range(len(muList)):
    
    num = 141 + i
    
    plt.subplot(num)
    plt.plot(tlist, Iw_list[i][1], color='blue', linestyle=symbols[0], linewidth=2)
    plt.plot(tlist, Iw_p_list[i][1], color='blue', linestyle=symbols[1], linewidth=2)
    plt.plot(tlist, Iw_e_list[i][1], color='blue', linestyle=symbols[2], linewidth=2)
    
    plt.xscale('log')
    plt.title(r'$\mu$ = '+f'{muList[i]}', fontsize=12)
    plt.xlabel('Time', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    if i == 0:
        plt.ylabel('Wigner Fisher Information', fontsize=12)
    if i == len(muList)-1:
        plt.legend(fontsize=12, loc='center left', bbox_to_anchor=(1, 0.5))
    
plt.suptitle('Cooling', fontsize=14)
plt.tight_layout()
plt.subplots_adjust(right=0.8)
plt.show()


### Entropy Production


plt.figure(figsize=(15,10))

for i in range(len(muList)):
    
    ##### Heating
    num_h = 241 + i
    
    plt.subplot(num)
    plt.plot(tlist, Sprod_list[i][0], color='red', linestyle=symbols[0], linewidth=2, label='Total')
    plt.plot(tlist, Sprod_p_list[i][0], color='red', linestyle=symbols[1], linewidth=2, label='Passive Contribution')
    plt.plot(tlist, Sprod_e_list[i][0], color='red', linestyle=symbols[2], linewidth=2, label='Ergotropic Contribution')
    
    plt.xscale('log')
    plt.title(r'$\mu$ = '+f'{muList[i]}', fontsize=12)
    plt.xlabel('Time', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    if i == 0:
        plt.ylabel('Entropy Production Rate', fontsize=12)
    
    ##### Cooling
    num_c = 244 + i
    
    plt.subplot(num)
    plt.plot(tlist, Sprod_list[i][1], color='blue', linestyle=symbols[0], linewidth=2, label='Total')
    plt.plot(tlist, Sprod_p_list[i][1], color='blue', linestyle=symbols[1], linewidth=2, label='Passive Contribution')
    plt.plot(tlist, Sprod_e_list[i][1], color='blue', linestyle=symbols[2], linewidth=2, label='Ergotropic Contribution')
    
    plt.xscale('log')
    plt.title(r'$\mu$ = '+f'{muList[i]}', fontsize=12)
    plt.xlabel('Time', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    if i == 0:
        plt.ylabel('Entropy Production Rate', fontsize=12)
    
#plt.suptitle('Cooling', fontsize=14)
#plt.suptitle('Heating', fontsize=14)
plt.tight_layout()
plt.subplots_adjust(right=0.8)
plt.show()'''

fig = plt.figure(figsize=(12, 8)) # Aumentei um pouco a altura para caber a legenda

for i in range(len(muList)):
    
    ##### HEATING
    ax_h = plt.subplot(2, 4, i + 1) # 2 linhas, 4 colunas, posição 5 a 8
    
    line_h1 = plt.plot(tlist, Sprod_list[i][0], color='red', linestyle=symbols[0], linewidth=2)
    line_h2 = plt.plot(tlist, Sprod_p_list[i][0], color='red', linestyle=symbols[1], linewidth=2)
    line_h3 = plt.plot(tlist, Sprod_e_list[i][0], color='red', linestyle=symbols[2], linewidth=2)
    
    plt.xscale('log')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.title(r'$\mu$ = ' + f'{muList[i]}', fontsize=12)

    if i == 0:
        plt.ylabel('Entropy Production Rate', fontsize=12)
        handles_red = [line_h1[0], line_h2[0], line_h3[0]]
        
    
    ##### COOLING
    ax_c = plt.subplot(2, 4, i + 5) # 2 linhas, 4 colunas, posição 1 a 4

    line_c1 = plt.plot(tlist, Sprod_list[i][1], color='blue', linestyle=symbols[0], linewidth=2)
    line_c2 = plt.plot(tlist, Sprod_p_list[i][1], color='blue', linestyle=symbols[1], linewidth=2)
    line_c3 = plt.plot(tlist, Sprod_e_list[i][1], color='blue', linestyle=symbols[2], linewidth=2)
    
    plt.xscale('log')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.title(r'$\mu$ = ' + f'{muList[i]}', fontsize=12)
    plt.xlabel('Time', fontsize=12)
    
    if i == 0:
        plt.ylabel('Entropy Production Rate', fontsize=12)
    
    if i == 0:
        handles_blue = [line_c1[0], line_c2[0], line_c3[0]]


# Legenda Heating (Vermelha) - Superior
leg_h = fig.legend(handles_red, labels, 
                   loc='lower center', 
                   ncol=3, 
                   title="Heating", 
                   title_fontproperties={'weight':'bold', 'size':12},
                   fontsize=12,
                   bbox_to_anchor=(0.5, 0.08), 
                   frameon=False)

# Legenda Cooling (Azul) - Inferior
leg_c = fig.legend(handles_blue, labels, 
                   loc='lower center', 
                   ncol=3, 
                   title="Cooling", 
                   title_fontproperties={'weight':'bold', 'size':12},
                   fontsize=12,
                   bbox_to_anchor=(0.5, 0.01), 
                   frameon=False)

# Ajusta o layout para os gráficos não baterem na legenda
plt.tight_layout()
plt.subplots_adjust(bottom=0.25) # Abre espaço no fundo para as duas legendas

plt.show()
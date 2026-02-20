import numpy as np
import matplotlib.pyplot as plt


muList = [0.1, 0.5, 1, 2]
symbols = ['-', '--', ':']

## Leitura Arquivos

Iw_p_list = []
Iw_e_list = []
Kevol_p_list = []
Kevol_e_list = []
Sprod_p_list = []
Sprod_e_list = []

Iw_list = []
Vw_list = []
Lw_list = []
completion_list = []
Kevol_list = []
Sprod_list = []

for k, mu in enumerate(muList):
    
    tlist, Iw_p, Iw_e, Kevol_p, Kevol_e, Sprod_p, Sprod_e = np.loadtxt(f'./ThermalKinematics/mu{mu}-heating-contributions.txt', unpack=True)
    
    
    Iw_p_list.append([Iw_p])
    Iw_e_list.append([Iw_e])
    Kevol_p_list.append([Kevol_p])
    Kevol_e_list.append([Kevol_e])
    Sprod_p_list.append([Sprod_p])
    Sprod_e_list.append([Sprod_e])
    
    tlist, Iw_p, Iw_e, Kevol_p, Kevol_e, Sprod_p, Sprod_e = np.loadtxt(f'./ThermalKinematics/mu{mu}-cooling-contributions.txt', unpack=True)
    
    Iw_p_list[k].append(Iw_p)
    Iw_e_list[k].append(Iw_e)
    Kevol_p_list[k].append(Kevol_p)
    Kevol_e_list[k].append(Kevol_e)
    Sprod_p_list[k].append(Sprod_p)
    Sprod_e_list[k].append(Sprod_e)
    
    
    Iw, Kevol, Sprod = np.loadtxt(f'./ThermalKinematics/mu{mu}-heating.txt', unpack=True, usecols=(1,5,6))
    
    Iw_list.append([Iw])
    Kevol_list.append([Kevol])
    Sprod_list.append([Sprod])
    
    Iw, Kevol, Sprod = np.loadtxt(f'./ThermalKinematics/mu{mu}-cooling.txt', unpack=True, usecols=(1,5,6))
    
    Iw_list[k].append(Iw)
    Kevol_list[k].append(Kevol)
    Sprod_list[k].append(Sprod)
    

## Plots

### Wigner Fisher Information

##### Heating

plt.figure(figsize=(10,8))

for i in range(len(muList)):
    
    num = 221 + i
    
    plt.subplot(num)
    plt.plot(tlist, Iw_list[i][0], color='red', linestyle=symbols[0], label='Total')
    plt.plot(tlist, Iw_p_list[i][0], color='red', linestyle=symbols[1], label='Passive')
    plt.plot(tlist, Iw_e_list[i][0], color='red', linestyle=symbols[2], label='Ergotropic')
    
    plt.xscale('log')
    plt.title(r'$\mu$ = '+f'{muList[i]}', fontsize=12)
    plt.ylabel('Wigner Fisher Information', fontsize=12)
    plt.xlabel('Time', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    
plt.suptitle('Heating', fontsize=14)
plt.tight_layout()
plt.show()

##### Cooling

plt.figure(figsize=(10,8))

for i in range(len(muList)):
    
    num = 221 + i
    
    plt.subplot(num)
    plt.plot(tlist, Iw_list[i][1], color='blue', linestyle=symbols[0], label='Total')
    plt.plot(tlist, Iw_p_list[i][1], color='blue', linestyle=symbols[1], label='Passive')
    plt.plot(tlist, Iw_e_list[i][1], color='blue', linestyle=symbols[2], label='Ergotropic')
    
    plt.xscale('log')
    plt.title(r'$\mu$ = '+f'{muList[i]}', fontsize=12)
    plt.ylabel('Wigner Fisher Information', fontsize=12)
    plt.xlabel('Time', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    
plt.suptitle('Cooling', fontsize=14)
plt.tight_layout()
plt.show()


### Relative Entropy

##### Heating

plt.figure(figsize=(10,8))

for i in range(len(muList)):
    
    num = 221 + i
    
    plt.subplot(num)
    plt.plot(tlist, Kevol_list[i][0], color='red', linestyle=symbols[0], label='Total')
    plt.plot(tlist, Kevol_p_list[i][0], color='red', linestyle=symbols[1], label='Passive')
    plt.plot(tlist, Kevol_e_list[i][0], color='red', linestyle=symbols[2], label='Ergotropic')
    
    plt.xscale('log')
    plt.title(r'$\mu$ = '+f'{muList[i]}', fontsize=12)
    plt.ylabel('Relative Entropy', fontsize=12)
    plt.xlabel('Time', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    
plt.suptitle('Heating', fontsize=14)
plt.tight_layout()
plt.show()

##### Cooling

plt.figure(figsize=(10,8))

for i in range(len(muList)):
    
    num = 221 + i
    
    plt.subplot(num)
    plt.plot(tlist, Kevol_list[i][1], color='blue', linestyle=symbols[0], label='Total')
    plt.plot(tlist, Kevol_p_list[i][1], color='blue', linestyle=symbols[1], label='Passive')
    plt.plot(tlist, Kevol_e_list[i][1], color='blue', linestyle=symbols[2], label='Ergotropic')
    
    plt.xscale('log')
    plt.title(r'$\mu$ = '+f'{muList[i]}', fontsize=12)
    plt.ylabel('Relative Entropy', fontsize=12)
    plt.xlabel('Time', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    
plt.suptitle('Cooling', fontsize=14)
plt.tight_layout()
plt.show()

### Entropy Production

##### Heating

plt.figure(figsize=(10,8))

for i in range(len(muList)):
    
    num = 221 + i
    
    plt.subplot(num)
    plt.plot(tlist, Sprod_list[i][0], color='red', linestyle=symbols[0], label='Total')
    plt.plot(tlist, Sprod_p_list[i][0], color='red', linestyle=symbols[1], label='Passive')
    plt.plot(tlist, Sprod_e_list[i][0], color='red', linestyle=symbols[2], label='Ergotropic')
    
    plt.xscale('log')
    plt.title(r'$\mu$ = '+f'{muList[i]}', fontsize=12)
    plt.ylabel('Entropy Production', fontsize=12)
    plt.xlabel('Time', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    
plt.suptitle('Heating', fontsize=14)
plt.tight_layout()
plt.show()

##### Cooling

plt.figure(figsize=(10,8))

for i in range(len(muList)):
    
    num = 221 + i
    
    plt.subplot(num)
    plt.plot(tlist, Sprod_list[i][1], color='blue', linestyle=symbols[0], label='Total')
    plt.plot(tlist, Sprod_p_list[i][1], color='blue', linestyle=symbols[1], label='Passive')
    plt.plot(tlist, Sprod_e_list[i][1], color='blue', linestyle=symbols[2], label='Ergotropic')
    
    plt.xscale('log')
    plt.title(r'$\mu$ = '+f'{muList[i]}', fontsize=12)
    plt.ylabel('Entropy Production', fontsize=12)
    plt.xlabel('Time', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)
    
plt.suptitle('Cooling', fontsize=14)
plt.tight_layout()
plt.show()


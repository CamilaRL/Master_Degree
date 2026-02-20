import numpy as np
import matplotlib.pyplot as plt


muList = [0.1, 0.5, 1, 2]
symbols = ['-', '--', '-.', ':']

## Leitura Arquivos

Iw_list = []
Vw_list = []
Lw_list = []
completion_list = []
Kevol_list = []
Sprod_list = []

for k, mu in enumerate(muList):
    
    tlist, Iw, Vw, Lw, completion, Kevol, Sprod = np.loadtxt(f'./ThermalKinematics/mu{mu}-heating.txt', unpack=True)
    
    Iw_list.append([Iw])
    Vw_list.append([Vw])
    Lw_list.append([Lw])
    completion_list.append([completion])
    Kevol_list.append([Kevol])
    Sprod_list.append([Sprod])
    
    tlist, Iw, Vw, Lw, completion, Kevol, Sprod = np.loadtxt(f'./ThermalKinematics/mu{mu}-cooling.txt', unpack=True)
    
    Iw_list[k].append(Iw)
    Vw_list[k].append(Vw)
    Lw_list[k].append(Lw)
    completion_list[k].append(completion)
    Kevol_list[k].append(Kevol)
    Sprod_list[k].append(Sprod)


## Plots

### Wigner Fisher Information

plt.figure(figsize=(10,5))

for i in range(len(muList)):
    
    plt.subplot(121)
    plt.plot(tlist, Iw_list[i][0], color='red', linestyle=symbols[i], label=r'$\mu$ = '+f'{muList[i]}')
    
    plt.subplot(122)
    plt.plot(tlist, Iw_list[i][1], color='blue', linestyle=symbols[i], label=r'$\mu$ = '+f'{muList[i]}')

plt.subplot(121)
plt.xscale('log')
plt.ylabel('Wigner Fisher Information', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.title('Heating', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.subplot(122)
plt.xscale('log')
plt.xlabel('Time', fontsize=12)
plt.title('Cooling', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.tight_layout()
plt.show()

### Velocity

plt.figure(figsize=(10,5))

for i in range(len(muList)):
    
    plt.subplot(121)
    plt.plot(tlist, Vw_list[i][0], color='red', linestyle=symbols[i], label=r'$\mu$ = '+f'{muList[i]}')
    
    plt.subplot(122)
    plt.plot(tlist, Vw_list[i][1], color='blue', linestyle=symbols[i], label=r'$\mu$ = '+f'{muList[i]}')

plt.subplot(121)
#plt.xscale('log')
plt.ylabel('Statistical Velocity', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.title('Heating', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.subplot(122)
#plt.xscale('log')
plt.xlabel('Time', fontsize=12)
plt.title('Cooling', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.tight_layout()
plt.show()

### Distance

plt.figure(figsize=(10,5))

for i in range(len(muList)):
    
    plt.subplot(121)
    plt.plot(tlist, Lw_list[i][0], color='red', linestyle=symbols[i], label=r'$\mu$ = '+f'{muList[i]}')
    
    plt.subplot(122)
    plt.plot(tlist, Lw_list[i][1], color='blue', linestyle=symbols[i], label=r'$\mu$ = '+f'{muList[i]}')

plt.subplot(121)
plt.ylabel('Statistical Distance', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.title('Heating', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.subplot(122)
plt.xlabel('Time', fontsize=12)
plt.title('Cooling', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.tight_layout()
plt.show()


### Completion

plt.figure(figsize=(10,5))

for i in range(len(muList)):
    
    plt.subplot(121)
    plt.plot(tlist[250:300], completion_list[i][0][250:300], color='red', linestyle=symbols[i], label=r'$\mu$ = '+f'{muList[i]}')
    
    plt.subplot(122)
    plt.plot(tlist[250:300], completion_list[i][1][250:300], color='blue', linestyle=symbols[i], label=r'$\mu$ = '+f'{muList[i]}')

plt.subplot(121)
plt.ylabel('Degree of Completion', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.title('Heating', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.xlim(left=25, right=30)

plt.subplot(122)
plt.xlabel('Time', fontsize=12)
plt.title('Cooling', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.xlim(left=25, right=30)

plt.tight_layout()
plt.show()

plt.figure(figsize=(10,8))

for i in range(len(muList)):
    
    num = 221 + i
    
    plt.subplot(num)
    plt.plot(tlist, completion_list[i][0], color='red', linestyle=symbols[i], label='Heating')
    plt.plot(tlist, completion_list[i][1], color='blue', linestyle=symbols[i], label='Cooling')
    
    plt.title(r'$\mu$ = '+f'{muList[i]}', fontsize=12)
    plt.ylabel('Degree of Completion', fontsize=12)
    plt.xlabel('Time', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=12)

plt.tight_layout()
plt.show()

### Relative Entropy

plt.figure(figsize=(10,5))

for i in range(len(muList)):
    
    plt.subplot(121)
    plt.plot(tlist, Kevol_list[i][0], color='red', linestyle=symbols[i], label=r'$\mu$ = '+f'{muList[i]}')
    
    plt.subplot(122)
    plt.plot(tlist, Kevol_list[i][1], color='blue', linestyle=symbols[i], label=r'$\mu$ = '+f'{muList[i]}')

plt.subplot(121)
plt.xscale('log')
plt.ylabel('Relative Entropy', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.title('Heating', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.subplot(122)
plt.xscale('log')
plt.xlabel('Time', fontsize=12)
plt.title('Cooling', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.tight_layout()
plt.show()

### Entropy Production

plt.figure(figsize=(10,5))

for i in range(len(muList)):
    
    plt.subplot(121)
    plt.plot(tlist, Sprod_list[i][0], color='red', linestyle=symbols[i], label=r'$\mu$ = '+f'{muList[i]}')
    
    plt.subplot(122)
    plt.plot(tlist, Sprod_list[i][1], color='blue', linestyle=symbols[i], label=r'$\mu$ = '+f'{muList[i]}')

plt.subplot(121)
plt.xscale('log')
plt.ylabel('Entropy Production', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.title('Heating', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.subplot(122)
plt.xscale('log')
plt.xlabel('Time', fontsize=12)
plt.title('Cooling', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)

plt.tight_layout()
plt.show()







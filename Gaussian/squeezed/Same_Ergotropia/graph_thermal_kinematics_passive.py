import numpy as np
import matplotlib.pyplot as plt

## rc = cooling | rh = heating
rcList, rhList = np.loadtxt('./ThermalKinematics/initial_temperatures.txt', unpack=True, usecols=(0, 2))
symbols = ['-', '--', '-.', ':']

## Leitura Arquivos

Iw_list = []
Vw_list = []
Lw_list = []
completion_list = []

for k in range(2):
    
    tlist, Iw, Vw, Lw, completion = np.loadtxt(f'./ThermalKinematics/r{rhList[k]}-heating-passive.txt', unpack=True)
    
    Iw_list.append([Iw])
    Vw_list.append([Vw])
    Lw_list.append([Lw])
    completion_list.append([completion])
    
    
    tlist, Iw, Vw, Lw, completion = np.loadtxt(f'./ThermalKinematics/r{rcList[k]}-cooling-passive.txt', unpack=True)
    
    Iw_list[k].append(Iw)
    Vw_list[k].append(Vw)
    Lw_list[k].append(Lw)
    completion_list[k].append(completion)


## Plots TOTAL

### Wigner Fisher Information

plt.figure(figsize=(10,5))

for i in range(2):
    
    plt.subplot(121)
    plt.plot(tlist, Iw_list[i][0], color='red', linestyle=symbols[i], linewidth=2, label=r'$r$ = '+f'{rhList[i]:.2f}')
    
    plt.subplot(122)
    plt.plot(tlist, Iw_list[i][1], color='blue', linestyle=symbols[i], linewidth=2, label=r'$r$ = '+f'{rcList[i]:.2f}')

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

for i in range(2):
    
    plt.subplot(121)
    plt.plot(tlist, Vw_list[i][0], color='red', linestyle=symbols[i], linewidth=2, label=r'$r$ = '+f'{rhList[i]:.2f}')
    
    plt.subplot(122)
    plt.plot(tlist, Vw_list[i][1], color='blue', linestyle=symbols[i], linewidth=2, label=r'$r$ = '+f'{rcList[i]:.2f}')

plt.subplot(121)
plt.ylabel('Statistical Velocity', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.title('Heating', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.xlim(right=50)
plt.yscale('log')

plt.subplot(122)
plt.xlabel('Time', fontsize=12)
plt.title('Cooling', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.xlim(right=80)
plt.yscale('log')

plt.tight_layout()
plt.show()

### Distance

plt.figure(figsize=(10,5))

for i in range(2):
    
    plt.subplot(121)
    plt.plot(tlist, Lw_list[i][0], color='red', linestyle=symbols[i], linewidth=2, label=r'$r$ = '+f'{rhList[i]:.2f}')
    
    plt.subplot(122)
    plt.plot(tlist, Lw_list[i][1], color='blue', linestyle=symbols[i], linewidth=2, label=r'$r$ = '+f'{rcList[i]:.2f}')

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

#plt.figure(figsize=(10,8))
plt.figure(figsize=(10,5))

for i in range(2):
    
    num = 121 + i #221 + i
    
    plt.subplot(num)
    plt.plot(tlist, completion_list[i][0], color='red', linestyle=symbols[i], label=r'$W_c$ - $r$ = '+f'{rhList[i]:.2f} ')
    plt.plot(tlist, completion_list[i][1], color='blue', linestyle=symbols[i], linewidth=2, label=r'$W_h$ - $r$ = '+f'{rcList[i]:.2f} ')
    plt.legend(fontsize=12)
    
    if i == 0:
        plt.ylabel('Degree of Completion', fontsize=12)
    
    plt.xlabel('Time', fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlim(right=80)

plt.tight_layout()
plt.show()


plt.plot(tlist, completion_list[0][0], color='red', linestyle=symbols[0], label=r'$r_c$ = '+f'{rhList[0]:.2f} ')
plt.plot(tlist, completion_list[0][1], color='blue', linestyle=symbols[0], linewidth=2, label=r'$r_h$ = '+f'{rcList[0]:.2f} ')
plt.plot(tlist, completion_list[1][0], color='red', linestyle=symbols[1], label=r'$r_c$ = '+f'{rhList[1]:.2f} ')
plt.plot(tlist, completion_list[1][1], color='blue', linestyle=symbols[1], linewidth=2, label=r'$r_h$ = '+f'{rcList[1]:.2f} ')

#plt.ylabel('Degree of Completion', fontsize=12)
plt.xlabel('Time', fontsize=12)
plt.title('Passive State', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(fontsize=12)
plt.xlim(right=80)

plt.show()

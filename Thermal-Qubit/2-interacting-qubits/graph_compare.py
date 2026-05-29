import numpy as np
import matplotlib.pyplot as plt
import os

def Read_Completion(dir, c, g):

    path = f'{dir}/completion_c{c}_g{g}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion


def Read_Velocidade(c, modo, g):

    path = f'./ThermalKinematics/{modo}/velocity_c{c}_g{g}.txt'
    
    tlist, completion = np.loadtxt(path, unpack=True)
    
    return tlist, completion




### MAIN ###

modoList = ['Cooling', 'Heating']

titulo = ['Zero Initial Coherence', 'Maximum Initial Coherence']

cList = ['min', 'max']

cores = ['blue', 'red']

J = [0.0, 0.8, 1.6]

# Leitura dos dados
## Sem EOF
tlist_r_J0_cmin, completion_r_J0_cmin = Read_Completion(f'./XX-model-noEoF/ThermalKinematics/{modoList[0]}', cList[0], J[0])
tlist_a_J0_cmin, completion_a_J0_cmin = Read_Completion(f'./XX-model-noEoF/ThermalKinematics/{modoList[1]}', cList[0], J[0])
tlist_r_J0_cmax, completion_r_J0_cmax = Read_Completion(f'./XX-model-noEoF/ThermalKinematics/{modoList[0]}', cList[1], J[0])
tlist_a_J0_cmax, completion_a_J0_cmax = Read_Completion(f'./XX-model-noEoF/ThermalKinematics/{modoList[1]}', cList[1], J[0])

tlist_r_J1_cmin, completion_r_J1_cmin = Read_Completion(f'./XX-model-noEoF/ThermalKinematics/{modoList[0]}', cList[0], J[1])
tlist_a_J1_cmin, completion_a_J1_cmin = Read_Completion(f'./XX-model-noEoF/ThermalKinematics/{modoList[1]}', cList[0], J[1])
tlist_r_J1_cmax, completion_r_J1_cmax = Read_Completion(f'./XX-model-noEoF/ThermalKinematics/{modoList[0]}', cList[1], J[1])
tlist_a_J1_cmax, completion_a_J1_cmax = Read_Completion(f'./XX-model-noEoF/ThermalKinematics/{modoList[1]}', cList[1], J[1])

## Com EOF (J=1.6 mas g=0.8)
tlist_rE_J2_cmin, completion_rE_J2_cmin = Read_Completion(f'./XX-model/ThermalKinematics/{modoList[0]}', cList[0], J[1])
tlist_aE_J2_cmin, completion_aE_J2_cmin = Read_Completion(f'./XX-model/ThermalKinematics/{modoList[1]}', cList[0], J[1])
tlist_rE_J2_cmax, completion_rE_J2_cmax = Read_Completion(f'./XX-model/ThermalKinematics/{modoList[0]}', cList[1], J[1])
tlist_aE_J2_cmax, completion_aE_J2_cmax = Read_Completion(f'./XX-model/ThermalKinematics/{modoList[1]}', cList[1], J[1])



############################################################################
## heating and cooling =  cmin

fig = plt.figure(figsize=(10,5))

plt.subplot(121)

a_J0_cmin = plt.plot(tlist_a_J0_cmin, completion_a_J0_cmin, color='red', linestyle=':', linewidth=2, label=f'J = {J[0]:.1f}')
a_J1_cmin = plt.plot(tlist_a_J1_cmin, completion_a_J1_cmin, color='red', linestyle='--', linewidth=2, label=f'J = {J[1]:.1f}')
aE_J2_cmin = plt.plot(tlist_aE_J2_cmin, completion_aE_J2_cmin, color='red', linestyle='-', linewidth=2, label=f'J = {J[2]:.1f}')

r_J0_cmin = plt.plot(tlist_r_J0_cmin, completion_r_J0_cmin, color='blue', linestyle=':', linewidth=2, label=f'J = {J[0]:.1f}')
r_J1_cmin = plt.plot(tlist_r_J1_cmin, completion_r_J1_cmin, color='blue', linestyle='--', linewidth=2, label=f'J = {J[1]:.1f}')
rE_J2_cmin = plt.plot(tlist_rE_J2_cmin, completion_rE_J2_cmin, color='blue', linestyle='-', linewidth=2, label=f'J = {J[2]:.1f}')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Degree of Completion', fontsize=12)
plt.title(f'{titulo[0]}', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)


plt.subplot(122)

a_J0_cmax = plt.plot(tlist_a_J0_cmax, completion_a_J0_cmax, color='red', linestyle=':', linewidth=2, label=f'J = {J[0]:.1f}')
a_J1_cmax = plt.plot(tlist_a_J1_cmax, completion_a_J1_cmax, color='red', linestyle='--', linewidth=2, label=f'J = {J[1]:.1f}')
aE_J2_cmax = plt.plot(tlist_aE_J2_cmax, completion_aE_J2_cmax, color='red', linestyle='-', linewidth=2, label=f'J = {J[2]:.1f}')

r_J0_cmax = plt.plot(tlist_r_J0_cmax, completion_r_J0_cmax, color='blue', linestyle=':', linewidth=2, label=f'J = {J[0]:.1f}')
r_J1_cmax = plt.plot(tlist_r_J1_cmax, completion_r_J1_cmax, color='blue', linestyle='--', linewidth=2, label=f'J = {J[1]:.1f}')
rE_J2_cmax = plt.plot(tlist_rE_J2_cmax, completion_rE_J2_cmax, color='blue', linestyle='-', linewidth=2, label=f'J = {J[2]:.1f}')

plt.xlabel('Time', fontsize=12)
plt.title(f'{titulo[1]}', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)

labels = [f'J = {J[0]:.1f}', f'J = {J[1]:.1f}', f'J = {J[2]:.1f}']
handles_red = [a_J0_cmin[0], a_J1_cmin[0], aE_J2_cmin[0]]
handles_blue = [r_J0_cmin[0], r_J1_cmin[0], rE_J2_cmin[0]]

# Legenda Heating (Vermelha) - Superior
leg_h = fig.legend(handles_red, labels, 
                   loc='lower center', 
                   ncol=3, 
                   title='Heating', 
                   title_fontproperties={'weight':'bold', 'size':12},
                   fontsize=12,
                   bbox_to_anchor=(0.5, 0.1), 
                   frameon=False)

# Legenda Cooling (Azul) - Inferior
leg_c = fig.legend(handles_blue, labels, 
                   loc='lower center', 
                   ncol=3, 
                   title='Cooling', 
                   title_fontproperties={'weight':'bold', 'size':12},
                   fontsize=12,
                   bbox_to_anchor=(0.5, 0.01), 
                   frameon=False)

plt.tight_layout()
plt.subplots_adjust(bottom=0.35)
plt.show()







############################################################################
## heating and cooling =  cmin
fig = plt.figure(figsize=(10,5))

plt.subplot(121)

plt.plot(tlist_a_J0_cmin, completion_a_J0_cmin, color='red', linestyle='--', linewidth=2, label=f'J = {J[0]:.1f}')
plt.plot(tlist_a_J1_cmin, completion_a_J1_cmin, color='red', linestyle='-', linewidth=2, label=f'J = {J[1]:.1f}')
plt.plot(tlist_aE_J2_cmin, completion_aE_J2_cmin, color='red', linestyle=':', linewidth=2, label=f'J = {J[2]:.1f}')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Degree of Completion', fontsize=12)
plt.title('Heating', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)

plt.subplot(122)

plt.plot(tlist_r_J0_cmin, completion_r_J0_cmin, color='blue', linestyle='--', linewidth=2, label=f'J = {J[0]:.1f}')
plt.plot(tlist_r_J1_cmin, completion_r_J1_cmin, color='blue', linestyle='-', linewidth=2, label=f'J = {J[1]:.1f}')
plt.plot(tlist_rE_J2_cmin, completion_rE_J2_cmin, color='blue', linestyle=':', linewidth=2, label=f'J = {J[2]:.1f}')

plt.xlabel('Time', fontsize=12)
plt.title('Cooling', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)

plt.suptitle(f'{titulo[0]}', fontsize=14)
plt.tight_layout()
plt.show()


## heating and cooling =  cmax
fig = plt.figure(figsize=(10,5))

plt.subplot(121)

plt.plot(tlist_a_J0_cmax, completion_a_J0_cmax, color='red', linestyle='--', linewidth=2, label=f'J = {J[0]:.1f}')
plt.plot(tlist_a_J1_cmax, completion_a_J1_cmax, color='red', linestyle='-', linewidth=2, label=f'J = {J[1]:.1f}')
plt.plot(tlist_aE_J2_cmax, completion_aE_J2_cmax, color='red', linestyle=':', linewidth=2, label=f'J = {J[2]:.1f}')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Degree of Completion', fontsize=12)
plt.title('Heating', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)

plt.subplot(122)

plt.plot(tlist_r_J0_cmax, completion_r_J0_cmax, color='blue', linestyle='--', linewidth=2, label=f'J = {J[0]:.1f}')
plt.plot(tlist_r_J1_cmax, completion_r_J1_cmax, color='blue', linestyle='-', linewidth=2, label=f'J = {J[1]:.1f}')
plt.plot(tlist_rE_J2_cmax, completion_rE_J2_cmax, color='blue', linestyle=':', linewidth=2, label=f'J = {J[2]:.1f}')

plt.xlabel('Time', fontsize=12)
plt.title('Cooling', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)

plt.suptitle(f'{titulo[1]}', fontsize=14)
plt.tight_layout()
plt.show()

############################################################################
## cmin and cmax = Heating
fig = plt.figure(figsize=(10,5))

plt.subplot(121)

plt.plot(tlist_a_J0_cmin, completion_a_J0_cmin, color='red', linestyle='--', linewidth=2, label=f'J = {J[0]:.1f}')
plt.plot(tlist_a_J1_cmin, completion_a_J1_cmin, color='red', linestyle='-', linewidth=2, label=f'J = {J[1]:.1f}')
plt.plot(tlist_aE_J2_cmin, completion_aE_J2_cmin, color='red', linestyle=':', linewidth=2, label=f'J = {J[2]:.1f}')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Degree of Completion', fontsize=12)
plt.title(f'{titulo[0]}', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)

plt.subplot(122)

plt.plot(tlist_a_J0_cmax, completion_a_J0_cmax, color='red', linestyle='--', linewidth=2, label=f'J = {J[0]:.1f}')
plt.plot(tlist_a_J1_cmax, completion_a_J1_cmax, color='red', linestyle='-', linewidth=2, label=f'J = {J[1]:.1f}')
plt.plot(tlist_aE_J2_cmax, completion_aE_J2_cmax, color='red', linestyle=':', linewidth=2, label=f'J = {J[2]:.1f}')

plt.xlabel('Time', fontsize=12)
plt.title(f'{titulo[1]}', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)

plt.suptitle(f'{modoList[1]}', fontsize=14)
plt.tight_layout()
plt.show()


## cmin and cmax = Cooling
fig = plt.figure(figsize=(10,5))

plt.subplot(121)

plt.plot(tlist_r_J0_cmin, completion_r_J0_cmin, color='blue', linestyle='--', linewidth=2, label=f'J = {J[0]:.1f}')
plt.plot(tlist_r_J1_cmin, completion_r_J1_cmin, color='blue', linestyle='-', linewidth=2, label=f'J = {J[1]:.1f}')
plt.plot(tlist_rE_J2_cmin, completion_rE_J2_cmin, color='blue', linestyle=':', linewidth=2, label=f'J = {J[2]:.1f}')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Degree of Completion', fontsize=12)
plt.title(f'{titulo[0]}', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)

plt.subplot(122)

plt.plot(tlist_r_J0_cmax, completion_r_J0_cmax, color='blue', linestyle='--', linewidth=2, label=f'J = {J[0]:.1f}')
plt.plot(tlist_r_J1_cmax, completion_r_J1_cmax, color='blue', linestyle='-', linewidth=2, label=f'J = {J[1]:.1f}')
plt.plot(tlist_rE_J2_cmax, completion_rE_J2_cmax, color='blue', linestyle=':', linewidth=2, label=f'J = {J[2]:.1f}')

plt.xlabel('Time', fontsize=12)
plt.title(f'{titulo[1]}', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)

plt.suptitle(f'{modoList[0]}', fontsize=14)
plt.tight_layout()
plt.show()


############################################################################
## J assimetria = cmin
fig = plt.figure(figsize=(10,5))

plt.subplot(131)

plt.plot(tlist_a_J0_cmin, completion_a_J0_cmin, color='red', linestyle='-', linewidth=2, label=f'Heating')
plt.plot(tlist_r_J0_cmin, completion_r_J0_cmin, color='blue', linestyle='-', linewidth=2, label=f'Cooling')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Degree of Completion', fontsize=12)
plt.title(f'J = {J[0]:.1f}', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)

plt.subplot(132)

plt.plot(tlist_a_J1_cmin, completion_a_J1_cmin, color='red', linestyle='-', linewidth=2, label=f'Heating')
plt.plot(tlist_r_J1_cmin, completion_r_J1_cmin, color='blue', linestyle='-', linewidth=2, label=f'Cooling')

plt.xlabel('Time', fontsize=12)
plt.title(f'J = {J[1]:.1f}', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)

plt.subplot(133)

plt.plot(tlist_aE_J2_cmin, completion_aE_J2_cmin, color='red', linestyle='-', linewidth=2, label=f'Heating')
plt.plot(tlist_rE_J2_cmin, completion_rE_J2_cmin, color='blue', linestyle='-', linewidth=2, label=f'Cooling')

plt.xlabel('Time', fontsize=12)
plt.title(f'J = {J[2]:.1f}', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)


plt.suptitle(f'{titulo[0]}', fontsize=14)
plt.tight_layout()
plt.show()



## J assimetria = cmax
fig = plt.figure(figsize=(10,5))

plt.subplot(131)

plt.plot(tlist_a_J0_cmax, completion_a_J0_cmax, color='red', linestyle='-', linewidth=2, label=f'Heating')
plt.plot(tlist_r_J0_cmax, completion_r_J0_cmax, color='blue', linestyle='-', linewidth=2, label=f'Cooling')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Degree of Completion', fontsize=12)
plt.title(f'J = {J[0]:.1f}', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)

plt.subplot(132)

plt.plot(tlist_a_J1_cmax, completion_a_J1_cmax, color='red', linestyle='-', linewidth=2, label=f'Heating')
plt.plot(tlist_r_J1_cmax, completion_r_J1_cmax, color='blue', linestyle='-', linewidth=2, label=f'Cooling')

plt.xlabel('Time', fontsize=12)
plt.title(f'J = {J[1]:.1f}', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)

plt.subplot(133)

plt.plot(tlist_aE_J2_cmax, completion_aE_J2_cmax, color='red', linestyle='-', linewidth=2, label=f'Heating')
plt.plot(tlist_rE_J2_cmax, completion_rE_J2_cmax, color='blue', linestyle='-', linewidth=2, label=f'Cooling')

plt.xlabel('Time', fontsize=12)
plt.title(f'J = {J[2]:.1f}', fontsize=12)
plt.legend(loc='lower right', fontsize=12)
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xscale('log')
plt.xlim(left=0.1)


plt.suptitle(f'{titulo[1]}', fontsize=14)
plt.tight_layout()
plt.show()
import numpy as np
import matplotlib.pyplot as plt


g = 0.8

cmodlist = np.loadtxt(f'./g-{g}/DensityMatrices/cmod.txt', unpack=True)
cmodlist.sort()

Slist1 = []
Slist2 = []

for c in cmodlist:
    
    tlist1, S1 = np.loadtxt(f'g-{g}/Entropy/entropy1-{c:.3f}.txt', unpack=True)
    tlist2, S2 = np.loadtxt(f'g-{g}/Entropy/entropy2-{c:.3f}.txt', unpack=True)

import matplotlib.pyplot as plt
import numpy as np

zeig, xeig, qreal, qimag = np.loadtxt('./KDQ_SG.txt', unpack=True)


pairs = np.arange(0, len(zeig), 1)

pairs_leg = []

for i in range(len(zeig)):

	pairs_leg.append(f'({zeig[i]}, {xeig[i]})')

fig = plt.figure(figsize=12)

plt.subplot(11)
plt.

#!/usr/bin/env python
import glob
import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

modes = ['0S2','0S3','0S4','0S5']
letters = ['(a)','(b)','(c)','(d)']
fig, ax = plt.subplots(2, 2, figsize=(16,12))
ax = ax.flatten()


for midx, mode in enumerate(modes):
    stas, sign, eve, colors, symbols, sensors = [],[],[],[],[], []
    files = glob.glob('*' + mode + '*.csv')

    for cfile in files:
        f = open(cfile, 'r')
        jday = cfile.split('_')[2]
        print(jday)
        for line in f:
            line = line.split(',')
            stas.append(line[0])
            snr = 20*np.log10(float(line[-2]))
            if snr > 20:
                 sign.append(20)
            elif snr < -80:
                sign.append(-80)
            else:
                sign.append(snr)
            if 'STS-1' in line[-1]:
                c='C0'
            elif 'STS-6' in line[-1]:
                c='C1'
            elif 'KS' in line[-1]:
                c = 'C2'
            elif '360' in line[-1]:
                c = 'C3'
            else:
                c='C4'
            if '360' in line[-1]:
                sensors.append('T-360GSN')
            else:
                sensors.append(line[-1])
            colors.append(c)
            eve.append(jday)
            if jday == '91':
                symbols.append('s')
            else:
                symbols.append('o')
        f.close()
    

    sensorsuni = list(set(sensors))

    #stasuni = list(set(stas))
    #stasuni.sort()
    #juststas =[]
    #for sta in stasuni:
    #    juststas.append(sta.split('.')[1])
    #juststas.sort()
    sign = np.array(sign)
    eve = np.array(eve)
    for idx, sen  in enumerate(sensorsuni):
        gidx = [i for i, x in enumerate(sensors) if sen in x]
        sigval0, sigval1 = [],[]
        for gidx2 in gidx:
            if eve[gidx2] == '91':
                sigval0.append(sign[gidx2])
            else:
                sigval1.append(sign[gidx2])

            ax[midx].plot(sign[gidx2], idx, symbols[gidx2], color=colors[gidx2], alpha=0.5)
        sigval0 = np.array(sigval0)
        sigval1 = np.array(sigval1)
        ax[midx].plot(np.mean(sigval0), idx, 's', color=colors[gidx2], markersize=20, alpha=0.5, markeredgecolor='k')
        ax[midx].plot(np.mean(sigval1), idx, 'o', color=colors[gidx2], markersize=20, alpha=0.5, markeredgecolor='k')
    ax[midx].set_title(letters[midx],loc='left')
    ax[midx].set_yticks(range(len(sensorsuni)))
    ax[midx].set_yticklabels(sensorsuni, fontsize=12)
    #ax[midx].set_xlim((-80,20))
    #ax[midx].set_ylim((0.5, len(sensorsuni)))
    ax[midx].set_xlabel('Noise of $_' + mode[0] + mode[1] + '_' + mode[2] + '$ (dB)')


plt.savefig('figure4again.pdf',format='PDF', dpi=400)
plt.savefig('figure4again.png', format='PNG', dpi=400)

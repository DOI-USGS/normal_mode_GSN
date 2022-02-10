#!/usr/bin/env python
import glob
import sys
import matplotlib.pyplot as plt
import numpy as np
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
import matplotlib as mpl
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import matplotlib.font_manager
import cartopy as cart

# Importing and applying font
mpl.rc('font', family = 'serif')
mpl.rc('font', serif = 'Times') 
mpl.rc('text', usetex = True)
mpl.rc('font', size=16)

import warnings
warnings.filterwarnings("ignore")


files = ['ALL_SNRresults_2014_91.csv', 'ALL_SNRresults_2021_63.csv','ALL_SNRresults_2021_210.csv','ALL_SNRresults_2021_224.csv']





def setupmap(handle):
    handle.add_feature(cfeature.LAND)
    handle.add_feature(cfeature.OCEAN)
    handle.add_feature(cfeature.COASTLINE)
    handle.add_feature(cfeature.BORDERS, linestyle=':')
    handle.add_feature(cfeature.LAKES, alpha=0.5)
    handle.add_feature(cfeature.RIVERS)
    handle.add_feature(cfeature.STATES, edgecolor='gray')
    return handle


elats =[-19.61, -29.735, 55.474, -58.416]
elons = [-70.769, -177.282, -157.917, -25.321]
letters = ['(a)', '(b)', '(c)', '(d)']
fig = plt.figure(1, figsize=(10,8))

sens = ['STS-1', 'STS-6', 'T-360GSN', 'KS-54000', 'Other']
markers = ['^', 's', 'o', 'X', 'd']
ax = [0,0,0,0]

for idx in range(4):
    ax[idx] = fig.add_subplot(2,2,idx+1, projection = ccrs.Robinson())
    ax[idx].set_title(letters[idx], loc='left')
    ax[idx].coastlines()
    ax[idx].add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
    ax[idx].set_global()
    ax[idx].coastlines()


for idx2, sen in enumerate(sens):
    for idx, cfile in enumerate(files):
        f = open(cfile,'r')

        alat, alon, val = [], [], []
        alat2, alon2, val2 =[],[],[]
        for line in f:
            line = line.rstrip()
            line = line.split(',')
            if line[-1] == sen:
                pass
            else:
                continue
            if float(line[3]) >= 2.5:
                alat.append(float(line[1]))
                alon.append(float(line[2]))
                val.append(float(line[3]))
                print(line)
                if (idx == 0):
                    print(line[0])
            else:
                alat2.append(float(line[1]))
                alon2.append(float(line[2]))
                val2.append(float(line[3]))


        
        #ax=setupmap(ax)
        

        ax[idx].scatter(elons[idx], elats[idx], c='r',marker='*',s= 200, transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10)
        
        im = ax[idx].scatter(alon2, alat2, c= val2, s=50, transform=ccrs.Geodetic(), zorder=3, vmin=1, vmax=5, marker=markers[idx2], alpha=0.2, facecolor="None")
        im = ax[idx].scatter(alon, alat, c=val, s=100, transform=ccrs.Geodetic(), zorder=3, vmin=1., vmax=5, marker=markers[idx2], label=sen)
        #if idx == 3:
        #    print(ax[idx].get_legend_handles_labels())
        f.close()
        mval = np.mean(val)
        if (len(val) + len(val2)) == 0:
            pro = 'None'
        else:
            pro = 100*len(val)/(len(val) + len(val2))
        print(sen + ' ' + str(len(val)+len(val2)))

        print(cfile + ' ' + sen +  ' ' + str(len(alat)) + ' ' + str(mval) + ' ' + str(pro))



    #ax.stock_img()

#for ran in [0, 5, 10, 15]:
#    im2 = ax.scatter([0],[-90], c='k', alpha=0.3, s=2**(ran/2.5)+5, label=str(2.5*ran) + ' \% Uncertainty' , transform=ccrs.Geodetic() )
    

#for trip in zip(alot, alon, ap):
#    print(trip)
fig2 = plt.figure(2, figsize=(10,8))
ax1 = fig2.add_subplot(2,2,1, projection = ccrs.Robinson())
sens.append('Earthquake')
markers.append('*')

for idx2, sen in enumerate(sens):
    if idx2 <=4:
        im3 = ax1.scatter(0,0, c=4, s=100, transform=ccrs.Geodetic(), zorder=3, vmin=1, vmax=5, marker=markers[idx2], label=sen)
    else:
        im3 = ax1.scatter(0,0, c='r', s=200, transform=ccrs.Geodetic(), zorder=3, vmin=1., vmax=5, marker=markers[idx2], label=sen)
handles, labels = ax1.get_legend_handles_labels()
#print(handles)

ax[3].legend(handles, labels, bbox_to_anchor=(1.1,-0.05), ncol=6, fontsize=12)
#fig.legend(ncol=4, )  
plt.figure(1)
cb_ax = fig.add_axes([0.2, 0.05, 0.6, 0.03])
cbar = fig.colorbar(im, cax=cb_ax, orientation='horizontal')
plt.close(fig2)

cbar.set_label('Signal-to-Noise Ratio') 
#plt.tight_layout()
plt.savefig('Map_of_usage_ALL.png', bbox_inches='tight', format='PNG', dpi=400)
plt.savefig('Map_of_usage_ALL.pdf', bbox_inches='tight', format='PDF', dpi=400)
plt.show()

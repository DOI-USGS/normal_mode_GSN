#!/usr/bin/env python
import cartopy as cart
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime

# Importing and applying font
mpl.rc('font', family='serif')
mpl.rc('font', serif='Times')
mpl.rc('text', usetex=True)
mpl.rc('font', size=16)

import warnings

warnings.filterwarnings("ignore")
client = Client("IRIS")

files = ['ALL_SNRresults_2014_91.csv', 'ALL_SNRresults_2021_63.csv',
         'ALL_SNRresults_2021_210.csv', 'ALL_SNRresults_2021_224.csv']


def setupmap(handle):
    handle.add_feature(cfeature.LAND)
    handle.add_feature(cfeature.OCEAN)
    handle.add_feature(cfeature.COASTLINE)
    handle.add_feature(cfeature.BORDERS, linestyle=':')
    handle.add_feature(cfeature.LAKES, alpha=0.5)
    handle.add_feature(cfeature.RIVERS)
    handle.add_feature(cfeature.STATES, edgecolor='gray')
    return handle


eves = [UTCDateTime('2014-04-01T23:46:47'), UTCDateTime('2021-03-04T17:41:25'),
        UTCDateTime('2021-08-12T18:35:20'), UTCDateTime('2021-07-29T06:15:49')]

elats = [-19.61, -29.735, 55.474, -58.416]
elons = [-70.769, -177.282, -157.917, -25.321]
letters = ['(a) April 1, 2014', '(b) March 4, 2021', '(c) July 29, 2021',
           '(d) August 21, 2021']
fig = plt.figure(1, figsize=(10, 8))

sens = ['STS-1', 'STS-6', 'Trillium 360', 'KS-54000', 'Other']
markers = ['^', 's', 'o', 'X', 'd']
ax = [0, 0, 0, 0]

for idx in range(4):
    ax[idx] = fig.add_subplot(2, 2, idx + 1, projection=ccrs.Robinson())
    ax[idx].set_title(letters[idx], loc='left')
    ax[idx].coastlines()
    ax[idx].add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
    ax[idx].set_global()
    ax[idx].coastlines()

for idx2, sen in enumerate(sens):

    for idx, cfile in enumerate(files):

        f = open(cfile, 'r')

        snclsgood = []
        alat, alon, val = [], [], []
        alat2, alon2, val2 = [], [], []
        for line in f:
            line = line.rstrip()
            line = line.split(',')
            if line[-1] == sen:
                snclsgood.append(line[0])
                pass

            else:
                continue

            if float(line[3]) >= 2.5:
                alat.append(float(line[1]))
                alon.append(float(line[2]))
                val.append(float(line[3]))

                if (idx == 0):
                    print(line[0])
            elif float(line[3]) >= 0:
                alat2.append(float(line[1]))
                alon2.append(float(line[2]))
                val2.append(float(line[3]))

        # whatisleft = list(set(sncls) - set(snclsgood))

        ax[idx].scatter(elons[idx], elats[idx], c='r', marker='*', s=200,
                        transform=ccrs.Geodetic(), zorder=3, vmin=0., vmax=10)

        im = ax[idx].scatter(alon2, alat2, c=val2, s=50,
                             transform=ccrs.Geodetic(), zorder=3, vmin=1,
                             vmax=5, marker=markers[idx2], alpha=0.3,
                             facecolor="None")
        if idx2 == 2:
            lab = 'T-360GSN'
        else:
            lab = sen
        im = ax[idx].scatter(alon, alat, c=val, s=100,
                             transform=ccrs.Geodetic(), zorder=3, vmin=1.,
                             vmax=5, marker=markers[idx2], label=lab)

        # if idx == 3:
        #    print(ax[idx].get_legend_handles_labels())
        f.close()
        mval = np.mean(val)
        if (len(val) + len(val2)) == 0:
            pro = 'None'
        else:
            pro = 100 * len(val) / (len(val) + len(val2))
        # print(sen + ' ' + str(len(val)+len(val2)))

        print(cfile + ' ' + sen + ' Good sens ' + str(
            len(alat)) + ' mean ' + str(mval) + ' percent ' + str(
            pro) + ' all sens ' + str(len(val) + len(val2)))

# for each map find stations not processed
for idx, cfile in enumerate(files):
    f = open(cfile, 'r')
    snclsgood = []
    for line in f:
        line = line.rstrip()
        line = line.split(',')
        snclsgood.append(line[0])
    f.close()

    inv = client.get_stations(network="II,IC,IU", channel="VHZ",
                              starttime=eves[idx], endtime=eves[idx] + 1,
                              location='00', level='channel')
    sncls = []
    for net in inv:
        for sta in net:
            sncls.append(net.code + '.' + sta.code + '.00.VHZ')

    hatisleft = list(set(sncls) - set(snclsgood))

    wlat, wlon = [], []

    for what in hatisleft:
        coors = inv.get_coordinates(what, eves[idx])
        wlat.append(coors['latitude'])
        wlon.append(coors['longitude'])
    im = ax[idx].scatter(wlon, wlat, c='k', s=100, transform=ccrs.Geodetic(),
                         zorder=3, vmin=1., vmax=5, marker='1')

fig2 = plt.figure(2, figsize=(10, 8))
ax1 = fig2.add_subplot(2, 2, 1, projection=ccrs.Robinson())
sens.append('Earthquake')
markers.append('*')

for idx2, sen in enumerate(sens):

    if idx2 <= 4:
        im3 = ax1.scatter(0, 0, c=4, s=100, transform=ccrs.Geodetic(), zorder=3,
                          vmin=1, vmax=5, marker=markers[idx2], label=sen)
    else:
        im3 = ax1.scatter(0, 0, c='r', s=200, transform=ccrs.Geodetic(),
                          zorder=3, vmin=1., vmax=5, marker=markers[idx2],
                          label=sen)
        im3 = ax1.scatter(0, 0, c='k', s=100, transform=ccrs.Geodetic(),
                          zorder=3, vmin=1., vmax=5, marker='1',
                          label='No Data')
handles, labels = ax1.get_legend_handles_labels()

ax[3].legend(handles, labels, bbox_to_anchor=(1.0, -0.05), ncol=7, fontsize=12)
plt.figure(1)
cb_ax = fig.add_axes([0.2, 0.05, 0.6, 0.03])
cbar = fig.colorbar(im, cax=cb_ax, orientation='horizontal')
plt.close(fig2)

cbar.set_label('Signal-to-Noise Ratio')
plt.savefig('Map_of_usage_ALL.png', bbox_inches='tight', format='PNG', dpi=400)
plt.savefig('Map_of_usage_ALL.pdf', bbox_inches='tight', format='PDF', dpi=400)

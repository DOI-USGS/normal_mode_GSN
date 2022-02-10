#!/usr/bin/env python
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib as mpl
import cartopy as cart
mpl.rc('font', family = 'serif')
mpl.rc('font', serif = 'Times') 
mpl.rc('text', usetex = True)
mpl.rc('font', size=18)

client = Client('IRIS')




fig = plt.figure(1, figsize=(12,12))

stime = UTCDateTime("2014-01-01")
etime = UTCDateTime("2014-01-01")


inv = client.get_stations(network="IU,II,IC", channel="LHZ",
        starttime=stime,endtime=etime, level='response', location='00')

lats, lons, sensors = [], [], []
for net in inv:
    for sta in net:
        for chan in sta:
            lats.append(chan.latitude)
            lons.append(chan.longitude)
            sensors.append(chan.sensor.description)
            #if 'Trillium' in chan.sensor.description:
            #print(sta.code + ' ' + chan.sensor.description)




ax = fig.add_subplot(2,1,1, projection = ccrs.Robinson())
ax.coastlines()
ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
ax.set_global()
ax.coastlines()
sts6, sts1, ks54000,t360, other = 0, 0, 0, 0, 0
sensors1 =[]
for lat, lon, sen in zip(lats, lons, sensors):

    if 'STS-1' in sen:
        ax.plot(lon, lat, label='STS-1', transform=ccrs.Geodetic(), marker="^", markersize=10, zorder=4, color='C0', linestyle='',alpha=0.7,)
        sts1 +=1
        sensors1.append('STS-1')
    elif 'KS' in sen:
        ax.plot(lon, lat, label='KS-54000', transform=ccrs.Geodetic(), marker="^", markersize=10, zorder=4, color='C2', linestyle='',alpha=0.7,)
        ks54000 += 1
        sensors1.append('KS')
    elif 'Trillium 360' in sen:
        ax.plot(lon, lat, label='T-360GSN', transform=ccrs.Geodetic(), marker="^", markersize=10, zorder=4, color='C3', linestyle='',alpha=0.7,)
        t360 +=1
        sensors1.append('Trillium 360')
    elif 'STS-6' in sen:
        ax.plot(lon, lat, label='STS-6', transform=ccrs.Geodetic(), marker="^", markersize=10, zorder=4, color='C1', linestyle='',alpha=0.7,)
        sts6 +=1
        sensors1.append('STS-6')
    else:
        print(sen)
        ax.plot(lon, lat, label='Other', transform=ccrs.Geodetic(), marker="^", markersize=10, zorder=4, color='C4', linestyle='',alpha=0.7,)
        other +=1
        sensors1.append('Other')

#ax.set_title('(a) 2014', loc='left')

lats1, lons1= lats, lons

stime = UTCDateTime("2021-01-01")
etime = UTCDateTime("2021-01-01")


inv = client.get_stations(network="IU,II,IC", channel="LHZ",
        starttime=stime,endtime=etime, level='response', location='00')

lats, lons, sensors = [], [], []
for net in inv:
    for sta in net:
        for chan in sta:
            lats.append(chan.latitude)
            lons.append(chan.longitude)
            sensors.append(chan.sensor.description)
            if 'Trillium' in chan.sensor.description:
                print(sta.code + ' ' + chan.sensor.description)





ax = fig.add_subplot(2,1,2, projection = ccrs.Robinson())
ax.coastlines()
ax.add_feature(cart.feature.LAND, zorder=2, edgecolor='k')
ax.set_global()
ax.coastlines()
sts6, sts1, ks54000,t360, other = 0, 0, 0, 0, 0
for lat, lon, sen in zip(lats, lons, sensors):
    try:
        idx = lons1.index(lon)
        sen1 = sensors1[idx]
    except:
        sen1 = 'Blah'
    if sen1 in sen:
        mark ="^"
    else:
        mark = 'o'
    if 'STS-1' in sen:
        ax.plot(lon, lat, label='STS-1', transform=ccrs.Geodetic(), marker=mark, markersize=10, zorder=4, color='C0', linestyle='',alpha=0.7,)
        sts1 +=1
    elif 'KS' in sen:
        ax.plot(lon, lat, label='KS-54000', transform=ccrs.Geodetic(), marker=mark, markersize=10, zorder=4, color='C2', linestyle='',alpha=0.7,)
        ks54000 += 1
    elif 'Trillium 360' in sen:
        ax.plot(lon, lat, label='T-360GSN', transform=ccrs.Geodetic(), marker=mark, markersize=10, zorder=4, color='C3', linestyle='',alpha=0.7,)
        t360 +=1
    elif 'STS-6' in sen:
        ax.plot(lon, lat, label='STS-6', transform=ccrs.Geodetic(), marker=mark, markersize=10, zorder=4, color='C1', linestyle='',alpha=0.7,)
        sts6 +=1
    else:
        print(sen)
        ax.plot(lon, lat, label='Other', transform=ccrs.Geodetic(), marker=mark, markersize=10, zorder=4, color='C4', linestyle='',alpha=0.7,)
        other +=1
#ax.set_title('(b) 2021', loc='left')
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(0.5, -0.2), loc='lower center', ncol=5, fontsize=18)


print('Here are the STS-6:' + str(sts6))
print('Here are the T-360:' + str(t360))


########################################################################

plt.savefig('figure1.pdf', format='PDF', dpi=500)
plt.savefig('figure1.png', format='PNG', dpi=500)
# print('STS-6s:' + str(sts6))
# print('STS-1s:' + str(sts1))
# print('KS-54000:' + str(ks54000))
# print('Other:' + str(other))



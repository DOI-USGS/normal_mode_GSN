#!/usr/bin/env python
from obspy.core import read, UTCDateTime, Stream
import matplotlib.pyplot as plt
from scipy.signal import periodogram, hilbert
from obspy.signal.invsim import evalresp
import numpy as np
from scipy.optimize import fmin
import math
import glob
import sys
from obspy.clients.fdsn import Client 
from scipy import signal
debug = True

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)



evenum = 0

for evenum in [0,1,2,4]:

    # network
    client = Client("IRIS")
    eves = [UTCDateTime('2021-08-12T18:35:20'), UTCDateTime('2021-07-29T06:15:49'), UTCDateTime('2021-03-04T17:41:25'), 
            UTCDateTime('2013-02-06T23:12:25'), UTCDateTime('2014-04-01T23:46:47')]



    eve = eves[evenum]


    inv = client.get_stations(network="II,IC,IU", channel="LHZ",
            starttime=eve,endtime=eve+1, location='00', level='channel')

    nsls = []
    for net in inv:
        for sta in net:
            for chan in sta:
                print(sta)
                nsls.append([net.code, sta.code, chan.location_code])



    # if evenum == 0:
    #     nsls= [['IU','ADK','00'],['IU','ANMO','00'],['IU','BBSR','00'],['IU','CCM','00'],['IU','CHTO','00'],['IU','DWPF','00'],
    #         ['IU','GNI','00'],['IU','GRFO','00'],['IU','HRV','00'],['IU','KEV','00'], ['IU','KIP','00'],['IU','MBWA','00'],
    #         ['IU','NWAO','00'],['IU','OTAV','00'], ['IU','PAB','00'],['IU','PAYG','00'],['IU','POHA','00'],['IU','RAR','00'],
    #         ['IU','RSSD','00'],['IU','SNZO','00'],['IU','SSPA','00'],['IU','ULN','00'],['IU','WVT','00'],['IU','XMAS','00'],
    #         ['II','ARTI','00'], ['II','BFO','00'], ['II','CMLA','00'], ['II','EFI','00'], ['II','KIV','00'], ['II','KURK','00'], 
    #         ['II','MBAR','00'], ['II','NNA','00'], ['II','OBN','00'], ['II','SIMI','00'], ['II','UOSS','00'], ['II','WRAB','00'],
    #         ['IC','HIA','00'],['IC','MDJ','00']]

    # elif evenum == 1:
    #    nsls = [['IC','HIA','00'],['IC','MDJ','00'],['IC','KMI','00'],['IU','ADK','00'],['IU','ANMO','00'],['IU','COLA','00'],
    #        ['IU','COR','00'],['IU','CTAO','00'],['IU','GNI','00'],['IU','KEV','00'],['IU','KIEV','00'],['IU','KIP','00'],
    #        ['IU','MBWA','00'],['IU','NWAO','00'],['IU','OTAV','00'],['IU','PAB','00'],['IU','PAYG','00'],['IU','SNZO','00'],
    #          ['IU','SSPA','00'],['IU','ULN','00'],['II','ARTI','00'],['II','BFO','00'],['II','CMLA','00'],['II','EFI','00'],
    #         ['II','ESK','00'],['II','KIV','00'],['II','MBAR','00'], ['II','NNA','00'], ['II','OBN','00'], ['II','PALK','00'],
    #            ['II','PFO','00'], ['II','SIMI','00'], ['II','SUR','00'], ['II','TLY','00'], ['II','UOSS','00'],['II','WRAB','00']]

    # elif evenum == 2:
    #     nsls = [['II','ALE','00'], ['II','BFO','00'], ['II','CMLA','00'],['II','JTS','00'],['II','KURK','00'],['II','MBAR','00'],
    #         ['II','NNA','00'],['II','PALK','00'],['II','PFO','00'],['II','RPN','00'],['II','SUR','00'],['II','UOSS','00'],
    #         ['IC','KMI','00'],['IC','LSA','00'],['IC','MDJ','00'],['IC','WMQ','00'],['IC','XAN','00'],['IU','ANMO','00'],
    #         ['IU','CCM','00'],['IU','CHTO','00'],['IU','COLA','00'],['IU','COR','00'],['IU','CTAO','00'],
    #         ['IU','DWPF','00'],['IU','GNI','00'],['IU','HRV','00'],['IU','KIEV','00'],['IU','KIP','00'],['IU','MAJO','00'],
    #         ['IU','MBWA','00'],['IU','NWAO','00'],['IU','PAB','00'],['IU','PAYG','00'],['IU','POHA','00'],['IU','RAR','00'],
    #         ['IU','RSSD','00'],['IU','SNZO','00'],['IU','TUC','00'],['IU','ULN','00'],['IU','WCI','00'],['IU','WVT','00'],['IU','XMAS','00']]
    # elif evenum == 3:
    #     nsls = [['II','NNA','00'], ['II','KURK','00'], ['II','ESK','00'], ['IC','KMI','00'], ['IU','KIP','00'], ['IU','LVC','00'], ['IU','WCI','00']]
    # elif evenum == 4:
    #     nsls = [['IU','AFI','00'],['IU','COR','00'] ,
    #         ['IU','CTAO','00'],['IU','GNI','00'],['IU','KIP','00'],['IU','PAB','00'],['IU','TUC','00'],['IU','ULN','00'],
    #         ['IC','KMI','00'],['IC','MDJ','00'],['IC','XAN','00'],['IC','WMQ','00'],['II','AAK','00'],['II','BFO','00'],
    #         ['II','BRVK','00'],['II','KURK','00'],['II','NNA','00'],['II','RPN','00'],['IU','TLY','00']]



    # min and max frequency in milliHz
    mf, Mf= 0.2, 1.0


    modeidx = 5

    # In days
    length = 400

    # Hours of data you want analyzed for sub 1 mHz you want to get 3 or more days
    #hours = int((1000.*1.1*195./0.76593)/(60.*60.))
    #hours = 126

    # List of reference modes to plot (only good for sub 1 mHz)
    modes = [[0.8140, '0S0', 5340.67285, 0.79], [0.3092, '0S2', 509.68240, 0.31], [0.46884, '0S3', 417.5490, 0.45],
                [0.6469, '0S4', 373.279, 0.6], [0.67934, '1S2', 310.27490, 0.65], [0.84022, '0S5', 355.695, 0.78]]



    for modeidx in [1,2,3,5]:

        sensors = {}



        hours = int((1000*1.1*modes[modeidx][2]/modes[modeidx][0])/(60.*60.))
        print(hours)

        numsts1, numsts6, num360, numks = 0, 0, 0, 0

        st = Stream()
        for nsl in nsls:
            if 'inv' not in vars():
                inv = client.get_stations(network=nsl[0], station=nsl[1], starttime=eve, 
                                  endtime = eve + hours*60*60, level="response", location='00', channel='VHZ')
            else:
                try:
                    inv += client.get_stations(network=nsl[0], station=nsl[1], starttime=eve, 
                                  endtime = eve + hours*60*60, level="response", location='00', channel='VHZ')
                except:
                    continue
            try:
                st += client.get_waveforms(nsl[0], nsl[1], nsl[2], 'VHZ', eve, eve + hours*60*60)
            except:
                continue

            for net in inv:
                for sta in net:
                    for chan in sta:
                        sensors[nsl[0] + '.' + nsl[1] + '.00.VHZ'] = chan.sensor.description


        for key in sensors:
            sen = sensors[key]
            if 'STS-1' in sen:
                numsts1 += 1
            elif 'STS-6' in sen:
                numsts6 += 1
            elif '360' in sen:
                print(sen)
                num360 += 1
            elif 'KS' in sen:
                numks += 1

        print('STS-6:' + str(numsts6))
        print('STS-1:' + str(numsts1))
        print('T-360:' + str(num360))
        print('KS-5400:' + str(numks))


        st.merge(fill_value=0)

        st.detrend('linear')
        st.detrend('constant')
        st.sort(['station'])

        fig, ax = plt.subplots(len(st), 1, figsize=(16,12))
        ax = ax.flatten()


        # Window
        for tr in st:
            tr.data *= signal.get_window(('kaiser', 2.*np.pi), tr.stats.npts)
        # Zero pad
        for tr in st:
            tr.data = np.pad(tr.data, (int((length*400*24/10. -tr.stats.npts)/2.),int((length*400*24/10. -tr.stats.npts)/2.)), 'edge') 
            print(tr)


        NFFT=2**(math.ceil(math.log(st[0].stats.npts, 2)))

        fhand = open('SNRresults_' +str(eve.year)+ '_' + str(eve.julday) + '_' + modes[modeidx][1] + '.csv','w')

        for idx2, tr in enumerate(st):
            
            
            if True:

                f,p= periodogram(tr.data, fs=tr.stats.sampling_rate, nfft= NFFT, scaling='spectrum')
                p, f = p[1:], f[1:]
                trid = tr.id
                inv_resp = inv.get_response(trid, tr.stats.starttime)
                resp, _ = inv_resp.get_evalresp_response(tr.stats.delta, NFFT, 'ACC')
                resp = resp[1:]
                # Convert units to nm/s/s
                p = np.sqrt(p/(np.abs(resp)**2))*10**9
                # Now have p in nm/s/s switch f to mHz
                f *= 1000.
                p = p[(f >= mf) & (f <= Mf)]
                f = f[(f >= mf) & (f <= Mf)]
                sig = np.max(p[(f >= modes[modeidx][0]-0.01) & (f <= modes[modeidx][0]+0.01)])
                #sig = np.max(p[(f >= 0.83) & (f <= 0.85)])
                noise = np.max(p[(f >= modes[modeidx][3]-0.01) & (f <= modes[modeidx][3]+0.01)])
                noisem = np.mean(p[(f >= modes[modeidx][3]-0.01) & (f <= modes[modeidx][3]+0.01)])

                snr = sig/noisem 
                coors = inv.get_coordinates(tr.id, tr.stats.starttime)

                

                p /= np.max(p)

                if 'STS-1' in sensors[tr.id]:
                    c='C0'
                    lab = 'STS-1'
                elif 'STS-6' in sensors[tr.id]:
                    c='C1'
                    lab ='STS-6'
                elif 'KS' in sensors[tr.id]:
                    c = 'C2'
                    lab = 'KS-54000'
                elif '360' in sensors[tr.id]:
                    c = 'C3'
                    lab = 'Trillium 360'
                else:
                    c='C4'
                    lab = 'Other'

                fhand.write(tr.id + ', ' + str(snr) + ',' + str(coors['latitude']) + ',' + str(coors['longitude']) + ',' + str(sig) + ',' + str(noisem) + ',' + lab + '\n')


                ax[idx2].plot(f,p, label=(tr.id).replace('.',' '), alpha=0.7, color=c)
                ax[idx2].fill_between(f, 0, p, alpha=0.7, color=c)
                #ax[idx2].legend(loc='upper center', fontsize=12, ncol=5)
                ax[idx2].set_xlim((mf,Mf))
                _, top = ax[idx2].get_ylim()

                

                #ax[idx2].set_ylabel('Amplitude ($nm/s^2$)')

            
                locp = top*.8
                
                for nidx, mode in enumerate(modes):
                    if (mode[0] <= max(f)) and (mode[0] >= min(f)):
                        label = mode[1]
                        label = '_' + label[0] + label[1] + '_' + label[2]
                        if len(label) == 4:
                            label += label[3]
                        ax[idx2].plot((mode[0], mode[0]), (-10., 500.), color='grey', alpha=0.5)
                        if idx2 == 0:
                            ax[idx2].text(mode[0] + .003, 1.3, '$' + label + '$' , fontsize=18)
                        locp -= top*0.15
                        if locp < 0.45*top:
                            locp = .8*top
                #, [0.93908, '1S3'], 
                #[0.94434, '3S1']
                ax[idx2].plot((0.93908,0.93908), (-10,500.), color='grey', alpha=0.5)
                if idx2 == 0:
                    ax[idx2].text(0.93908+ .003, 1.3, '$_1S_3$/$_3S_1$', fontsize=18)
                ax[idx2].set_ylim((0,1))
            if idx2 < len(st)-1:
                ax[idx2].set_xticks([])
            
            ax[idx2].set_yticks([0.5])
            ax[idx2].set_yticklabels([tr.stats.station])
        ax[idx2].set_xlabel('Frequency ($mHz$)')
        plt.savefig('BIGONE_' + str(eve.julday).zfill(3) + '_' + modes[modeidx][1]  + '.jpg', format='JPEG', dpi=400)
        plt.close('all') 
        fhand.close() 
            #plt.savefig('Test' + '.pdf', format='PDF', dpi=400) 
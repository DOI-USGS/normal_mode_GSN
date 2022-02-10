#!/usr/bin/env python
import glob


# # List of reference modes to plot (only good for sub 1 mHz)
# modes = [[0.8140, '0S0', 5340.67285, 0.79], [0.3092, '0S2', 509.68240, 0.31], [0.46884, '0S3', 417.5490, 0.45],
#         [0.6469, '0S4', 373.279, 0.6], [0.67934, '1S2', 310.27490, 0.65], [0.84022, '0S5', 355.695, 0.78]]

files = ['SNRresults_2014_91.csv', 'SNRresults_2021_63.csv','SNRresults_2021_210.csv','SNRresults_2021_224.csv']



for cfile in files:
    # Now we open each file for that event
    evefiles = glob.glob(cfile.replace('.csv','_*'))
    lat, lon, sig, noise, sta, sensor = [], [], [], [], [], []
    for evefile in evefiles:
        f = open(evefile, 'r')
        for line in f:
            line = (line.rstrip()).split(',')
            lat.append(float(line[2]))
            lon.append(float(line[3]))
            sig.append(float(line[4]))
            noise.append(float(line[5]))
            sta.append(line[0])
            sensor.append(line[-1])
        f.close()
    asta = list(set(sta))
    asnr, asensor, alat, alon = [], [], [], []
    f = open('ALL_' + cfile,'w')
    for csta in asta:
        print(csta)
        idx = [index for index, element in enumerate(sta) if element == csta]
        anoise, asig = 0, 0
        for cidx in idx:
            print(cidx)
            print(noise[cidx])
            anoise += float(noise[cidx])
            asig += float(sig[cidx])
        f.write(csta + ',' + str(lat[cidx]) + ',' + str(lon[cidx]) + ',' + str(asig/anoise) + ',' + sensor[cidx] + '\n')
    f.close()

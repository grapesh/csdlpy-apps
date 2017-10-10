# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 20:13:19 2017

@author: svinogra
"""
import argparse
import sys
import glob
from datetime import datetime as dt
import numpy as np
import shutil

import matplotlib
matplotlib.use('Agg',warn=False)
import matplotlib.pyplot as plt

#==============================================================================
def read_metrics_cfg (cfgFile):
    params = dict()
    f = open(cfgFile,'r')
    params['domain']       = f.readline().split(',')[0]
    params['modelPath']    = f.readline().split(',')[0]
    params['fileMask']     = f.readline().split(',')[0]
    
    params['stats']        = f.readline().split(',')
    params['stats']        = params['stats'][:-1]
    params['stats']        = map(str.strip, params['stats'])
    
    params['units']        = f.readline().split(',')
    params['units']        = params['units'][:-1]
    params['units']        = map(str.strip, params['units'])
    
    params['outputPath']   = f.readline().split(',')[0]
    params['archivePath']  = f.readline().split(',')[0]
    
    f.close()
    return params

#==============================================================================
def run_metrics (argv):

    toolkitPath = ''
    cfgFile     = ''
    cycle       = ''

    parser = argparse.ArgumentParser()
    parser.add_argument('-t','--toolkitPath', required=True)
    parser.add_argument('-c','--cfgFile',     required=True)
    parser.add_argument('-f','--fctCycle',    required=False)
    args = parser.parse_args()
    if args.toolkitPath:
        toolkitPath = args.toolkitPath
    if args.cfgFile:
        cfgFile     = args.cfgFile
    if args.fctCycle:
        cycle    = args.fctCycle

    print 'event.py configured with :'
    print 'cfgFile     =', cfgFile
    print 'fctCycle    =', cycle

    sys.path.insert(0, toolkitPath )
    import csdlpy

    if cycle == '':
        latest = csdlpy.estofs.latestForecast ()
    else:        
        latest               = dict()
        latest['yyyymmdd']   = cycle[0:8]
        latest['tHHz'] = 't' + cycle[8:10] + 'z'
    print '[info]: requesting ', latest

    params   = read_metrics_cfg (cfgFile)
    compute_metrics(params, toolkitPath, latest)

#==============================================================================
def compute_metrics(params, toolkitPath, latest):

    sys.path.insert(0, toolkitPath )
    import csdlpy
    
    stats = params['stats']
    units = params['units']
    YYYYMMDD = latest['yyyymmdd']
    ncFiles  = glob.glob(params['modelPath'] + params['fileMask'])
    
    datespan = [dt.strptime(YYYYMMDD,"%Y%m%d"), dt.strptime(YYYYMMDD,"%Y%m%d")]

    forecasts = list()
    for ncFile in ncFiles:
        forecasts.append (csdlpy.estofs.getPointsWaterlevel (ncFile) )

    #Get all stations IDs from the first forecast
    station_id  = []
    station_lon = []
    station_lat = []
    for n in range(len(forecasts[0]['stations'])):
        station_id.append (forecasts[0]['stations'][n].split()[2])
        station_lon.append(forecasts[0]['lon'][n])
        station_lat.append(forecasts[0]['lat'][n])

    #Set forecast lead times
    fxLead = np.array(range(-144,12,6),dtype='int')

    statvals = np.empty( (len(fxLead), len(stats), len(station_id)) ,dtype='float')
    statavgs = np.empty( (len(fxLead), len(stats)), dtype='float')

    #Set output files
    fout = list()
    for s in stats:
        fout.append(open('metrics-'+ YYYYMMDD +'-'+ s +'.csv','w',0))
        plt.figure(figsize=(10,4.5))

    # Collect figure handles
    #figs = list(map(plt.figure, plt.get_fignums()))

    #Write the header
    header = 'coops_id, lon, lat, ' + ','.join(str(e) for e in fxLead) + ','
    for f in fout:
        f.write(header + '\n')

    #Write rows station by station
    for n in range(len(station_id)):

        idrow  = station_id[n] +','+ str(station_lon[n]) +','+ str(station_lat[n])
        for f in fout:
            f.write(idrow)

        mx, fx   = csdlpy.estofs.metrics(station_id[n],  datespan, forecasts)
        fx       = np.array(np.round(fx),dtype='int')

        n_fx = 0
        for f in fxLead:
            idx = np.where(fx==f)
            n_st = 0
            if np.array(idx).size:
                for s in stats:
                    statvals[n_fx, n_st, n] = mx[idx[0][0]][s]
                    n_st += 1
            else:
                statvals[n_fx, :] = np.nan
            
            for n_st in range(len(stats)):
                fout[n_st].write(',' + str(statvals[n_fx, n_st, n]))

            n_fx += 1
        for f in fout:
            f.write('\n')

        for n_st in range(len(stats)):
            plt.figure(n_st)
            plt.plot(fxLead, statvals[:,n_st, n],c='lightgray')

    #Compute averages
    for n_st in range(len(stats)):
        fout[n_st].write('Mean,over,stations')
        for n_fx in range(len(fxLead)):
            statavgs[n_fx, n_st] = np.nanmean(statvals[n_fx, n_st, :])
            fout[n_st].write(',' + str(statavgs[n_fx, n_st]))

        plt.figure(n_st)
        plt.title (YYYYMMDD + ' ' + stats[n_st] + ' (' + units[n_st] + ')')

        if stats[n_st]=='vexp':
            plt.ylim([0,100])
        if stats[n_st]=='plag':
            plt.ylim([-720,720])
        if stats[n_st]=='npts':
            plt.ylim([0, 240])
        if stats[n_st]=='rmsd':
            plt.ylim([0, 0.6])
        if stats[n_st]=='rval':
            plt.ylim([0., 1.])
        if stats[n_st]=='skil':
            plt.ylim([0, 1.])
        if stats[n_st]=='peak':
            plt.ylim([-.5, .5])
        if stats[n_st]=='bias':
            plt.ylim([-.5, .5])
        
        
        plt.plot(fxLead, statavgs[:,n_st],c='k',linewidth=3.0)
        plt.xlabel('LEAD TIME, HOURS')
        plt.grid()
        plt.savefig('ts-' + YYYYMMDD + '-' + stats[n_st] + '.png')

    for f in fout:
        f.write('\n')
        f.close()

#Copy 
    for s in stats:
        shutil.copy('metrics-'+ YYYYMMDD +'-'+ s +'.csv', 'metrics-'+s+'.csv')
        shutil.copy('ts-'+ YYYYMMDD +'-'+ s +'.png', 'ts-'+s+'.png')
    
#==============================================================================
if __name__ == "__main__":

    run_metrics (sys.argv[1:])



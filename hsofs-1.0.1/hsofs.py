# -*- coding: utf-8 -*-
"""
@author: Sergey.Vinogradov
"""
import sys
import argparse
import matplotlib
matplotlib.use('Agg',warn=False)
import matplotlib.pyplot as plt
from datetime import datetime as dt
from datetime import timedelta as td
import numpy as np
import copy

#==============================================================================
def run_hsofs (argv):
   
    inputPath   = ''
    stormID     = ''
    inputCycle  = ''
    outputPath  = ''
    toolkitPath = ''
    cfgFile     = ''

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-i','--inputPath',  required=True)
    parser.add_argument('-s','--stormID',    required=True)
    parser.add_argument('-z','--inputCycle',  required=True)    
    parser.add_argument('-o','--outputPath',  required=True)
    parser.add_argument('-t','--toolkitPath', required=True)
    parser.add_argument('-c','--cfgFile',     required=True)
    args = parser.parse_args()
    if args.inputPath:
        inputPath = args.inputPath
    if args.stormID:
        stormID = args.stormID
    if args.inputCycle:
        inputCycle = args.inputCycle
    if args.outputPath:
        outputPath = args.outputPath
    if args.toolkitPath:
        toolkitPath = args.toolkitPath
    if args.cfgFile:
        cfgFile     = args.cfgFile
        
    print 'hsofs.py configured with :'
    print 'inputPath   =', inputPath
    print 'stormID     =', stormID
    print 'outputCycle =', inputCycle
    print 'outputPath  =', outputPath
    print 'toolkitPath =', toolkitPath
    print 'cfgFile     =', cfgFile
    
    params = read_event_cfg (cfgFile)
    hsofs_plots (params, inputPath, stormID, inputCycle, outputPath, toolkitPath)
    
#==============================================================================
#Harvey,                                        # Event Name
#-100., -80., 20., 30.,                         # lonlim, latlim 
#0., 4.,                                        # clim
#http://ftp.nhc.noaa.gov/atcf/btk/bal092017.dat,# Best Track
#http://ftp.nhc.noaa.gov/atcf/fst/al092017.fst, # Fcst
#estofs_atl,                                    # Domain 
#estofs.atl,                                    # Prefix      
#/gpfs/hps/nco/ops/com/estofs/prod/,            # Path
#ftp://ocsftp.ncd.noaa.gov/estofs/atl/,         # Grid 
#-98., -89., 26., 30.,                          # bbox for statons output requested:

#TODO: cfg API across CSDLPY 
    
def read_event_cfg (cfgFile):
    params = dict() 
    f = open(cfgFile,'r')
    params['eventName'] = f.readline().split(',')[0]
    line = f.readline().split(',')
    params['lonlim'] = float(line[0]), float(line[1])
    params['latlim'] = float(line[2]), float(line[3])
    line = f.readline().split(',')
    params['clim'] = (float(line[0]), float(line[1]))
    params['bestTrackURL'] = f.readline().split(',')[0]
    params['advTrackURL']  = f.readline().split(',')[0]
    params['domain']       = f.readline().split(',')[0]
    params['prefix']       = f.readline().split(',')[0]
    params['prodPath']     = f.readline().split(',')[0]
    params['gridPath']     = f.readline().split(',')[0]
    line = f.readline().split(',')
    params['xlim'] = float(line[0]), float(line[1])
    params['ylim'] = float(line[2]), float(line[3])    
    f.close()
    return params

#==============================================================================
def hsofs_plots (params, inputPath, stormID, inputCycle, outputPath, toolkitPath):

    sys.path.insert(0, toolkitPath )
    from csdlpy import estofs
    from csdlpy import adcirc
    from csdlpy import transfer
    from csdlpy import plotter
    from csdlpy import atcf
    from csdlpy.obs import coops
    
    ens = "nhctrk", "higherSpeed", "lowerSpeed","shiftLeft","shiftRight "
    ens_col = "k", "r","b","c","m"
    
    maxPSfile = 'inputPath' + '/hsofs.' + stormID + '.' + inputCycle + '.fields.maxPS.nc'
    maxPS = estofs.getFieldsWaterlevel (maxPSfile, 'zeta_max')
    
    fcst = dict()
    counter = 0
    for member in ens:
        print '[info]: working on ' + member
        trackFile  = 'inputPath' + '/hsofs.' + stormID + '.' + inputCycle + '.' + member + '.surfaceforcing' 
        maxeleFile = 'inputPath' + '/hsofs.' + stormID + '.' + inputCycle + '.' + member + '.fields.maxele.nc'
        pointsFile = 'inputPath' + '/hsofs.' + stormID + '.' + inputCycle + '.' + member + '.points.waterlevel.nc'
        
        fcst['ens'][counter] = member
        fcst['track'][counter]  = atcf.readTrack(trackFile)
        fcst['maxele'][counter] = estofs.getFieldsWaterlevel (maxeleFile, 'zeta_max')
        fcst['points'][counter] = estofs.getPointsWaterlevel (pointsFile )            
        counter += 1
    
    gridFile = 'fort.14'
    trkFile = 'trk.dat'
    advFile = 'adv.dat'

    transfer.download (params['gridPath'], gridFile)
    transfer.download (params['bestTrackURL'], trkFile)
    transfer.download (params['advTrackURL'],  advFile)    

    grid = adcirc.readGrid (gridFile)
    trk = atcf.readTrack(trkFile)
    adv = atcf.readTrack(advFile)
    
    cwl = fcst['points'][0] # nhctrk stations
    
    mod_dates = cwl['time'] 
    #          
    # Plot maxeles maxPS
    f = plotter.plotMap    (params['lonlim'], params['latlim'], fig_w=10.)
    plotter.addSurface (grid, maxPS['value'],clim=params['clim'])
    
    plt.plot(trk['lon'], trk['lat'],'o-k',markersize=1,zorder=10)    
    plt.plot(adv['lon'], adv['lat'],'o--r',markersize=1,zorder=11)    
        
    title = 'HSOFS experimental ' + stormID + '.' + inputCycle + '.maxPS'
    plt.text (params['lonlim'][0]+0.03, \
              params['latlim'][0]+0.03, \
                  title )    
    plotter.save(title, outputPath + '/' + stormID + '.' + inputCycle + '.maxPS.maxele.png')
    plt.close(f)
    #

    #          
    # Plot ens maxeles
    counter = 0
    for member in ens:
        
        f = plotter.plotMap    (params['lonlim'], params['latlim'], fig_w=10.)
        plotter.addSurface (grid, fcst['maxele'][counter]['value'],clim=params['clim'])
        
        plt.plot(trk['lon'], trk['lat'],'o-k',markersize=1,zorder=10)    
        plt.plot(adv['lon'], adv['lat'],'o--r',markersize=1,zorder=11)    
        plt.plot(fcst['track'][counter]['lon'], fcst['track'][counter]['lat'],'o-k',markersize=2,zorder=12)
        
        counter += 1
        
        # Add stations locations:                      
        for n in range(len(cwl['lon'])):
            if params['xlim'][0] <= cwl['lon'][n] and \
                     cwl['lon'][n] <= params['xlim'][1] and \
                     params['ylim'][0] <= cwl['lat'][n] and \
                     cwl['lat'][n] <= params['ylim'][1]:
                         
               plt.plot(cwl['lon'][n], cwl['lat'][n],'wo', \
                        markeredgecolor='k', zorder=15)
        title = 'HSOFS experimental ' + stormID + '.' + inputCycle + '.' + member
        plt.text (params['lonlim'][0]+0.03, \
                  params['latlim'][0]+0.03, \
                        title )    
        plotter.save(title, outputPath + '/' + stormID + '.' + inputCycle + '.' + member + '.maxele.png')
        plt.close(f)
    #
    
    # Plot time series
    utcnow = dt.utcnow()
    daterange = (utcnow-td(days=1), utcnow)
        
    figureCounter = 0
    for n in range(len(cwl['lon'])):
        
        if params['xlim'][0] <= cwl['lon'][n] and \
           cwl['lon'][n] <= params['xlim'][1] and \
           params['ylim'][0] <= cwl['lat'][n] and \
           cwl['lat'][n] <= params['ylim'][1]:
               
               mod_vals = np.squeeze( cwl['zeta'][:,n] )
    
               obs_vals = copy.deepcopy(mod_vals) #in case there is no obs
               obs_dates = copy.deepcopy(mod_dates) # in case if there is no obs
               
               coops_id = cwl['stations'][n].split()[2]
               print coops_id
               obs_coops = coops.getData ( coops_id, daterange)
               if len(obs_coops['values'])>0:
                   obs_dates = obs_coops['dates']
                   obs_vals  = obs_coops['values']
                                        
               # Plot observations
               f = plt.figure(figsize=(20,4.5))               
               plt.plot(obs_dates, obs_vals, '.',color='g',label='OBS')
               
               # Plot model
               counter = 0
               for member in ens:
                   mod_vals = fcst['points'][counter]['zeta'][:,n]
                   plt.plot(mod_dates, mod_vals, '.', color=ens_col[counter],label=member)
                   counter += 1 
                   
               plt.ylim([-1.5,params['clim'][1]])
               plt.title(cwl['stations'][n])
               plt.xlim([daterange[0],mod_dates[-1]])
               plt.grid()
               plt.xlabel('DATE UTC')
               plt.ylabel('WL, meters LMSL')
               plt.legend(bbox_to_anchor=(0.9, 0.35))

               figFile=outputPath + '/ts-' + str(figureCounter).zfill(3) + '.png'
               plt.savefig(figFile)
               plt.close(f)
               figureCounter += 1
               

#==============================================================================
if __name__ == "__main__":

    run_hsofs (sys.argv[1:])

    


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
def run_event (argv):
   
    outputPath  = ''
    toolkitPath = ''
    cfgFile     = ''

    parser = argparse.ArgumentParser()
    parser.add_argument('-o','--outputPath',  required=False)
    parser.add_argument('-t','--toolkitPath', required=True)
    parser.add_argument('-c','--cfgFile',     required=True)
    args = parser.parse_args()
    if args.outputPath:
        outputPath = args.outputPath
    if args.toolkitPath:
        toolkitPath = args.toolkitPath
    if args.cfgFile:
        cfgFile     = args.cfgFile
        
    print 'event.py configured with :'
    print 'outputPath  =', outputPath
    print 'toolkitPath =', toolkitPath
    print 'cfgFile     =', cfgFile
    
    params = read_event_cfg (cfgFile)
    event_maxele (params, outputPath, toolkitPath)
    
#==============================================================================
#Henry,                                         # Event Name
#-100., -80., 20., 30.,                         # lonlim, latlim 
#0., 4.,                                        # clim
#http://ftp.nhc.noaa.gov/atcf/btk/bal092017.dat,# Best Track
#http://ftp.nhc.noaa.gov/atcf/fst/al092017.fst, # Fcst
#estofs_atl,                                    # Domain 
#estofs.atl,                                    # Prefix      
#/gpfs/hps/nco/ops/com/estofs/prod/,            # Path
#ftp://ocsftp.ncd.noaa.gov/estofs/atl/,         # Grid 
#-98., -89., 26., 30.,                          # bbox for statons output requested:

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
def event_maxele (params, outputPath, toolkitPath):

    sys.path.insert(0, toolkitPath )
    from csdlpy import estofs
    from csdlpy import adcirc
    from csdlpy import transfer
    from csdlpy import plotter
    from csdlpy import atcf
    from csdlpy.obs import coops
    
    latest = estofs.latestForecast ()
    
    maxeleFile = params['prodPath'] + params['domain'] + '.' + \
                       latest['yyyymmdd'] + '/' + params['prefix'] + '.' + \
                       latest['tHHz'] + '.fields.cwl.maxele.nc'
                         
    maxele = estofs.getFieldsWaterlevel (maxeleFile, 'zeta_max')    
    
    gridFile = 'fort.14'
    trkFile = 'trk.dat'
    advFile = 'adv.dat'

    transfer.download (params['gridPath'], gridFile)
    transfer.download (params['bestTrackURL'], trkFile)
    transfer.download (params['advTrackURL'],  advFile)    

    grid = adcirc.readGrid (gridFile)
    trk = atcf.readTrack(trkFile)
    adv = atcf.readTrack(advFile)
    
    plotter.plotMap    (params['lonlim'], params['latlim'], fig_w=10.)   
    try:
        plotter.addSurface (grid, maxele['value'],clim=params['clim'])
    except:
        print '[warn]: cannot plot surface'
    plt.plot(trk['lon'], trk['lat'],'o-k',markersize=2,zorder=10)    
    plt.plot(adv['lon'], adv['lat'],'o-r',markersize=2,zorder=11)    
    
    cwlFile = params['prodPath'] + params['domain'] + '.' + \
                       latest['yyyymmdd'] + '/' + params['prefix'] + '.' + \
                       latest['tHHz'] + '.points.cwl.nc'
    htpFile = params['prodPath'] + params['domain'] + '.' + \
                       latest['yyyymmdd'] + '/' + params['prefix'] + '.' + \
                       latest['tHHz'] + '.points.htp.nc'
                       
    cwl = estofs.getPointsWaterlevel (cwlFile)
    htp = estofs.getPointsWaterlevel (htpFile)
    mod_dates = cwl['time']
    htp_dates = htp['time']
    
    for n in range(len(cwl['lon'])):        
           
        if params['xlim'][0] <= cwl['lon'][n] and \
           cwl['lon'][n] <= params['xlim'][1] and \
           params['ylim'][0] <= cwl['lat'][n] and \
           cwl['lat'][n] <= params['ylim'][1]:
               plt.plot(cwl['lon'][n], cwl['lat'][n],'wo', \
                        markeredgecolor='k', zorder=15)
    title = 'ESTOFS (GFS) ' + latest['yyyymmdd'] + '.' + latest['tHHz']
    plt.text (params['lonlim'][0]+0.02, \
              params['latlim'][0]+0.02, \
              title )    
    plotter.save(title, outputPath + '/maxele.png')        

    utcnow = dt.utcnow()
    daterange = (utcnow-td(days=2), utcnow)
        
    figureCounter = 0
    for n in range(len(cwl['lon'])):
        
        if params['xlim'][0] <= cwl['lon'][n] and \
           cwl['lon'][n] <= params['xlim'][1] and \
           params['ylim'][0] <= cwl['lat'][n] and \
           cwl['lat'][n] <= params['ylim'][1]:
               
               mod_vals = np.squeeze( cwl['zeta'][:,n] )
               htp_vals = np.squeeze( htp['zeta'][:,n] )
    
               obs_vals = copy.deepcopy(mod_vals) #in case there is no obs
               obs_dates = copy.deepcopy(mod_dates) # in case if there is no obs
               
               coops_id = ''               
               keepSearch = True
               maybe_ids = cwl['stations'][n].split()               
               # Find which one is coops_id               
               for ids in maybe_ids:
                  if keepSearch:
                      try:
                         int_id = int(ids)
                         if len(str(int_id)) == 7:
                             coops_id = str(int_id)               
                             keepSearch = False
                      except:
                          pass
                  
               print coops_id
               obs_coops = coops.getData ( coops_id, daterange)
               if len(obs_coops['values'])>0:
                   obs_dates = obs_coops['dates']
                   obs_vals  = obs_coops['values']
                                        
               f = plt.figure()               
               figureCounter += 1
               plotter.plot_estofs_timeseries(obs_dates, obs_vals, \
                                              mod_dates, mod_vals, \
                                              figFile=outputPath + '/ts-' + \
                                              str(figureCounter).zfill(3) + '.png', \
                                              stationName=latest['yyyymmdd']+'.t' +latest['tHHz']+'z '+ cwl['stations'][n], \
                                              htp_dates=htp_dates, htp_vals=htp_vals, \
                                              daterange = daterange, \
                                              ylim = (-1.5,params['clim'][1]))
               plt.close(f)
    

#==============================================================================
if __name__ == "__main__":

    run_event (sys.argv[1:])

    


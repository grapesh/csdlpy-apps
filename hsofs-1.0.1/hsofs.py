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

    gridFile = 'fort.14'
    trkFile = 'trk.dat'
    advFile = 'adv.dat'

    transfer.download (params['gridPath'], gridFile)
    transfer.download (params['bestTrackURL'], trkFile)
    transfer.download (params['advTrackURL'],  advFile)    

    grid = adcirc.readGrid (gridFile)
    trk = atcf.readTrack(trkFile)
    adv = atcf.readTrack(advFile)
  
    # nhctrk
    ens = 'nhctrk'
    ncfile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.' + ens + '.fields.maxele.nc'
    trkFile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.'+ ens + '.surfaceforcing' 
    pointsFile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.'+ ens + '.points.waterlevel.nc'
    
    maxele = estofs.getFieldsWaterlevel (ncfile, 'zeta_max')    
    track_nhctrk  = atcf.readTrack(trkFile)
    points_nhctrk = estofs.getPointsWaterlevel (pointsFile )

#    # Plot maxeles maxele_nhctrk
#    f = plotter.plotMap    (params['lonlim'], params['latlim'], fig_w=10.)
#    plotter.addSurface (grid, maxele['value'],clim=params['clim'])
#    plt.plot(track_nhctrk['lon'], track_nhctrk['lat'],'o-k',markersize=1,zorder=10)    
#    title = 'HSOFS experimental ' + stormID + '.' + inputCycle + '.' + ens
#    plt.text (params['lonlim'][0]+0.03, \
#              params['latlim'][0]+0.03, \
#                  title )    
#    plotter.save(title, outputPath + '/' + stormID + '.' + inputCycle + '.'+ ens + '.maxele.png')
#    plt.close(f)
    
    # 2
    ens = 'lowerSpeed'
    ncfile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.' + ens + '.fields.maxele.nc'
    trkFile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.'+ ens + '.surfaceforcing' 
    pointsFile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.'+ ens + '.points.waterlevel.nc'
    
    maxele = estofs.getFieldsWaterlevel (ncfile, 'zeta_max')    
    track_lowerSpeed  = atcf.readTrack(trkFile)
    points_lowerSpeed = estofs.getPointsWaterlevel (pointsFile )

#    # Plot maxeles maxele_nhctrk
#    f = plotter.plotMap    (params['lonlim'], params['latlim'], fig_w=10.)
#    plotter.addSurface (grid, maxele['value'],clim=params['clim'])
#    plt.plot(track_lowerSpeed['lon'], track_lowerSpeed['lat'],'o-k',markersize=1,zorder=10)    
#    title = 'HSOFS experimental ' + stormID + '.' + inputCycle + '.' + ens
#    plt.text (params['lonlim'][0]+0.03, \
#              params['latlim'][0]+0.03, \
#                  title )    
#    plotter.save(title, outputPath + '/' + stormID + '.' + inputCycle + '.'+ ens + '.maxele.png')
#    plt.close(f)
    
    # 3
    ens = 'shiftRight'
    ncfile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.' + ens + '.fields.maxele.nc'
    trkFile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.'+ ens + '.surfaceforcing' 
    pointsFile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.'+ ens + '.points.waterlevel.nc'
    
    maxele = estofs.getFieldsWaterlevel (ncfile, 'zeta_max')    
    track_shiftRight  = atcf.readTrack(trkFile)
    points_shiftRight = estofs.getPointsWaterlevel (pointsFile )

#    # Plot maxeles maxele_nhctrk
#    f = plotter.plotMap    (params['lonlim'], params['latlim'], fig_w=10.)
#    plotter.addSurface (grid, maxele['value'],clim=params['clim'])
#    plt.plot(track_shiftRight['lon'], track_shiftRight['lat'],'o-k',markersize=1,zorder=10)    
#    title = 'HSOFS experimental ' + stormID + '.' + inputCycle + '.' + ens
#    plt.text (params['lonlim'][0]+0.03, \
#              params['latlim'][0]+0.03, \
#                  title )    
#    plotter.save(title, outputPath + '/' + stormID + '.' + inputCycle + '.'+ ens + '.maxele.png')
#    plt.close(f)

    # 4
    ens = 'shiftLeft'
    ncfile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.' + ens + '.fields.maxele.nc'
    trkFile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.'+ ens + '.surfaceforcing' 
    pointsFile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.'+ ens + '.points.waterlevel.nc'
    
    maxele = estofs.getFieldsWaterlevel (ncfile, 'zeta_max')    
    track_shiftLeft  = atcf.readTrack(trkFile)
    points_shiftLeft = estofs.getPointsWaterlevel (pointsFile )

#    # Plot maxeles maxele_nhctrk
#    f = plotter.plotMap    (params['lonlim'], params['latlim'], fig_w=10.)
#    plotter.addSurface (grid, maxele['value'],clim=params['clim'])
#    plt.plot(track_shiftLeft['lon'], track_shiftLeft['lat'],'o-k',markersize=1,zorder=10)    
#    title = 'HSOFS experimental ' + stormID + '.' + inputCycle + '.' + ens
#    plt.text (params['lonlim'][0]+0.03, \
#              params['latlim'][0]+0.03, \
#                  title )    
#    plotter.save(title, outputPath + '/' + stormID + '.' + inputCycle + '.'+ ens + '.maxele.png')
#    plt.close(f)
    
    # 5
    ens = 'higherSpeed'
    ncfile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.' + ens + '.fields.maxele.nc'
    trkFile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.'+ ens + '.surfaceforcing' 
    pointsFile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.'+ ens + '.points.waterlevel.nc'
    
    maxele = estofs.getFieldsWaterlevel (ncfile, 'zeta_max')    
    track_higherSpeed  = atcf.readTrack(trkFile)
    points_higherSpeed = estofs.getPointsWaterlevel (pointsFile )

#    # Plot maxeles maxele_nhctrk
#    f = plotter.plotMap    (params['lonlim'], params['latlim'], fig_w=10.)
#    plotter.addSurface (grid, maxele['value'],clim=params['clim'])
#    plt.plot(track_higherSpeed['lon'], track_higherSpeed['lat'],'o-k',markersize=1,zorder=10)    
#    title = 'HSOFS experimental ' + stormID + '.' + inputCycle + '.' + ens
#    plt.text (params['lonlim'][0]+0.03, \
#              params['latlim'][0]+0.03, \
#                  title )    
#    plotter.save(title, outputPath + '/' + stormID + '.' + inputCycle + '.'+ ens + '.maxele.png')
#    plt.close(f)

#    # MaxPS
#    maxPSfile = inputPath + '/hsofs.' + stormID + '.' + inputCycle + '.fields.maxPS.nc'
#    maxPS = estofs.getFieldsWaterlevel (maxPSfile, 'zeta_max')
#
#    # Plot maxeles maxPS
#    f = plotter.plotMap    (params['lonlim'], params['latlim'], fig_w=10.)
#    plotter.addSurface (grid, maxPS['value'],clim=params['clim'])
#    
#    plt.plot(trk['lon'], trk['lat'],'o-k',markersize=1,zorder=10)    
#    plt.plot(track_nhctrk['lon'],      track_nhctrk['lat'],     'o-k',markersize=1,zorder=10)    
#    plt.plot(track_higherSpeed['lon'], track_higherSpeed['lat'],'o-k',markersize=1,zorder=10)    
#    plt.plot(track_lowerSpeed['lon'],  track_lowerSpeed['lat'], 'o-k',markersize=1,zorder=10)    
#    plt.plot(track_shiftLeft['lon'],   track_shiftLeft['lat'],  'o-k',markersize=1,zorder=10)    
#    plt.plot(track_shiftRight['lon'],  track_shiftRight['lat'], 'o-k',markersize=1,zorder=10)    
#        
#    title = 'HSOFS experimental ' + stormID + '.' + inputCycle + '.maxPS'
#    plt.text (params['lonlim'][0]+0.03, \
#              params['latlim'][0]+0.03, \
#                  title )    
#    plotter.save(title, outputPath + '/' + stormID + '.' + inputCycle + '.maxPS.maxele.png')
#    plt.close(f)



    cwl = points_nhctrk # nhctrk stations
    mod_dates = cwl['time'] 
    
    # Plot time series
    utcnow = dt.utcnow()
    daterange = (utcnow-td(days=1), utcnow)
        
    figureCounter = 0
    print '[info]: Plotting time series...'
    
    for n in range(len(cwl['lon'])):
        
        print str(n).zfill(3) 
        
        if params['xlim'][0] <= cwl['lon'][n] and \
           cwl['lon'][n] <= params['xlim'][1] and \
           params['ylim'][0] <= cwl['lat'][n] and \
           cwl['lat'][n] <= params['ylim'][1]:
               
               mod_vals = np.squeeze( cwl['zeta'][:,n] )
    
               obs_vals = copy.deepcopy(mod_vals) #in case there is no obs
               obs_dates = copy.deepcopy(mod_dates) # in case if there is no obs
               
               coops_id = cwl['stations'][n].split()[0] #[2]
               print coops_id
               obs_coops = coops.getData ( coops_id, daterange)
               if len(obs_coops['values'])>0:
                   obs_dates = obs_coops['dates']
                   obs_vals  = obs_coops['values']
                                        
               # Plot observations
               f = plt.figure(figsize=(20,4.5))               
               plt.plot(obs_dates, obs_vals, '.',color='g',label='OBS')
               
               # Plot model
               plt.plot(points_nhctrk['time'],      points_nhctrk['zeta'][:,n],     '.', c='k',label='nhctrk')
               plt.plot(points_higherSpeed['time'], points_higherSpeed['zeta'][:,n],'.', c='r',label='higherSpeed')
               plt.plot(points_lowerSpeed['time'],  points_lowerSpeed['zeta'][:,n], '.', c='m',label='lowerSpeed')
               plt.plot(points_shiftRight['time'],  points_shiftRight['zeta'][:,n], '.', c='b',label='shiftRight')
               plt.plot(points_shiftLeft['time'],   points_shiftLeft['zeta'][:,n],  '.', c='b',label='shiftLeft')
                   
               plt.ylim([-1.5,params['clim'][1]])
               plt.title(cwl['stations'][n])
               plt.xlim([daterange[0],mod_dates[-1]])
               plt.grid()
               plt.xlabel('DATE UTC')
               plt.ylabel('WL, meters LMSL')
               plt.legend(bbox_to_anchor=(0.9, 0.35))

               figFile=outputPath + '/' + stormID + '.' + inputCycle + '.ts-' + str(figureCounter).zfill(3) + '.png'
               plt.savefig(figFile)
               plt.close(f)
               figureCounter += 1
               

#==============================================================================
if __name__ == "__main__":

    run_hsofs (sys.argv[1:])

    

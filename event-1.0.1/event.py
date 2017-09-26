# -*- coding: utf-8 -*-
"""
@author: Sergey.Vinogradov
"""
import sys, os
import argparse
import matplotlib
matplotlib.use('Agg',warn=False)
import matplotlib.pyplot as plt
import datetime
from datetime import datetime as dt
from datetime import timedelta as td
import numpy as np
import copy

#==============================================================================
def run_event (argv):

    outputPath  = ''
    toolkitPath = ''
    cfgFile     = ''
    cycle       = ''

    parser = argparse.ArgumentParser()
    parser.add_argument('-o','--outputPath',  required=False)
    parser.add_argument('-t','--toolkitPath', required=True)
    parser.add_argument('-c','--cfgFile',     required=True)
    parser.add_argument('-f','--fctCycle',    required=False)
    args = parser.parse_args()
    if args.outputPath:
        outputPath = args.outputPath
    if args.toolkitPath:
        toolkitPath = args.toolkitPath
    if args.cfgFile:
        cfgFile     = args.cfgFile
    if args.fctCycle:
        cycle    = args.fctCycle

    print 'event.py configured with :'
    print 'outputPath  =', outputPath
    print 'toolkitPath =', toolkitPath
    print 'cfgFile     =', cfgFile
    print 'fctCycle    =', cycle

    sys.path.insert(0, toolkitPath )
    import csdlpy

    if cycle == '':
        latest = csdlpy.estofs.latestForecast ()
    else:
        retroCycle = datetime.datetime(int(cycle[0:4]), int(cycle[4:6]), \
                                       int(cycle[6:8]), int(cycle[8:10]))
        latest = csdlpy.estofs.latestForecast (retroCycle)
    print '[info]: requesting ', latest

    params   = read_event_cfg (cfgFile)
    #stations = event_timeseries (params, outputPath, latest)
    #event_maxele      (params, outputPath, stations, latest)
    event_inundation  (params, outputPath, latest)

#==============================================================================
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
def event_timeseries(params, outputPath, latest):

    import csdlpy

    cwlFile = params['prodPath'] + params['domain'] + '.' + \
                       latest['yyyymmdd'] + '/' + params['prefix'] + '.' + \
                       latest['tHHz'] + '.points.cwl.nc'
    htpFile = params['prodPath'] + params['domain'] + '.' + \
                       latest['yyyymmdd'] + '/' + params['prefix'] + '.' + \
                       latest['tHHz'] + '.points.htp.nc'

    cwl = csdlpy.estofs.getPointsWaterlevel (cwlFile)
    htp = csdlpy.estofs.getPointsWaterlevel (htpFile)
    mod_dates = cwl['time']
    htp_dates = htp['time']

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
               obs_coops = csdlpy.obs.coops.getData ( coops_id, daterange)
               if len(obs_coops['values'])>0:
                   obs_dates = obs_coops['dates']
                   obs_vals  = obs_coops['values']

               f = plt.figure()
               figureCounter += 1
               csdlpy.plotter.plot_estofs_timeseries(obs_dates, obs_vals, \
                                              mod_dates, mod_vals, \
                                              figFile=outputPath + '/ts-' + \
                                              str(figureCounter).zfill(3) + '.png', \
                                              stationName=latest['yyyymmdd']+'.' +  \
                                                          latest['tHHz']+' for ' + \
                                                          cwl['stations'][n], \
                                              htp_dates=htp_dates, htp_vals=htp_vals, \
                                              daterange = daterange, \
                                              ylim = (-1.5,params['clim'][1]))
               plt.close(f)
    
    return {'lon' : cwl['lon'], 'lat' : cwl['lat']}

#==============================================================================
def event_maxele (params, outputPath, stations, latest):

    import csdlpy

    maxeleFile = params['prodPath'] + params['domain'] + '.' + \
                       latest['yyyymmdd'] + '/' + params['prefix'] + '.' + \
                       latest['tHHz'] + '.fields.cwl.maxele.nc'
    if os.path.exists (maxeleFile):
        maxele = csdlpy.estofs.getFieldsWaterlevel (maxeleFile, 'zeta_max')
    else:
        hourlyFields = params['prodPath'] + params['domain'] + '.' + \
                       latest['yyyymmdd'] + '/' + params['prefix'] + '.' + \
                       latest['tHHz'] + '.fields.cwl.nc'
        maxele = csdlpy.adcirc.computeMaxele(hourlyFields)

    gridFile = 'fort.14'
    trkFile = 'trk.dat'
    advFile = 'adv.dat'

    csdlpy.transfer.download (params['gridPath'], gridFile)
    csdlpy.transfer.download (params['bestTrackURL'], trkFile)
    csdlpy.transfer.download (params['advTrackURL'],  advFile)

    grid = csdlpy.adcirc.readGrid (gridFile)
    trk  = csdlpy.atcf.readTrack(trkFile)
    adv =  csdlpy.atcf.readTrack(advFile)

    csdlpy.plotter.plotMap    (params['lonlim'], params['latlim'], fig_w=10.)
    try:
        csdlpy.plotter.addSurface (grid, maxele['value'],clim=params['clim'])
    except:
        print '[warn]: cannot plot surface'
    plt.plot(trk['lon'], trk['lat'],'o-k',markersize=2,zorder=10)
    plt.plot(adv['lon'], adv['lat'],'o-r',markersize=2,zorder=11)


    for n in range(len(stations['lon'])):
        if params['xlim'][0] <= stations['lon'][n] and \
           stations['lon'][n] <= params['xlim'][1] and \
           params['ylim'][0] <= stations['lat'][n] and \
           stations['lat'][n] <= params['ylim'][1]:
               plt.plot(stations['lon'][n], stations['lat'][n],'wo', \
                        markeredgecolor='k', zorder=15)
    title = 'ESTOFS (GFS) ' + latest['yyyymmdd'] + '.' + latest['tHHz']
    plt.text (params['lonlim'][0]+0.02, \
              params['latlim'][0]+0.02, \
              title )
    csdlpy.plotter.save(title, outputPath + '/maxele.png')

#==============================================================================
def event_inundation (params, outputPath, latest):

    import csdlpy

    maxeleFile = params['prodPath'] + params['domain'] + '.' + \
                       latest['yyyymmdd'] + '/' + params['prefix'] + '.' + \
                       latest['tHHz'] + '.fields.cwl.maxele.nc'
    if os.path.exists (maxeleFile):
        maxele = csdlpy.estofs.getFieldsWaterlevel (maxeleFile, 'zeta_max')
    else:
        hourlyFields = params['prodPath'] + params['domain'] + '.' + \
                       latest['yyyymmdd'] + '/' + params['prefix'] + '.' + \
                       latest['tHHz'] + '.fields.cwl.nc'
        maxele = csdlpy.adcirc.computeMaxele(hourlyFields)

    gridFile = 'fort.14'

    csdlpy.transfer.download (params['gridPath'], gridFile)

    grid = csdlpy.adcirc.readGrid (gridFile)

    csdlpy.plotter.plotMap    (params['lonlim'], params['latlim'], fig_w=10.)
    field = maxele['value'] + grid['depth']
    
    field[np.where(np.isnan(maxele['value']))]=np.nan    
    field[np.where(grid['depth']>0.)]=np.nan
    #mask = np.any(np.where(np.isnan(field)))
    field = np.ma.masked_where(np.isnan(field), field, copy=False)
    
    try:
        csdlpy.plotter.addSurface (grid, 3.28*field, clim=[0.,6.0], zorder = 100)
    except:
        print '[warn]: cannot plot surface'

    title = 'ESTOFS (GFS) ' + latest['yyyymmdd'] + '.' + latest['tHHz']
    plt.text (params['lonlim'][0]+0.02, \
              params['latlim'][0]+0.02, \
              title )
    csdlpy.plotter.save(title, outputPath + '/inundation_agl_ft.png')


#==============================================================================
if __name__ == "__main__":

    run_event (sys.argv[1:])





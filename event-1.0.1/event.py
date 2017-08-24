# -*- coding: utf-8 -*-
"""
@author: Sergey.Vinogradov
"""
import sys
import argparse
import matplotlib
matplotlib.use('Agg',warn=False)
import matplotlib.pyplot as plt

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
    
    latest = estofs.latestForecast ()
    
    maxeleFile = params['prodPath'] + params['domain'] + '.' + \
                       latest['yyyymmdd'] + '/' + params['prefix'] + '.' + \
                       latest['tHHz'] + '.fields.cwl.maxele.nc'
    
    #maxeleFile = "C:\Users\sergey.vinogradov\Python\maxele.63.nc"
                          
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
    
    plotter.plotMap    (params['lonlim'], params['latlim'])   
    plotter.addSurface (grid, maxele['value'],clim=params['clim'])
    plt.plot(trk['lon'], trk['lat'],'o-k',markersize=2,zorder=10)    
    plt.plot(adv['lon'], adv['lat'],'o-r',markersize=2,zorder=11)    
    
    ncFile = params['prodPath'] + params['domain'] + '.' + \
                       latest['yyyymmdd'] + '/' + params['prefix'] + '.' + \
                       latest['tHHz'] + '.points.cwl.nc'
    #ncFile = "C:/Users/sergey.vinogradov/Downloads/estofs.atl.t12z.points.cwl.nc"
    print '[info]: reading ', ncFile
           
    points = estofs.getPointsWaterlevel (ncFile)
    
    for n in range(len(points['lon'])):
        if params['xlim'][0] <= points['lon'][n] and \
           points['lon'][n] <= params['xlim'][1] and \
           params['ylim'][0] <= points['lat'][n] and \
           points['lat'][n] <= params['ylim'][1]:
               plt.plot(points['lon'][n], points['lat'][n],'wo', \
                        markeredgecolor='k', zorder=15)
    title = 'ESTOFS-ATL HSOFS/GFS ' + latest['yyyymmdd'] + '.' + latest['tHHz']
    plt.text (params['lonlim'][0]+0.01, \
              params['latlim'][0]+0.01, \
              title )
    
    plotter.save(title, outputPath + '/maxele.png')        

#==============================================================================
if __name__ == "__main__":

    run_event (sys.argv[1:])

    


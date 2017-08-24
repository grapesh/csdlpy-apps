# -*- coding: utf-8 -*-
"""
@author: Sergey.Vinogradov
"""
import sys, os
from datetime import datetime as dt
import shutil
from datetime import timedelta
import argparse

##==============================================================================
def run_bias_table (argv):
   
    outputPath  = ''
    archivePath = ''
    toolkitPath = ''
    avgDays     = '5.0'

    parser = argparse.ArgumentParser()
    parser.add_argument('-o','--outputPath',  required=False)
    parser.add_argument('-a','--archivePath', required=False)
    parser.add_argument('-t','--toolkitPath', required=True)
    parser.add_argument('-d','--avgDays',     required=False)
    args = parser.parse_args()
    if args.outputPath:
        outputPath = args.outputPath
    if args.archivePath:
        archivePath = args.archivePath
    if args.toolkitPath:
        toolkitPath = args.toolkitPath
    if args.avgDays:
        avgDays     = int(round(float(args.avgDays)+0.499999))
        
    print 'bias_table.py configured with :'
    print 'outputPath  =', outputPath
    print 'archivePath =', archivePath
    print 'toolkitPath =', toolkitPath
    print 'avgDays     =', str(avgDays)
    
    csvFile = bias_table(outputPath, archivePath, toolkitPath, avgDays)
    bias_map (csvFile, outputPath, toolkitPath, avgDays)
    
#===============================================================================
def bias_table (outputPath, archivePath, toolkitPath, avgDays): 

    sys.path.insert(0, toolkitPath )
    from csdlpy.obs import coops

    now = dt.now()
    startDate = now-timedelta(days=avgDays-1)
    endDate   = now
    dates = (startDate, endDate)
        
    datestr  = dt.strftime(dt.now(),'%Y%m%d%H')
    csvFile  = os.path.join( outputPath, \
                               'biases.' + str(avgDays).zfill(3) + 'days.csv')

    print '[info]: Creating table: ' + csvFile
    coops.bias_table( csvFile, dates )

    if not archivePath == '':
        archiveFile = os.path.join( archivePath, 'biases.' + datestr + '.csv')
        print '[info]: Archiving to ' + archiveFile
        shutil.copy ( csvFile, archiveFile )

    return csvFile

#==============================================================================
def bias_map (csvFile, outputPath, toolkitPath, avgDays): 

    sys.path.insert(0, toolkitPath )
    from csdlpy.obs import coops
    from csdlpy import plotter
    
    x,y,z,span = coops.read_bias_table (csvFile)
    plotter.plotMap (x,y,fig_w=16.0,lonlim=(-180, 180),latlim=(-15,75))
    plotter.addTriangles((x,y,z))
    plotter.save    ( span, os.path.join(outputPath, \
                         'map-biases-' + str(avgDays).zfill(3) +'days.png') )
   
#==============================================================================
if __name__ == "__main__":

    run_bias_table (sys.argv[1:])

    


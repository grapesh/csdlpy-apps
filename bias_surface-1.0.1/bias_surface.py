# -*- coding: utf-8 -*-
"""
@author: Sergey.Vinogradov
"""
import sys, os
import numpy as np
import argparse
import matplotlib
matplotlib.use('Agg',warn=False)
import matplotlib.pyplot as plt

##==============================================================================
def run_bias_surface (argv):
   
    outputFile  = ''
    inputFile = ''
    toolkitPath = ''
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-o','--outputFile',  required=True)
    parser.add_argument('-a','--inputFile',   required=True)
    parser.add_argument('-t','--toolkitPath', required=True)
    parser.add_argument('-g','--gridFile',    required=True)
    parser.add_argument('-m','--mapFile',     required=True)
    
    args = parser.parse_args()
    if args.outputFile:
        outputFile = args.outputFile
    if args.inputFile:
        inputFile = args.inputFile
    if args.toolkitPath:
        toolkitPath = args.toolkitPath
    if args.gridFile:
        gridFile    = args.gridFile
    if args.mapFile:
        mapFile    = args.mapFile
        
    print 'bias_surface.py configured with :'
    print 'outputFile  =', outputFile
    print 'inputFile   =', inputFile
    print 'gridFile   =',  gridFile
    print 'mapFile    =',  mapFile
    print 'toolkitPath =', toolkitPath
            
    v = bias_surface (inputFile, outputFile, mapFile, toolkitPath, gridFile)
    
#===============================================================================
def bias_surface (inputFile, outputFile, mapFile, toolkitPath, gridFile): 

    sys.path.insert(0, toolkitPath )
    from csdlpy.obs import coops
    from csdlpy import transfer
    from csdlpy import adcirc
    from csdlpy import interp
    
    # Read bias table
    xo,yo,vo,datespan = coops.read_bias_table(inputFile)
    
    # Read / download grid
    if "ftp://" in gridFile or "http://" in gridFile or "https://" in gridFile:
        transfer.download (gridFile, 'fort.14')
        gridFile = 'fort.14'
    grid = adcirc.readGrid (gridFile)
            
    print '[test]: # Compute bias surface on the grid...'    
    vg = np.zeros(len(grid['depth']), dtype=float)
    
    # Interpolation parameters...............................................
    z_full  = 0.   # meters
    z_zero  = 200. # meters
    p       = 2.0  # scalar 
    R       = 2.0  # degrees
    #........................................................................
        
    print '[test]: ## Interpolate on the shelf...'
    ind_shelf     = np.where(grid['depth'] < z_zero)[0]
    vg[ind_shelf] = interp.shepard_idw (xo, yo, vo, grid['lon'][ind_shelf], 
                                                    grid['lat'][ind_shelf], p)

    print '[test]: ## Taper by depth...'
    ind_taper     = np.where (np.logical_and(z_full <= grid['depth'], 
                                             z_zero >= grid['depth']))[0]    
    vg[ind_taper] = interp.taper_linear (z_full, z_zero, 
                                          grid['depth'][ind_taper], 
                                          vg[ind_taper])
    
    print '[test]: ## Zero out the results that are too distant from data'
    dist = interp.distance_matrix(xo, yo, grid['lon'], grid['lat'])
    for n in range(len(grid['lon'])):
        if np.min(dist[:,n]) > R:
            vg[n] = 0.

    print '[test] :## Write out Offset file'
    adcirc.writeOffset63 (vg, outputFile)
    
    # Plot
    surface_map (grid, vg, toolkitPath, datespan, mapFile)
    
    return vg

#==============================================================================
def surface_map (grid, vg, toolkitPath, datespan, mapFile): 

    sys.path.insert(0, toolkitPath )
    from csdlpy import plotter
    
    print '[test]: ## Plot bias surface and data'
    plotter.plotMap    (grid['lon'], grid['lat'],   fig_w=10.0)
    plotter.addSurface (grid,    vg, clim=[-0.3, 0.3])
    plt.title(datespan,fontsize=8)
    plt.savefig(mapFile)
    plt.close()   

#==============================================================================
if __name__ == "__main__":

    run_bias_surface (sys.argv[1:])

    


######################################################################
# MOOSE script to assemble parts of Figure 6
# First, the data for two sequences is needed. Run the following scripts:
# python Fig6_generate_data.py --fnumber 0 --sequence 1234
# python Fig6_generate_data.py --fnumber 2 --sequence 40312
# Then, run this script
# python Fig6_assemble_data.py
######################################################################
import pylab
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import scipy.stats as ss

gridShape = (6, 2)

def plotBoilerplate( panelTitle, plotPos, xlabel = '', ylabel = '', xticks = [] ):
    ax = plt.subplot2grid( gridShape, plotPos, colspan =2 )
    if len( xticks ) > 0:
        ax.xaxis.set_ticks( xticks )
    #ax.locator_params( 
    ax.spines['top'].set_visible( False )
    ax.spines['right'].set_visible( False )
    ax.tick_params( direction = 'out' )
    #ax.set_xticklabels([])
    ax.set_xticklabels( xticks )
    ax.set_xlabel( xlabel )
    for tick in ax.xaxis.get_major_ticks():
        tick.tick2On = False
    for tick in ax.yaxis.get_major_ticks():
        tick.tick2On = False
    
    plt.ylabel( ylabel, fontsize = 14 )
    # alternate way of doing this separately.
    #plt.yaxis.label.size_size(16)
    #plt.title( 'B' )
    ax.text( -0.15, 1.05, panelTitle, fontsize = 18, weight = 'bold',
    transform=ax.transAxes )
    return ax

def panelB( fig ):
    dt = 0.0001
    npts = int(1 + 25/dt)
    tab = readXplot('Fig6_data.0.1234.xplot', 'Soma_Vm.0')
    #tab = tab[:50001]
    tab = tab[:npts]
    #ax = plotBoilerplate( 'B', (0,0), ylabel = 'Vm (mV)', xticks=range(6) )
    ax = plotBoilerplate( 'B', (0,0), ylabel = 'Vm (mV)', xticks=range(0,26,5) )
    t = np.array( range( len( tab ) ) ) * dt
    print "LEN = ", len( tab )
    plt.plot( t, tab )
    plt.xlabel( "Time (s)", fontsize = 14, horizontalalignment='left', position = (0.02,25) )

def panelC( fig ):
    dt = 0.05
    npts = int(1 + 26/dt)
    ax = plotBoilerplate( 'C', (1,0), ylabel = 'Ca ($\mu$M)', xticks=range(0, 26, 5) )
    tab = readXplot('Fig6_data.0.1234.xplot', 'psd3589_Ca.0')[:npts]
    t = np.array( range( len( tab ) ) ) * dt
    plt.plot( t, tab )
    tab = readXplot('Fig6_data.0.1234.xplot', 'spine3589_Ca.0')[:npts]
    plt.plot( t, tab )
    tab = readXplot('Fig6_data.0.1234.xplot', 'dend70_8_Ca.0')[:npts]
    plt.plot( t, tab )
    plt.xlabel( "Time (s)", fontsize = 14, horizontalalignment='left', position = (0.02,25) )

def dendZone( idx ):
    #idx = 9
    return 'dend70_' + str( idx ) + '_P_MAPK.0'

def panelD( fig ):
    dt = 0.05
    #npts = 421
    npts = int( 1 + 25/dt )
    ax = plotBoilerplate( 'D', (2,0), ylabel = 'MAPK-P ($\mu$M)', xticks = range( 0, 26, 5 ) )
    tab = readXplot('Fig6_data.0.1234.xplot', dendZone( 4 ) )[:npts]
    t = np.array( range( len( tab ) ) ) * dt
    plt.plot( t, tab, linewidth = 3 )
    tab = readXplot('Fig6_data.0.1234.xplot', dendZone( 5 ))[:npts]
    plt.plot( t, tab, linewidth = 3 )
    tab = readXplot('Fig6_data.0.1234.xplot', dendZone( 6 ))[:npts]
    plt.plot( t, tab, linewidth = 3 )
    tab = readXplot('Fig6_data.0.1234.xplot', dendZone( 8 ))[:npts]
    plt.plot( t, tab, linewidth = 3 )
    tab = readXplot('Fig6_data.0.1234.xplot', dendZone( 9 ))[:npts]
    plt.plot( t, tab, linewidth = 3 )
    ################# scrambled points here #################

    tab = readXplot('Fig6_data.2.40312.xplot', dendZone( 4 ))[:npts]
    plt.plot( t, tab, 'b--' )
    tab = readXplot('Fig6_data.2.40312.xplot', dendZone( 5 ))[:npts]
    plt.plot( t, tab, 'g--' )
    tab = readXplot('Fig6_data.2.40312.xplot', dendZone( 6 ))[:npts]
    plt.plot( t, tab, 'r--' )
    tab = readXplot('Fig6_data.2.40312.xplot', dendZone( 8 ))[:npts]
    plt.plot( t, tab, 'c--' )
    tab = readXplot('Fig6_data.2.40312.xplot', dendZone( 9 ))[:npts]
    plt.plot( t, tab, 'm--' )
    plt.xlabel( "Time (s)", fontsize = 14, horizontalalignment='left', position = (0.02,25) )


def panelEFGH( fig ):
    paoc = getSeqMatrix()
    #valmatrix = []
    for site, panelTitle in zip( range( 6 ), ['E','F','G','H','I','J'] ):
        valvector = paoc[site]
        '''
        valvector = []
        for i in range( 25 ):
            values = paoc[:,i]
            #valvector.append( ( values[0] - np.mean( values )) / max( values ) )
            valvector.append( values[0] )
        '''
        #valmatrix.append( valvector )
        print site, max(valvector)
        plotMatrix( fig, site, valvector, panelTitle )

def plotMatrix( fig, site, vec, panelTitle ):
    print 'COORDS = ', 3 + site%3 , site / 3
    ax = plt.subplot2grid( gridShape, (3 + site%3, site/3 ) )
    d = np.array( vec )
    d.shape = (5,5)
    temp = np.array( d[0] )
    d[0] = d[4]
    d[4] = temp
    temp = np.array( d[1] )
    d[1] = d[3]
    d[3] = temp
    #cax = plt.imshow( d,cmap=cm.cool, vmin = 0, vmax = 0.68, interpolation = 'none' )
    #cax = plt.imshow( d,cmap=cm.cool, interpolation = 'none' )
    cax = plt.imshow( d,cmap=cm.cool, vmin = 0, vmax = 0.9, interpolation = 'none', aspect = 'auto' )
    if site > 2:
        cbar = fig.colorbar( cax )
    if site %3 == 2: # bottom row
        plt.xlabel( "Interval (s)", fontsize = 14 )
    xticks = [str(x) for x in range( 1, 6 )]
    plt.xticks( range( len(xticks ) ), xticks )
    if ( site / 3) == 0:    # left row
        plt.ylabel( r"Spacing ($\mu$m)", fontsize = 14)
        #yticks = ['1', '2', '3', '4', '6']
    yticks = ['6', '4', '3', '2', '1']
    plt.yticks( range( len(yticks ) ), yticks )

    ax.text( -0.15, 1.1, panelTitle, fontsize = 18, weight = 'bold',
        transform=ax.transAxes )

def main():
    fig = plt.figure(figsize = (6,12), facecolor='white')
    panelB( fig )
    panelC( fig )
    panelD( fig )
    panelEFGH( fig )
    plt.tight_layout()
    plt.show()

def readXplot( fname, plotname ):
    tab = []
    flag = False
    try:
        with open( fname ) as fd:
            for line in fd:
                if len( line ) > 1:
                    if line[0:9] == '/plotname':
                        if line[10:-1] == plotname:
                            flag = True
                        else:
                            flag = False
                        #print line[11:20], flag
                    if  flag and not( line[0] == '%' or line[0] == '/'):
                        tab.append( float( line[:-1] ) * 1000 )
    except IOError:
        print "NOTAB"
        return 0
    return tab

def readVm( fname ):
    tab = []
    flag = False
    try:
        with open( fname ) as fd:
            for line in fd:
                if len( line ) > 1:
                    if line[0:9] == '/plotname':
                        if line[10:19] == 'Soma_Vm.0':
                            flag = True
                        else:
                            flag = False
                        #print line[11:20], flag
                    if  flag and not( line[0] == '%' or line[0] == '/'):
                        tab.append( float( line[:-1] ) * 1000 )
    except IOError:
        print "NOTAB"
        return 0
    return tab

########################################################################
# Grid setup below.
########################################################################
def getSeqMatrix():
    # Order will be: EFG (vertically): 3 runs in Zone 2, show variability
    # HIJ: Vertically, 2nd column: Zone 0, Zone 1 and Zone 4
    # Zone 2: From homework/SingleNeuronSeqLearn/2017Jan20_smallerSpine, ran
    # python 52matrix.py 2017Jan20_smallerSpine.[123]/ 2
    # to get the first 3 entries in the seq matrix.
    # For the next 3 entries it is:
    
    # Zone 0: ~/homework/SingleNeuronSeqLearn/2017Jan17_matrix_repeats
    # Ran python 52matrix.py 2017Jan17_MatrixRepeats.4/ 0
    # Zone 1: ~/homework/SingleNeuronSeqLearn/2017Jan19_smallSpine
    # Ran: python 52matrix.py 2017Jan19_smallSpine.1/ 1                 
    # Zone 3: ~/homework/SingleNeuronSeqLearn/2017Jan17_matrix_repeats
    # Ran python 52matrix.py 2017Jan17_MatrixRepeats.4/ 3

    sortedAoc = [
        [0.06, 0.01, 0.03, 0.00, -0.10, -0.04, -0.06, 0.11, -0.02, -0.09, 
            0.03, 0.30, 0.49, 0.46, 0.08, 0.03, 0.27, 0.62, 0.22, -0.05, 
            0.11, 0.44, 0.69, 0.13, 0.16],
        [0.07, -0.07, -0.14, -0.07, -0.05, -0.03, 0.14, 0.07, -0.08, 
            -0.03, 0.03, -0.07, 0.42, -0.08, 0.04, 0.07, 0.01, 0.57, 
            0.80, 0.10, 0.12, 0.59, 0.37, 0.35, -0.03],
        [-0.03, -0.08, -0.14, 0.34, 0.04, 0.00, 0.07, -0.13, 0.67, 0.05, 
            -0.01, -0.03, -0.02, 0.74, 0.08, -0.00, 0.23, 0.73, 0.71, 
            -0.03, 0.07, 0.59, 0.78, 0.17, 0.10],
        [ 0.00, 0.01, -0.09, -0.42, 0.01, -0.07, 0.02, 0.05, 0.08, -0.01,
            0.01, -0.00, 0.25, 0.33, 0.26, -0.04, -0.02, 0.14, 0.42, -0.05,
            0.02, 0.18, 0.15, 0.67, 0.27],
        [-0.06, 0.05, 0.06, 0.24, 0.01, 0.03, 0.08, 0.33, 0.55, -0.07, 0.09, 0.02, 0.11, 0.30, 0.17, 0.10, 0.29, -0.12, 0.82, 0.15, 0.24, 0.33, 0.75, -0.08, -0.08],
        [0.06, -0.01, -0.04, 0.03, 0.03, -0.02, 0.12, 0.03, -0.02, -0.02, 0.01, 0.05, 0.07, 0.03, -0.03, 0.05, -0.02, 0.04, 0.02, 0.02, 0.06, 0.10, -0.02, 0.08, -0.01]
        ]
    # Order is Zone2, Zone 2, Zone 2, Zone 0, Zone 1, Zone 4.
    return np.array( sortedAoc )

########################################################################
if __name__ == '__main__':
    main()

import moose
import pylab
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import cm
import rdesigneur as rd
import xml.etree.ElementTree as ET
import itertools
from scipy import stats

params = {
    'diffusionLength':0.5e-6,  # Diffusion characteristic length, used as voxel length too.
    'dendDiameter': 1.0e-6, # Diameter of section of dendrite in model
    'dendLength': 0.5e-6,    # Length of section of dendrite in model
    'diffConstCa':100e-12,  # Diffusion constant of Ca
    'stimAmplitude': 0.005, # Ca Stimulus amplitude, mM
    'baseCa':2.5e-4,        # Base Ca level, mM.
    'BAPCa':0.002,          # Dend Ca spike amplitude
    'BAPwidth':0.1,         # Dend Ca spike width.
    'blankVoxelsAtEnd':10,  # of voxels to leave blank at end of cylinder
    'preStimTime':10.0,     # Time to run before turning on stimulus.
    'postStimTime':40.0,    # Time to run after stimulus.
    'stimWidth': 2.9,       # Duration of Ca influx for each stimulus.
    'spineSpacing':1.1e-6,  # Spacing between spines.
    'diffConstMAPK': 5e-12, # Diffusion constant for MAPK
    'diffConstPP': 2e-12,   # Diff constant for MAPK-activated phosphatase
    'CaActivateRafKf': 6e6, # 1/sec/mM^2: rate for activation of Raf by Ca
    'cellModel':'PassiveSoma',  # Cell morphology script
    'chemModel':'NN_mapk14.g',  # Chem model definition
    'seqDt': 3.0,           # Time interval between successive inputs in seq
    'seqDx': 3.0e-6,        # Distance between successive inputs in seq.
    #'seed': 3,            # Seed for random number generator
    'seed': 12345,           # Seed for random number generator
    'sequence': '01234',    # Sequence of spines, spaced by seqDx microns, 
                            # activated every seqDt seconds
    'fnumber': 0,           # identifier for run
}
numSpine = 5

def makePassiveSoma( name, length, diameter ):
    elecid = moose.Neuron( '/library/' + name )
    dend = moose.Compartment( elecid.path + '/soma' )
    dend.diameter = diameter
    dend.length = length
    dend.x = length
    return elecid

def setDiffConst( element, paramName ):
    e = moose.element( '/library/chem/kinetics/DEND/' + element )
    e.diffConst = params[ paramName ]

def buildStimulusQ( baseCa, sequence ):
    blanks = params['blankVoxelsAtEnd']
    step = int( round( params['seqDx'] / params['spineSpacing'] ) )
    sequence = [ int(i) for i in sequence ]
    stimulusQ = {}
    onCa = params['stimAmplitude']
    stimStart = params['preStimTime']
    for i in sequence:
        stimulusQ[ stimStart ] = [ blanks + i * step, onCa ]
        stimEnd = stimStart + params['stimWidth']
        stimulusQ[ stimEnd ] = [ blanks + i * step, baseCa ]
        stimStart += params['seqDt']
    return stimulusQ

def panelBCsingleCompt( fig ):
    rdes = rd.rdesigneur(
        useGssa = False,
        turnOffElec = True,
        chemPlotDt = 1.0,
        diffusionLength = params['dendLength'],
        spineProto = [['makePassiveSpine()', 'spine']],
        spineDistrib = [['spine', '#', str(0.8*params['dendLength']),'1e-7','1.4','0']],
        cellProto = [['cell', 'soma']],
        chemProto = [[params['chemModel'], 'chem']],
        chemDistrib = [['chem', 'soma', 'install', '1' ]],
        plotList = [
            ['soma', '1', 'dend/DEND/P_MAPK', 'conc', '[dend P_MAPK]'],
            ['#head#', '1', 'psd/Ca', 'conc', '[PSD Ca]'],
            ['#head#', '1', 'spine/Ca', 'conc', '[spine Ca]'],
            ['soma', '1', 'dend/DEND/Ca', 'conc', '[dend Ca]'],
            ],
    )
    #moose.le( '/library/chem/kinetics' )
    moose.element( '/library/chem/kinetics/DEND/Ca_activate_Raf' ).Kf = params['CaActivateRafKf']
    rdes.buildModel()
    baseCa = params['baseCa'] / 50.0
    Ca_input = moose.vec( '/model/chem/psd/Ca_input' )
    Ca_input.concInit = baseCa
    moose.vec( '/model/chem/dend/DEND/Ca_input' ).concInit = baseCa
    moose.reinit()
    moose.start(50)
    Ca_input.concInit = params['stimAmplitude']
    moose.start(1)
    Ca_input.concInit = baseCa
    moose.start(49)
    # Here is another pulse stimulus
    Ca_input.concInit = params['stimAmplitude']
    moose.start(1)
    Ca_input.concInit = baseCa
    moose.start(99)
    # Here is the step stimulus
    Ca_input.concInit = params['stimAmplitude']
    moose.start(50)
    Ca_input.concInit = baseCa
    moose.start(50)
    mapkPvec = moose.element( '/model/graphs/plot0' ).vector * 1000
    cavec = moose.element( '/model/graphs/plot3' ).vector * 1000
    t = np.arange( 0, len( cavec ) * rdes.chemPlotDt, rdes.chemPlotDt )
    xticks = ['0', '', '100', '', '200', '', '300' ]
    ax = plotBoilerplate( 'B', (2,0), 'Time (s)', 'Ca ($\mu$M)', xticks )
    plt.plot( t, cavec ) 
    ax = plotBoilerplate( 'C', (2,1), 'Time (s)', 'MAPK-P ($\mu$M)', xticks )
    plt.plot( t, mapkPvec ) 
    print "Finished panelBC for single Compt dynamics"
    moose.delete( '/model' )

def runStimulus( sequence ):
    baseCa = params['baseCa'] / 50.0
    stimQ = buildStimulusQ( baseCa, sequence )
    Ca_input = moose.vec( '/model/chem/psd/Ca_input' )
    Ca_input.concInit = baseCa
    moose.vec( '/model/chem/dend/DEND/Ca_input' ).concInit = baseCa
    moose.reinit()

    clock = moose.element( '/clock' )
    for t in sorted( stimQ ):
        [index, conc] = stimQ[t]
        currt = clock.currentTime
        if ( t > currt ):
            moose.start( t - currt )
            print "At t = ", t, "; assigning CaInput[", index, "] = ", conc
            Ca_input[ index ].concInit = conc
    moose.start(params['postStimTime'] )
    print "Finished stimulus run at t = ", clock.currentTime

def runAndDisplaySumPlot( seq, label, pos, chemPlotDt ):
    runStimulus( seq )
    mapk = moose.vec( '/model/graphs/plot0' )
    mapkPvec = np.zeros( len( mapk[0].vector ) )
    for i in mapk:
        mapkPvec += i.vector
    mapkPvec *= 1000
    t = np.arange( 0, len( mapkPvec ) * chemPlotDt, chemPlotDt )
    xt = np.arange( 0, len( mapkPvec ) * chemPlotDt, 50 )
    xticks = [ str(i) for i in xt ]
    print "XT = ", xt, xticks
    ax = plotBoilerplate( label, pos, 'Time (s)', 'MAPK-P ($\mu$M)', xticks )
    plt.plot( t, mapkPvec ) 

def runAndDisplay( seq, label, pos, chemPlotDt, maxy ):
    runStimulus( seq )
    mapk = moose.vec( '/model/graphs/plot0' )
    parentVoxel = moose.element( '/model/chem/spine' ).neuronVoxel
    blanks = params['blankVoxelsAtEnd']
    step = int( round( params['seqDx'] / params['spineSpacing'] ) )
    plotIndices = [int(parentVoxel[blanks + i * step]) for i in range(len(seq)) ]
    print plotIndices
    t = np.arange( 0, len( mapk[0].vector ) * chemPlotDt, chemPlotDt )
    ax = plotBoilerplate( label, pos, 'Time (s)', 'MAPK-P ($\mu$M)')
    maxx = max( t )
    xt = np.arange( 0, maxx, 20 )
    xticks = [ str(int(i)) for i in xt ]
    print "XT = ", xt, xticks

    ax.xaxis.set_ticks( xt )
    ax.set_xticklabels( xticks )
    ax.set_ylim( [0.0, maxy] )
    for i in plotIndices:
        plt.plot( t, mapk[i].vector * 1000 )

def panelEFspatialSeq( fig ):
    print "Starting Panel EF"
    moose.seed( int(params['seed']) )
    rdes = rd.rdesigneur(
        useGssa = False,
        turnOffElec = True,
        chemPlotDt = 0.02,
        diffusionLength = params['diffusionLength'],
        spineProto = [['makePassiveSpine()', 'spine']],
        spineDistrib = [['spine', '#', str(params['spineSpacing']),'1e-7','1.4','0']],
        cellProto = [['cell', 'soma']],
        chemProto = [[params['chemModel'], 'chem']],
        chemDistrib = [['chem', 'soma', 'install', '1' ]],
        plotList = [
            ['soma', '1', 'dend/DEND/P_MAPK', 'conc', '[dend P_MAPK]'],
            ['#head#', '1', 'psd/Ca', 'conc', '[PSD Ca]'],
            ['#head#', '1', 'spine/Ca', 'conc', '[spine Ca]'],
            ['soma', '1', 'dend/DEND/Ca', 'conc', '[dend Ca]'],
            ],
    )
    # Assign parameters to the prototype model.
    setDiffConst( 'Ca', 'diffConstCa' )
    setDiffConst( '../../compartment_1/Ca', 'diffConstCa' )
    setDiffConst( '../../compartment_2/Ca', 'diffConstCa' )
    setDiffConst( 'P_MAPK', 'diffConstMAPK' )
    setDiffConst( 'MAPK', 'diffConstMAPK' )
    setDiffConst( 'reg_phosphatase', 'diffConstPP' )
    setDiffConst( 'inact_phosphatase', 'diffConstPP' )
    moose.element( '/library/chem/kinetics/DEND/Ca_activate_Raf' ).Kf = params['CaActivateRafKf']
    print "Set up rdesigneur"

    rdes.buildModel()
    print "MODEL BUILT"

    ################################################################
    # Run and display the stimulus
    runAndDisplay( '01234', 'E', (4,0), rdes.chemPlotDt, 0.6 )
    runAndDisplay( '40312', 'F', (4,1), rdes.chemPlotDt, 0.6 )

    moose.delete( '/model' )

def plotBoilerplate( panelTitle, plotPos, xlabel = '', ylabel = '', xticks = [] ):
    ax = plt.subplot2grid( (5,2), plotPos )
    #ax.xaxis.set_ticks( i[1] )
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
    ax.text( -0.3, 1.1, panelTitle, fontsize = 18, weight = 'bold',
    transform=ax.transAxes )
    return ax

def main():
    global params
    fig = plt.figure(figsize = (6,10), facecolor='white')

    library = moose.Neutral( '/library' )

    for ii in range( len( sys.argv ) ):
        if sys.argv[ii][:2] == '--':
            argName = sys.argv[ii][2:]
            if argName in params:
                params[argName] = float( sys.argv[ii+1] )
                if argName == 'sequence':
                    params[argName] = sys.argv[ii+1] # Leave it as a str.

    moose.seed( int(params['seed']) )
    '''
    '''
    makePassiveSoma( 'cell', params['dendLength'], params['dendDiameter'] )
    moose.le( '/library' )
    panelBCsingleCompt( fig );
    moose.le( '/library' )
    moose.delete( '/library/soma' )
    params['dendLength'] = 60e-6
    makePassiveSoma( 'cell', params['dendLength'], params['dendDiameter'] )
    panelEFspatialSeq( fig );
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()

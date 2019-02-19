################################################################
# This version generates a cylinder with the correct # of spines
# As of the current MOOSE, the random numbers aren't being handled
# independently, so I need to run a dummy single-compartment model
# to match the simulations in the rest of Figure 5.
################################################################
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
    'dendLength': 60e-6,    # Length of section of dendrite in model
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
    'seed': 12345,            # Seed for random number generator
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
    e = moose.element( '/library/chem/dend/DEND/' + element )
    e.diffConst = params[ paramName ]

def buildStimulusQ( sequence ):
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

def dummy():
    print "Starting Dummy"
    makePassiveSoma( 'cell', 0.5e-6, params['dendDiameter'] )
    moose.reinit()
    moose.seed( int(params['seed']) )
    rdes = rd.rdesigneur(
        useGssa = False,
        turnOffElec = True,
        chemPlotDt = 0.02,
        diffusionLength = params['diffusionLength'],
        spineProto = [['makePassiveSpine()', 'spine']],
        spineDistrib = [['spine', '#', '0.4e-6','1e-7','1.4','0']],
        cellProto = [['cell', 'soma']],
        chemProto = [[params['chemModel'], 'chem']],
        chemDistrib = [['chem', 'soma', 'install', '1' ]],
        #moogList = [ ['soma', '1', '.', 'Vm', 'Vm', -0.1, 0.05], ]
    )
    rdes.buildModel()
    moose.le( '/library' )
    moose.delete( '/library/soma' )
    moose.delete( '/library/chem' )
    moose.le( '/' )
    moose.delete( '/model' )
    print "Finsihed dummy "

def paneDlayout():
    print "Starting Panel D"
    moose.reinit()
    moose.seed( int(params['seed']) )

    print "MODEL PREBUILT"
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
        moogList = [ ['#', '1', '.', 'Vm', 'Vm', -0.1, 0.05], ]
    )
    rdes.buildModel()
    print "MODEL BUILT"
    sequence = [ int(i) for i in params['sequence']]
    blanks = params['blankVoxelsAtEnd']
    step = int( round( params['seqDx'] / params['spineSpacing'] ) )
    for i in moose.wildcardFind( '/model/elec/head#,/model/elec/shaft#' ):
        i.initVm = 0.0
    for i in sequence:
        e = moose.element( '/model/elec/head' + str(blanks + i * step) )
        e.initVm = 0.1
        e = moose.element( '/model/elec/shaft' + str(blanks + i * step) )
        e.initVm = 0.1

    ################################################################
    # Run and display the stimulus
    moose.reinit()
    rdes.displayMoogli( 0.0001, 0.001, 0.0 )

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

    dummy()
    makePassiveSoma( 'cell', params['dendLength'], params['dendDiameter'] )
    paneDlayout()

if __name__ == '__main__':
    main()

import moose
import pylab
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import rdesigneur as rd
import xml.etree.ElementTree as ET
#import itertools
from scipy import stats
import time

baseFname = 'Fig6_data'
displayMoogli = False
displayScatter = False
displayRuns = False
dumpSpatialPlotsFlag = False
scatterParamToUseInStats = 0
realStartTime = time.time()
rdes = 0
spikeStartIndices = [1131, 2000, 500, 1490]
# The dataRange maps to the spikeStartIndices, but is a dendrite index.
dataRange = [[1810,1850],[770,820],[3620,3680],[2510,2550]]
aoc = np.zeros( len( dataRange ) )

# Stim amplitude is unitless, defined in units of A mol conc.
# Stim Width is unitless, defined as multiple of diffusion length.
# Stim Vel is unitless, defined in terms of diffusion length by  time units of diff constt.
# diffConst is defined in terms of square of diffusion length by unit time.
# diffLength here is in SI units: m^2/s
#
params = {
    'diffusionLength':1.0e-6,  # Diffusion characteristic length, used as voxel length too.
    'dendDiameter': 1e-6,  # Diameter of section of dendrite in model
    'dendLength': 100e-6,   # Length of section of dendrite in model
    'spineSizeScale': 1.0,  # Length scaling for spines. Vol wil be x^3.
    'diffConstCa':100.0e-12,      # Diffusion constant of Ca
    'diffConstMAPK': 5e-12, # Diffusion constant for MAPK
    'diffConstPP': 2e-12, # Diff constant for MAPK-activated phosphatase
    'CaActivateRafKf': 6e6,	# 1/sec.mM: rate for activation of Raf by Ca
    'blankVoxelsAtEnd':10,  # of voxels to leave blank at end of cylinder
    'preStimTime':1.0,      # Time to run before turning on stimulus.
    'stimBurstTime':2.0,    # Time for a stimulus burst
    'postStimTime':10.0,    # Time to run after stimulus. ~3x decay time
    'runtime':20.0,         # Simulation run time
    'checkPoint':1.0,       # How often to do a checkpoint.
    'chemPlotDt':0.05,      # Plotting timestep for chemical signaling.
    'elecPlotDt':0.1e-3,    # Plotting timestep for electrical signaling.
    'spineSpacing':1.0e-6,  # Spacing between spines.
    'stimSpacing':4,        # Stimulus spacing, in terms of # of spines.
    'meanSpikeRate':0.1,    # Basal mean rate for all synapses.
    'activeSpikeRate':20.0,  # Active input rate on specified synapses.
    'baseGabaSpikeRate':1.0,    # 1 Hz.
    'thetaGabaSpikeRate':0.5, # This is the peak, not the average theta.
    'thetaFreq':8.0, 		# 
    'amparSynapseWeight': 30.0,    #
    'nmdarSynapseWeight': 30.0,   #
    'gabarSynapseWeight': 30.0,    #
    'LCaDensity': 0.0,     # Channel density for LCa channels.
    'adaptorScale':60.0e3, # Adaptor scale factor from conc to Na density in Seimens/m^2.
    'CaPsdScale': 0.08,    # Adaptor scale from psd elec Ca to chem conc.
    'Em': -60.0e-3,     #   Resting potential of neuron
    'refractoryPeriod':0.010,  # 10 ms refractory time.
    'cellModel': 'VHC-neuron.CNG.swc',  # Cell morphology file
    #'cellModel': 'h10.CNG.swc',  # Cell morphology file
    'chemModel': 'NN_mapk14.g',  # Chemical model file.
    'fnumber': 0,  # Output file index
    'seed': 1234,   # Seeds random numbers
    'sequence': 1234,  # stimulus sequence. 
    'seqDx': 4.0e-6,    # Sequence spatial interval
    'seqDt': 3.0,       # Sequence time interval
}

def attachStimulus( fname, rdes ):
    compts = rdes.elecid.compartments
    heads = []
    lookupSpineIndexFromCompt = []
    for i in compts:
        sl = rdes.elecid.spinesOnCompartment[i]
        lookupSpineIndexFromCompt.append( len( heads ) )
        heads.extend( sl[1::2] )

    ampar = moose.wildcardFind( '/model/elec/#/glu' )
    assert( len(ampar) == len( heads ) )
    moose.RandSpike( '/model/elec/spike', len( heads ) )
    spikes = moose.vec( '/model/elec/spike' )
    spikes.rate = params['meanSpikeRate']
    spikes.refractT = params['refractoryPeriod']
    amparSynWeight = params['amparSynapseWeight']
    nmdarSynWeight = params['nmdarSynapseWeight']

    for j in zip( heads, range( len( heads ) ) ):
        sh = moose.element( j[0].path + '/glu/sh' )
        sh.numSynapses = 1
        sh.synapse[0].weight = amparSynWeight
        moose.connect( spikes[j[1]], 'spikeOut', sh.synapse[0], 'addSpike' )
        sh = moose.element( j[0].path + '/NMDA/sh' )
        sh.numSynapses = 1
        sh.synapse[0].weight = nmdarSynWeight
        moose.connect( spikes[j[1]], 'spikeOut', sh.synapse[0], 'addSpike' )

    thetaGaba = moose.Function( '/model/elec/thetaGabaRhythm' )
    timeExpr = 'cos( 6.283 * t * ' + str(params['thetaFreq']/2.0) + ' )^2'
    rateExpr = str( 2.0 * params['thetaGabaSpikeRate'])
    thetaGaba.expr = str( params['baseGabaSpikeRate']) + ' + ' + rateExpr + ' * ' +  timeExpr
    print thetaGaba.expr,
    gabar = moose.wildcardFind( '/model/elec/#/GABA' )
    print "     NUM-GABA = ", len(gabar)
    rsGabar = moose.RandSpike( '/model/elec/gabaSpike', len( gabar ) )
    moose.connect( thetaGaba, 'valueOut', rsGabar, 'setRate', 'OneToAll' )
    gabaSpikes = moose.vec( '/model/elec/gabaSpike' )
    #spikes.rate = params['meanSpikeRate']
    gabaSpikes.refractT = params['refractoryPeriod']
    gabarSynWeight = params['gabarSynapseWeight']
    for i in range( len(gabar) ):
        sh = moose.element( gabar[i].path + '/sh' )
        sh.numSynapses = 1
        sh.synapse[0].weight = gabarSynWeight
        moose.connect( gabaSpikes[i], 'spikeOut', sh.synapse[0], 'addSpike' )

    return lookupSpineIndexFromCompt
    
def attachStimulus2( fname ):
    ampar = moose.wildcardFind( '/model/elec/#/glu' )
    nmdar = moose.wildcardFind( '/model/elec/#/NMDA' )
    numSyn = len( ampar )
    assert( numSyn == len( nmdar ) );
    moose.le( '/model' )

    moose.RandSpike( '/model/elec/spike', len( ampar ) )
    spikes = moose.vec( '/model/elec/spike' )
    spikes.rate = params['meanSpikeRate']
    spikes.refractT = params['refractoryPeriod']
    amparSynWeight = params['amparSynapseWeight']
    nmdarSynWeight = params['nmdarSynapseWeight']
    for i in range( numSyn ):
        sh = moose.element( ampar[i].path + '/sh' )
        sh.numSynapses = 1
        sh.synapse[0].weight = amparSynWeight
        moose.connect( spikes[i], 'spikeOut', sh.synapse[0], 'addSpike' )
        sh = moose.element( nmdar[i].path + '/sh' )
        sh.numSynapses = 1
        sh.synapse[0].weight = nmdarSynWeight
        moose.connect( spikes[i], 'spikeOut', sh.synapse[0], 'addSpike' )

def dumpPlots( fname ):
    currt = moose.element( '/clock' ).currentTime
    #print "os.rename( " + fname + ".xplot", fname + ".old )"
    os.rename( fname + ".xplot", fname + ".old" )
    print 'Current sim time = ', currt, ", user time = ", time.time() - realStartTime
    sys.stdout.flush()
    fd = open( fname + ".xplot", "w" )
    for ii in params:
        fd.write('% ' + ii + '  ' + str( params[ii] ) + '\n' )
    fd.flush()
    fd.close()

    for i in rdes.plotNames:
        # Args are tabname, plot title, index, scale, units 
        for j in moose.vec( i[0] ):
            plotLabel = i[1].replace( " ", "_" ) + "." + str(j.index)
            j.xplot( fname + ".xplot", plotLabel )

def dumpSpatialPlots( fname ):
    global aoc
    fd2 = open( fname + "_spatial.xplot", "a" )
    currt = moose.element( '/clock' ).currentTime
    elist = {}
    elist['spineCa'] = moose.vec( '/model/chem/spine/Ca' )
    elist['Ca4_CaM'] = moose.vec( '/model/chem/dend/DEND/CaM/CaM_dash_Ca4' )
    elist['dendCa'] = moose.vec( '/model/chem/dend/DEND/Ca' )
    #elist['dendKA_p'] = moose.vec( '/model/chem/dend/DEND/K_A_p' )
    elist['dendMAPK_p'] = moose.vec( '/model/chem/dend/DEND/P_MAPK' )
    #elist['dendPP'] = moose.vec( '/model/chem/dend/DEND/phosphatase' )
    #elist['dend_regPP'] = moose.vec( '/model/chem/dend/DEND/reg_phosphatase' )

    if dumpSpatialPlotsFlag:
        for key, elm in elist.items():
            fd2.write( "\n/newplot" )
            fd2.write( "\n/plotname " + key + "_" + str(currt) + "\n" )
            for v in elm.conc:
                fd2.write( str( v ) + "\n" )
            fd2.flush()
        fd2.close()

    ### Here we dump the running total for Area of Curve.
    conc = elist['dendMAPK_p'].conc
    fd3 = open( fname + ".aoc", "a" )
    fd3.write( str(currt) )
    for i in range( len( dataRange ) ):
        aoc[i] += sum( conc[ dataRange[i][0]: dataRange[i][1] ] )
        fd3.write( '\t' + str(aoc[i] ) )
    fd3.write( '\n' )
    fd3.flush()
    fd3.close()


def convertSeq( arg ):
    x = int( arg )
    ret = [0] * 5
    for i in range( 5 ):
        ret[-i-1] = x % 10
        x /= 10
    return ret

def printStuff( sequence, firstSpineOnCompt, rdes ):
    for i in sequence:
        for j in spikeStartIndices:
            spikeIdx = firstSpineOnCompt[j] + i * params['stimSpacing']
            headPath = '/model/elec/head' + str( spikeIdx )
            head = moose.element( headPath )
            #headCa = moose.element( headPath + '/Ca_conc' )
            #headGlu = moose.element( headPath + '/glu' )
            #headNmda = moose.element( headPath + '/NMDA' )
            print i, j, spikeIdx, head.name, rdes.elecid.parentCompartmentOfSpine[head].name

def main():
    global rdes
    global params
    sequence = [0,1,2,3,4]

    #print params['sequence'], sequence
    for ii in range( len( sys.argv ) ):
        if sys.argv[ii][:2] == '--':
            argName = sys.argv[ii][2:]
            if argName in params:
                if type( params[argName] ) is str:
                    params[argName] = sys.argv[ii+1]
                elif type( params[argName] ) is int:
                    params[argName] = int( sys.argv[ii+1] )
                else:
                    params[argName] = float(sys.argv[ii+1])

    params[ 'stimSpacing' ] = int( round(params['seqDx']/params['spineSpacing']) )
    if params['stimBurstTime'] >= params['seqDt']:
        params['stimBurstTime'] = params['seqDt']
        
    # Type correction hack
    #params['stimSpacing'] = int( params['stimSpacing'] ) 
    #params['fnumber'] = int( params['fnumber'] ) 
    #params['sequence'] = int( params['sequence'] ) 

    sequence = convertSeq( params['sequence'] )
    fname = baseFname + '.' +  str( int( params['fnumber'] ) ) + '.' + str( int( params['sequence'] ) )
    print params['sequence'], sequence, fname

    diffusionLength = params['diffusionLength']
    dendLength = params['dendLength']
    diffusionLength = params['diffusionLength']
    library = moose.Neutral( '/library' )
    if displayMoogli:
        ml = [
                ['#', '1', '.', 'Vm', 'Memb potential'],
                ['#', '1', 'Ca_conc', 'Ca', '[Ca]', 0, 10, 0.001],
            ]
    else:
        ml = []

    chanpath = os.path.dirname( os.path.realpath(__file__)) + '/proto21.'
    moose.seed( params['seed'] )
    rdes = rd.rdesigneur(
        useGssa = False,
        turnOffElec = False,
        chemPlotDt = params['chemPlotDt'],
        diffusionLength = diffusionLength,
        spineProto = [['makeExcSpine()', 'spine']],
        chanProto = [
            [ chanpath + 'make_K_AHP()', 'K_AHP' ],
            [ chanpath + 'make_K_A()', 'K_A' ],
            [ chanpath + 'make_K_C()', 'K_C' ],
            [ chanpath + 'make_K_DR()', 'K_DR' ],
            [ chanpath + 'make_Na()', 'Na' ],
            [ chanpath + 'make_Ca_conc()', 'Ca_conc' ],
            [ chanpath + 'make_Ca()', 'Ca' ],
            [ chanpath + 'make_NMDA()', 'NMDA' ],
            [ chanpath + 'make_glu()', 'glu' ],
            [ chanpath + 'make_GABA()', 'GABA' ],
        ],
        chemProto = [[params['chemModel'], 'chem']],
        cellProto = [[params['cellModel'], 'elec']],
        chanDistrib = [
            ["Ca_conc", "#", "tau", "0.0133" ],
            ["Ca", "#dend#,#basal#,#apical#", "Gbar", str( params["LCaDensity"] ) ],
            ["Ca", "#soma#", "Gbar", "40" ],
            ["Na", "#dend#,#basal#", "Gbar", "60" ],
            ["Na", "#soma#", "Gbar", "600" ],
            ["Na", "#apical#", "Gbar", "40+40*exp(-p/200e-6)" ],
            ["K_DR", "#dend#,#basal#", "Gbar", "(p < 400e-6)*200" ],
            ["K_DR", "#soma#", "Gbar", "250" ],
            ["K_DR", "#apical#", "Gbar", "60+40*(p < 125e-6)" ],
            ["K_AHP", "#", "Gbar", "8" ],
            ["K_C", "#basal#,#dend#,#apical#", "Gbar", "50+150*exp(-p/200e-6)" ],
            ["K_C", "#soma#", "Gbar", "100" ],
            ["K_A", "#soma#", "Gbar", "50" ],
            ["K_A", "#dend#,#apical#", "Gbar", "50*(1 + 2.0e-6/(dia + 0.1e-6))" ],
            ["GABA", "#apical#,#dend#,#basal#", "Gbar", "10 + 30*(p < 125e-6)" ],
        ],
        spineDistrib = [['spine','#dend#,#apical#', str(params['spineSpacing']),'1e-7', str( params['spineSizeScale'] ), '0.0' ]],
        chemDistrib = [['chem', '#', 'install', '1' ]],
        adaptorList = [
            [ 'Ca_conc', 'Ca', 'psd/Ca_input', 'concInit', 2e-6, params['CaPsdScale'] ],
            [ 'Ca_conc', 'Ca','dend/DEND/Ca_input','concInit',2.0e-6,0.0001],
            [ 'dend/DEND/channel_p', 'conc', 'Na', 'modulation', 1.0, params['adaptorScale']],
        ],
        plotList = [
            ['soma', '1', '.', 'Vm', 'Soma Vm'],
            ['head1732', '1', 'Ca_conc', 'Ca', 'head1732 eCa'],
            ['apical_36_2', '1', 'dend/DEND/P_MAPK', 'conc', 'dend36_2_P_MAPK'],
            ['apical_36_2', '1', 'dend/DEND/Ca', 'conc', 'dend36_2 Ca'],
            ['head1732', '1', 'spine/Ca', 'conc', 'spine1732_Ca'],
            ['head1732', '1', 'psd/Ca', 'conc', 'psd1732_Ca'],

            ['apical_36_3', '1', 'dend/DEND/P_MAPK', 'conc', 'dend36_3_P_MAPK'],
            ['apical_36_3', '1', 'dend/DEND/Ca', 'conc', 'dend36_3_Ca'],
            ['head1735', '1', 'spine/Ca', 'conc', 'spine1735_Ca'],
            ['head1735', '1', 'psd/Ca', 'conc', 'psd1735_Ca'],

            ['apical_36_5', '1', 'dend/DEND/P_MAPK', 'conc', 'dend36_5_P_MAPK'],
            ['apical_36_5', '1', 'dend/DEND/Ca', 'conc', 'dend36_5_Ca'],
            ['head1738', '1', 'spine/Ca', 'conc', 'spine1738_Ca'],
            ['head1738', '1', 'psd/Ca', 'conc', 'psd1738_Ca'],

            ['apical_36_6', '1', 'dend/DEND/P_MAPK', 'conc', 'dend36_6_P_MAPK'],
            ['apical_36_6', '1', 'dend/DEND/Ca', 'conc', 'dend36_6_Ca'],
            ['head1741', '1', 'spine/Ca', 'conc', 'spine1741_Ca'],
            ['head1741', '1', 'psd/Ca', 'conc', 'psd1741_Ca'],

            ['apical_36_7', '1', 'dend/DEND/P_MAPK', 'conc', 'dend36_7_P_MAPK'],
            ['apical_36_7', '1', 'dend/DEND/Ca', 'conc', 'dend36_7_Ca'],
            ['head1744', '1', 'spine/Ca', 'conc', 'spine1744_Ca'],
            ['head1744', '1', 'psd/Ca', 'conc', 'psd1744_Ca'],


            ['apical_70_4', '1', 'dend/DEND/P_MAPK', 'conc', 'dend70_4_P_MAPK'],
            ['apical_70_4', '1', 'dend/DEND/Ca', 'conc', 'dend70_4_Ca'],
            ['head3573', '1', 'spine/Ca', 'conc', 'spine3573_Ca'],
            ['head3573', '1', 'psd/Ca', 'conc', 'psd3573_Ca'],

            ['apical_70_5', '1', 'dend/DEND/P_MAPK', 'conc', 'dend70_5_P_MAPK'],
            ['apical_70_5', '1', 'dend/DEND/Ca', 'conc', 'dend70_5_Ca'],
            ['head3577', '1', 'spine/Ca', 'conc', 'spine3577_Ca'],
            ['head3577', '1', 'psd/Ca', 'conc', 'psd3577_Ca'],

            ['apical_70_6', '1', 'dend/DEND/P_MAPK', 'conc', 'dend70_6_P_MAPK'],
            ['apical_70_6', '1', 'dend/DEND/Ca', 'conc', 'dend70_6_Ca'],
            ['head3581', '1', 'spine/Ca', 'conc', 'spine3581_Ca'],
            ['head3581', '1', 'psd/Ca', 'conc', 'psd3581_Ca'],

            ['apical_70_8', '1', 'dend/DEND/P_MAPK', 'conc', 'dend70_8_P_MAPK'],
            ['apical_70_8', '1', 'dend/DEND/Ca', 'conc', 'dend70_8_Ca'],
            ['head3585', '1', 'spine/Ca', 'conc', 'spine3585_Ca'],
            ['head3585', '1', 'psd/Ca', 'conc', 'psd3585_Ca'],

            ['apical_70_9', '1', 'dend/DEND/P_MAPK', 'conc', 'dend70_9_P_MAPK'],
            ['apical_70_9', '1', 'dend/DEND/Ca', 'conc', 'dend70_9_Ca'],
            ['head3589', '1', 'spine/Ca', 'conc', 'spine3589_Ca'],
            ['head3589', '1', 'psd/Ca', 'conc', 'psd3589_Ca'],


            ['apical_16_15', '1', 'dend/DEND/P_MAPK', 'conc', 'dend16_15_P_MAPK'],
            ['apical_16_15', '1', 'dend/DEND/Ca', 'conc', 'dend16_15_Ca'],
            ['head733', '1', 'spine/Ca', 'conc', 'spine733_Ca'],
            ['head733', '1', 'psd/Ca', 'conc', 'psd733_Ca'],

            ['apical_52_12', '1', 'dend/DEND/P_MAPK', 'conc', 'dend52_12_P_MAPK'],
            ['apical_52_12', '1', 'dend/DEND/Ca', 'conc', 'dend52_12_Ca'],
            ['head2425', '1', 'spine/Ca', 'conc', 'spine2425_Ca'],
            ['head2425', '1', 'psd/Ca', 'conc', 'psd2425_Ca'],

        ],
        moogList = ml
    )

    ############## Set spine dimensions ##########################
    # This is now done in the rdesigneur spineDistrib
    #moose.element( '/library/spine/shaft' ).diameter *= params['shaftDiaScale']
    #moose.element( '/library/spine/head' ).diameter *= params['headDiaScale']
    ############## Set Ca diffusion const ##########################
    for ca in moose.wildcardFind( '/library/##/Ca[ISA=PoolBase]' ):
        ca.diffConst = params['diffConstCa']

    ############## Set MAPK diffusion const ##########################
    temp = params['diffConstMAPK']
    moose.element( '/library/chem/dend/DEND/P_MAPK' ).diffConst = temp
    moose.element( '/library/chem/dend/DEND/MAPK' ).diffConst = temp
    ############## Set PP diffusion const ##########################
    temp = params['diffConstPP']
    moose.element( '/library/chem/dend/DEND/reg_phosphatase' ).diffConst = temp
    moose.element( '/library/chem/dend/DEND/inact_phosphatase' ).diffConst = temp
    ############## Set resting potential ##########################
    for i in moose.wildcardFind( "/library/##[][ISA=CompartmentBase]" ):
        i.Em = params[ 'Em' ]
        i.initVm = params[ 'Em' ]
    ############## Set sensitivity to Ca ##########################
    moose.element( '/library/chem/dend/DEND/Ca_activate_Raf' ).Kf = params['CaActivateRafKf']

    #################### Build the model ##########################
    rdes.buildModel()
    #moose.showfields( '/model/chem/dend/DEND/PKA/KA' )
    #moose.le( '/model/elec' )
    #print "NUM PLOTS = ", len( rdes.plotNames)
    firstSpineOnCompt = attachStimulus( fname, rdes )

    fd = open( fname + ".xplot", "w" )
    fd.flush()
    fd2 = open( fname + "_spatial.xplot", "w" )
    fd.close()

    checkPoint = params['checkPoint']
    moose.setClock( 19, checkPoint )
    periodicOutput = moose.PyRun( '/model/graphs/periodicOutput' )
    periodicOutput.tick = 19
    periodicOutput.runString = 'dumpPlots( "' + fname + '" )'
    periodicOutput.initString = 'print "Starting PeriodicOutput"'

    spatOutput = moose.PyRun( '/model/graphs/spatOutput' )
    spatOutput.tick = 18
    spatOutput.runString = 'dumpSpatialPlots( "' + fname + '" )'
    spatOutput.initString = 'print "Starting Spatial Output"'
    moose.reinit()
    moose.seed( 1 )

    printStuff( sequence, firstSpineOnCompt, rdes )

    spikes = moose.vec( '/model/elec/spike' )

    runtime = params['runtime']
    if displayMoogli:
        rdes.displayMoogli( 0.05, runtime, 0.0 )
    moose.start( params['preStimTime'] )
    '''
    ni = len( spikes ) / 8.0
    spikeStartIndices = [ int( (i + 0.5) * ni ) for i in range( 8 ) ]
    print "########## len(spikes) = ", len(spikes)
    print "########## spikeStartIndices = ", spikeStartIndices
    '''

    stimIntervalTime = params['seqDt'] - params['stimBurstTime']
    assert( stimIntervalTime >= 0 )
    for i in sequence:
        for j in spikeStartIndices:
            spikeIdx = firstSpineOnCompt[j] + i * params['stimSpacing']
            spikes[spikeIdx].rate = params['activeSpikeRate']
        moose.start( params['stimBurstTime'] )
        for j in spikeStartIndices:
            spikeIdx = firstSpineOnCompt[j] + i * params['stimSpacing']
            spikes[spikeIdx].rate = params['meanSpikeRate']
        if ( stimIntervalTime > 0.0 ):
            moose.start( stimIntervalTime )
    moose.start( params['postStimTime'] )

    quit()

if __name__ == '__main__':
    main()





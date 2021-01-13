import FWCore.ParameterSet.Config as cms
#
# attention: default is changed to work on unsuppressed digis!! ##############
#
simEcalTriggerPrimitiveDigis = cms.EDProducer("EcalTrigPrimProducer",
    BarrelOnly = cms.bool(False),
    InstanceEB = cms.string(''),
    InstanceEE = cms.string(''),
    binOfMaximum = cms.int32(6), ## optional from release 200 on, from 1-10
    Famos = cms.bool(False),
    TcpOutput = cms.bool(False),
    Debug = cms.bool(False),
    Label = cms.string('simEcalUnsuppressedDigis'),
    oddWeightsTxtFile = cms.string(''),
    TPinfoPrintout = cms.bool(False),
    TPmode = cms.string('Run2') ##-- Mode to run in. Default to Run2 
)

# print "Ecal debug: leaving ecalTriggerPrimitiveDigis_readDBOffline_cfi.py"
# print "simEcalTriggerPrimitiveDigis",simEcalTriggerPrimitiveDigis
# Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
from AthenaCommon.SystemOfUnits import GeV
import sys
from argparse import ArgumentParser


ComponentAccumulator.debugMode="trackCA trackEventAlgo ..."

def TrackingHitInputToolCfg(flags):

    acc = ComponentAccumulator()
    TrackingHitInputTool = CompFactory.TrackingHitInputTool()
    acc.setPrivateTools(TrackingHitInputTool)

    return acc

def TrackingHitMapToolCfg(flags):

    acc = ComponentAccumulator()
    TrackingHitMapTool = CompFactory.TrackingHitMapTool()
    acc.setPrivateTools(TrackingHitMapTool)

    return acc


def TrackConversionToolCfg(flags):

    acc = ComponentAccumulator()

    from ActsConfig.ActsEventCnvConfig import ActsToTrkConverterToolCfg

    TrackConversionTool = CompFactory.TrackConversionTool('TrackConversionTool',
                                                    HitMapTool=acc.popToolsAndMerge(TrackingHitMapToolCfg(flags)),
                                                    converterTool=acc.popToolsAndMerge(ActsToTrkConverterToolCfg(flags))
                                                    )

    acc.setPrivateTools(TrackConversionTool)

    return acc

def TrackingRecoToolCfg(flags):


    acc = ComponentAccumulator()
    #TrackConversionTool = acc.popToolsAndMerge(TrackingHitInputToolCfg(flags))
    TrackingRecoTool = CompFactory.TrackingRecoTool('TrackingRecoTool')

    acc.setPrivateTools(TrackingRecoTool)

    return acc


def TrackingAlgCfg(flags):

    acc = ComponentAccumulator()

    from PixelGeoModelXml.ITkPixelGeoModelConfig import ITkPixelReadoutGeometryCfg
    acc.merge(ITkPixelReadoutGeometryCfg(flags))

    alg = CompFactory.TrackingAlg("TrackingAlg",
                                  HitInputTool=acc.popToolsAndMerge(TrackingHitInputToolCfg(flags)),
                                  HitMapTool=acc.popToolsAndMerge(TrackingHitMapToolCfg(flags)),
                                  RecoTool=acc.popToolsAndMerge(TrackingRecoToolCfg(flags)),
                                  TrackConversionTool=acc.popToolsAndMerge(TrackConversionToolCfg(flags))
                                 )

    acc.addEventAlgo(alg)

    ################################################################################
    # Convert ActsTrk::TrackContainer to xAOD::TrackParticleContainer
    #prefix = flags.Tracking.ActiveConfig.extension
    acts_tracks="ActsTracccTracks"
    from ActsConfig.ActsTrackFindingConfig import ActsTrackToTrackParticleCnvAlgCfg
    acc.merge(ActsTrackToTrackParticleCnvAlgCfg(flags, "TracccTrackToAltTrackParticleCnvAlg",
                                                ACTSTracksLocation=[acts_tracks],
                                                TrackParticlesOutKey="TracccTrackParticles"))
    """
    from ActsConfig.ActsTruthConfig import ActsTrackParticleTruthDecorationAlgCfg
    acc.merge(ActsTrackParticleTruthDecorationAlgCfg(flags,
                                                    "TracccActsTrackParticleTruthDecorationAlg",
                                                    TrackToTruthAssociationMaps=[acts_tracks+"TracccToTruthParticleAssociation"],
                                                    TrackParticleContainerName="TracccTrackParticles",
                                                    TruthParticleHitCounts="TracccTruthParticleHitCounts",
                                                    ComputeTrackRecoEfficiency=True))
    """
    #Adding the output to the AOD file
    inputList = []
    inputList.append("xAOD::TrackStatesContainer#*")
    inputList.append("xAOD::TrackStatesAuxContainer#*")
    inputList.append("xAOD::TrackJacobianContainer#*")
    inputList.append("xAOD::TrackJacobianAuxContainer#*")
    inputList.append("xAOD::TrackMeasurementContainer#*")
    inputList.append("xAOD::TrackMeasurementAuxContainer#*")
    inputList.append("xAOD::TrackSurfaceContainer#*")
    inputList.append("xAOD::TrackSurfaceAuxContainer#*")
    inputList.append("xAOD::TrackParticleContainer#*")
    inputList.append("xAOD::TrackParticleAuxContainer#*")
    acc.merge(OutputStreamCfg(flags, 'AOD', ItemList=inputList))


    return acc

if __name__=="__main__":
    from AthenaCommon.Logging import log
    from AthenaCommon.Constants import INFO
    log.setLevel(INFO)

    from AthenaConfiguration.AllConfigFlags import initConfigFlags
    flags = initConfigFlags()
    flags.Input.Files = ['/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/PhaseIIUpgrade/RDO/ATLAS-P2-RUN4-01-01-00/mc21_14TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.recon.RDO.e8481_s4038_r14365/RDO.32520484._000001.pool.root.1']
    flags.Input.isMC = True
    #flags.Input.ProjectName = "mc16_13TeV"
    #flags.Input.RunNumbers = [350200] # MC23 PhaseII mu=200 run number
    #flags.Input.TimeStamps = [1625130000] # MC23 PhaseII mu=200 time stamp
    #flags.IOVDb.GlobalTag = "PixelXDD-02-00"
    #from AthenaConfiguration.TestDefaults import defaultGeometryTags
    #flags.GeoModel.AtlasVersion = "ATLAS-P2-RUN4-03-00-00"

    flags.Detector.GeometryITkPixel = True
    flags.Detector.GeometryITkStrip = True
    flags.lock()

    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    cfg = MainServicesCfg(flags)

    from AtlasGeoModel.GeoModelConfig import GeoModelCfg
    cfg.merge(GeoModelCfg(flags))

    cfg.merge(TrackingAlgCfg(flags))

    cfg.run(maxEvents=1)

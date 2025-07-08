# Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg

inputDirectory = '/eos/project/a/atlas-eftracking/GPU/ITk_data/ATLAS-P2-RUN4-03-00-01/'
doTruth = True
doTruthSeeding = False
checkSeeds = False
doRecoFromClusters = False
makeStandaloneFiles = False
outputDirectory = ''
trackingRecoToolImpl = "CUDA"
writeOutput = True

ComponentAccumulator.debugMode="trackCA trackEventAlgo ..."

def TrackingHitInputToolCfg(flags):

    acc = ComponentAccumulator()
    TrackingHitInputTool = CompFactory.TrackingHitInputTool('TrackingHitInputTool',
                                                            filesDir=inputDirectory)
    acc.setPrivateTools(TrackingHitInputTool)

    return acc


def TrackConversionToolCfg(flags):

    acc = ComponentAccumulator()

    TrackConversionTool = CompFactory.TrackConversionTool('TrackConversionTool',
                                                        doTruth=doTruth,
                                                        recoFromHits=(not doRecoFromClusters),
                                                        filesDir=inputDirectory)

    acc.setPrivateTools(TrackConversionTool)

    return acc

def TrackingRecoToolCfg(flags):
    acc = ComponentAccumulator()

    if trackingRecoToolImpl.lower() == 'cuda':
        TrackingRecoTool = CompFactory.TrackingRecoTool_TrackingRecoToolTraitCUDA_(
            'TrackingRecoTool',
            filesDir=inputDirectory,
            TrackConversionTool=acc.popToolsAndMerge(TrackConversionToolCfg(flags)),
            doTruthSeeding=doTruthSeeding,
            HitInputTool=acc.popToolsAndMerge(TrackingHitInputToolCfg(flags)),
            checkSeeds=checkSeeds)
    elif trackingRecoToolImpl.lower() == 'sycl':
        TrackingRecoTool = CompFactory.TrackingRecoTool_TrackingRecoToolTraitSYCL_(
            'TrackingRecoTool',
            filesDir=inputDirectory,
            TrackConversionTool=acc.popToolsAndMerge(TrackConversionToolCfg(flags)),
            doTruthSeeding=doTruthSeeding,
            HitInputTool=acc.popToolsAndMerge(TrackingHitInputToolCfg(flags)),
            checkSeeds=checkSeeds)

    acc.setPrivateTools(TrackingRecoTool)
    return acc

def TrackingAlgCfg(flags):

    acc = ComponentAccumulator()

    from PixelGeoModelXml.ITkPixelGeoModelConfig import ITkPixelReadoutGeometryCfg
    acc.merge(ITkPixelReadoutGeometryCfg(flags))

    # if we want to do truth seeding, then we need to make the cluster - particle maps beforehand
    ################################################################################
    # Create truth links by associating clusters to hits
    from ActsConfig.ActsTruthConfig import ActsTruthParticleHitCountAlgCfg, ActsPixelClusterToTruthAssociationAlgCfg,ActsStripClusterToTruthAssociationAlgCfg

    if(doTruthSeeding):

        # Schedule Cluster conversion
        from InDetConfig.InDetPrepRawDataFormationConfig import ITkInDetToXAODClusterConversionCfg
        acc.merge(ITkInDetToXAODClusterConversionCfg(flags,
                                                    OutputPixelClustersName="xAODPixelClustersFromInDetCluster",
                                                    OutputStripClustersName="xAODStripClustersFromInDetCluster"))

        acc.merge(ActsPixelClusterToTruthAssociationAlgCfg(flags,
                                                    name="ActsTracccPixelClusterToTruthAssociationAlg",
                                                    InputTruthParticleLinks="xAODTruthLinks",
                                                    AssociationMapOut="ITkTracccPixelClustersToTruthParticles",
                                                    Measurements="xAODPixelClustersFromInDetCluster"))
        acc.merge(ActsStripClusterToTruthAssociationAlgCfg(flags,
                                                        name="ActsTracccStripClusterToTruthAssociationAlg",
                                                        InputTruthParticleLinks="xAODTruthLinks",
                                                        AssociationMapOut="ITkTracccStripClustersToTruthParticles",
                                                        Measurements="xAODStripClustersFromInDetCluster"))

        acc.merge(ActsTruthParticleHitCountAlgCfg(flags,
                                                name="ActsTracccTruthParticleHitCountAlg",
                                                PixelClustersToTruthAssociationMap="ITkTracccPixelClustersToTruthParticles",
                                                StripClustersToTruthAssociationMap="ITkTracccStripClustersToTruthParticles",
                                                TruthParticleHitCountsOut="TracccTruthParticleHitCounts"))

    alg = CompFactory.TrackingAlg("TrackingAlg",
                                  HitInputTool=acc.popToolsAndMerge(TrackingHitInputToolCfg(flags)),
                                  RecoTool=acc.popToolsAndMerge(TrackingRecoToolCfg(flags)),
                                  doTruth=doTruth,
                                  recoFromHits= (not doRecoFromClusters),
                                  makeStandaloneFiles=makeStandaloneFiles,
                                  outputDir=outputDirectory,
                                  recoFromClusters=(doRecoFromClusters or doTruthSeeding),
                                  filesDir=inputDirectory
                                 )

    acc.addEventAlgo(alg)
    ################################################################################
    # Convert ActsTrk::TrackContainer to xAOD::TrackParticleContainer

    from ActsConfig.ActsTrackFindingConfig import ActsTrackToTrackParticleCnvAlgCfg
    acc.merge(ActsTrackToTrackParticleCnvAlgCfg(flags, name="ActsTracccTrackToTrackParticleCnvAlg",
                                                ACTSTracksLocation=["ActsTracccTracks"],
                                                TrackParticlesOutKey="TracccTrackParticles"))

    if(doTruth):

        if(not doTruthSeeding):

            if(doRecoFromClusters):

                acc.merge(ActsPixelClusterToTruthAssociationAlgCfg(flags,
                                                    name="ActsTracccPixelClusterToTruthAssociationAlg",
                                                    InputTruthParticleLinks="xAODTruthLinks",
                                                    AssociationMapOut="ITkTracccPixelClustersToTruthParticles",
                                                    Measurements="xAODPixelClustersFromInDetCluster"))
                acc.merge(ActsStripClusterToTruthAssociationAlgCfg(flags,
                                                                name="ActsTracccStripClusterToTruthAssociationAlg",
                                                                InputTruthParticleLinks="xAODTruthLinks",
                                                                AssociationMapOut="ITkTracccStripClustersToTruthParticles",
                                                                Measurements="xAODStripClustersFromInDetCluster"))

            else:

                acc.merge(ActsPixelClusterToTruthAssociationAlgCfg(flags,
                                                            name="ActsTracccPixelClusterToTruthAssociationAlg",
                                                            InputTruthParticleLinks="xAODTruthLinks",
                                                            AssociationMapOut="ITkTracccPixelClustersToTruthParticles",
                                                            Measurements="xAODPixelClustersFromTracccCluster"))
                acc.merge(ActsStripClusterToTruthAssociationAlgCfg(flags,
                                                                name="ActsTracccStripClusterToTruthAssociationAlg",
                                                                InputTruthParticleLinks="xAODTruthLinks",
                                                                AssociationMapOut="ITkTracccStripClustersToTruthParticles",
                                                                Measurements="xAODStripClustersFromTracccCluster"))
        acc.merge(ActsTruthParticleHitCountAlgCfg(flags,
                                            name="ActsTracccTruthParticleHitCountAlg",
                                            PixelClustersToTruthAssociationMap="ITkTracccPixelClustersToTruthParticles",
                                            StripClustersToTruthAssociationMap="ITkTracccStripClustersToTruthParticles",
                                            TruthParticleHitCountsOut="TracccTruthParticleHitCounts"))

        from ActsConfig.ActsTruthConfig import ActsTrackToTruthAssociationAlgCfg, ActsTrackFindingValidationAlgCfg
        acc.merge(ActsTrackToTruthAssociationAlgCfg(flags,
                                                name="TracccToTruthAssociationAlg",
                                                PixelClustersToTruthAssociationMap="ITkTracccPixelClustersToTruthParticles",
                                                StripClustersToTruthAssociationMap="ITkTracccStripClustersToTruthParticles",
                                                ACTSTracksLocation="ActsTracccTracks",
                                                AssociationMapOut="TracccToTruthParticleAssociation"))

        acc.merge(ActsTrackFindingValidationAlgCfg(flags,
                                                    name="ActsTracccTrackFindingValidationAlg",
                                                    TrackToTruthAssociationMap="TracccToTruthParticleAssociation",
                                                    TruthParticleHitCounts="TracccTruthParticleHitCounts"
                                                    ))


        from ActsConfig.ActsTruthConfig import ActsTrackParticleTruthDecorationAlgCfg
        acc.merge(ActsTrackParticleTruthDecorationAlgCfg(flags, name="ActsTracccTrackParticleTruthDecorationAlg",
                                                        TrackToTruthAssociationMaps=["TracccToTruthParticleAssociation"],
                                                        TrackParticleContainerName="TracccTrackParticles",
                                                        TruthParticleHitCounts="TracccTruthParticleHitCounts",
                                                        ComputeTrackRecoEfficiency=True))

    if(checkSeeds):

        from ActsConfig.ActsSeedingConfig import ActsSeedToTrackCnvAlgCfg
        from ActsConfig.ActsAnalysisConfig import ActsPixelSeedsToTrackParamsAlgCfg
        flags = flags.cloneAndReplace("Tracking.ActiveConfig","Tracking.ITkActsPass")
        acc.merge(ActsPixelSeedsToTrackParamsAlgCfg(flags,
                                                name = "ActsTracccSeedsToTrackParamsAlgCfg",
                                                extension = "Traccc",
                                                InputSeedContainerKey = "ActsTracccSeedTracks",
                                                OutputTrackParamsCollectionKey = "ActsTracccEstimatedTrackParams"))

        acc.merge(ActsSeedToTrackCnvAlgCfg(flags,
                                        name= "ActsTracccSeedToTrackCnvAlg",
                                        EstimatedTrackParametersKey = ["ActsTracccEstimatedTrackParams"],
                                        SeedContainerKey = ["ActsTracccSeedTracks"],
                                        ACTSTracksLocation = "ActsTracccSeedTracks"))


        from ActsConfig.ActsTrackFindingConfig import ActsTrackToTrackParticleCnvAlgCfg
        acc.merge(ActsTrackToTrackParticleCnvAlgCfg(flags, name="ActsTracccSeedTrackToTrackParticleCnvAlg",
                                                    ACTSTracksLocation=["ActsTracccSeedTracks"],
                                                    TrackParticlesOutKey="TracccSeedTrackParticles"))


        if(doTruth):
            from ActsConfig.ActsTruthConfig import ActsTrackToTruthAssociationAlgCfg, ActsTrackFindingValidationAlgCfg
            acc.merge(ActsTrackToTruthAssociationAlgCfg(flags,
                                                    name="TracccSeedToTruthAssociationAlg",
                                                    PixelClustersToTruthAssociationMap="ITkTracccPixelClustersToTruthParticles",
                                                    StripClustersToTruthAssociationMap="ITkTracccStripClustersToTruthParticles",
                                                    ACTSTracksLocation="ActsTracccSeedTracks",
                                                    AssociationMapOut="TracccSeedToTruthParticleAssociation"))

            acc.merge(ActsTrackFindingValidationAlgCfg(flags,
                                                        name="ActsTracccSeedTrackFindingValidationAlg",
                                                        TrackToTruthAssociationMap="TracccSeedToTruthParticleAssociation",
                                                        TruthParticleHitCounts="TracccTruthParticleHitCounts"
                                                        ))

            from ActsConfig.ActsTruthConfig import ActsTrackParticleTruthDecorationAlgCfg
            acc.merge(ActsTrackParticleTruthDecorationAlgCfg(flags, name="ActsTracccSeedTrackParticleTruthDecorationAlg",
                                                            TrackToTruthAssociationMaps=["TracccSeedToTruthParticleAssociation"],
                                                            TrackParticleContainerName="TracccSeedTrackParticles",
                                                            TruthParticleHitCounts="TracccTruthParticleHitCounts",
                                                            ComputeTrackRecoEfficiency=True))

    if writeOutput:
        #Adding the output to the AOD file
        inputList = []

        inputList.append("xAOD::TruthParticleContainer#*")
        inputList.append("xAOD::TruthParticleAuxContainer#*")
        inputList.append("xAOD::TrackJacobianContainer#*")
        inputList.append("xAOD::TrackJacobianAuxContainer#*")
        inputList.append("xAOD::TrackMeasurementContainer#*")
        inputList.append("xAOD::TrackMeasurementAuxContainer#*")
        inputList.append("xAOD::TrackSurfaceContainer#*")
        inputList.append("xAOD::TrackSurfaceAuxContainer#*")
        inputList.append("xAOD::TrackParticleContainer#*")
        inputList.append("xAOD::TrackParticleAuxContainer#*")
        inputList.append("xAOD::PixelClusterContainer#*")
        inputList.append("xAOD::StripClusterContainer#*")
        inputList.append("xAOD::PixelClusterAuxContainer#*")
        inputList.append("xAOD::StripClusterAuxContainer#*")

        acc.merge(OutputStreamCfg(flags, 'AOD', ItemList=inputList))

    return acc

if __name__=="__main__":
    from AthenaCommon.Logging import log
    from AthenaCommon.Constants import INFO
    log.setLevel(INFO)

    from AthenaConfiguration.AllConfigFlags import initConfigFlags
    flags = initConfigFlags()
    flags.Input.Files = ['/cvmfs/atlas-nightlies.cern.ch/repo/data/data-art/PhaseIIUpgrade/RDO/ATLAS-P2-RUN4-03-00-01/mc21_14TeV.601229.PhPy8EG_A14_ttbar_hdamp258p75_SingleLep.recon.RDO.e8514_s4345_r15583_tid39626672_00/RDO.39626672._001121.pool.root.1']
    flags.Input.isMC = True
    flags.Output.AODFileName = "AOD.test.root"

    # flags.IOVDb.GlobalTag="OFLCOND-MC21-SDR-RUN4-03"

    flags.Detector.GeometryITkPixel = True
    flags.Detector.GeometryITkStrip = True
    flags.Detector.GeometryHGTD = False

    flags.Common.MsgSuppression = False

    # flags.Exec.FPE = -1 # FPE check mode: -2 (no FPE check), -1 (abort with core-dump), 0 (FPE Auditor w/o stack-tace) , >0 (number of stack-traces printed by the job)')

    flags.lock()

    from AthenaConfiguration.MainServicesConfig import MainServicesCfg
    cfg = MainServicesCfg(flags)
    from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
    cfg.merge(PoolReadCfg(flags))

    cfg.merge(TrackingAlgCfg(flags))

    cfg.getService('MessageSvc').Format = "%t % F%{:d}W%C%7W%R%T %0W%M".format(flags.Common.MsgSourceLength)

    cfg.run(maxEvents=-1)

# Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration

from AthenaConfiguration.ComponentFactory import CompFactory
from AthenaConfiguration.ComponentAccumulator import ComponentAccumulator
from OutputStreamAthenaPool.OutputStreamConfig import OutputStreamCfg
from AthenaCommon.SystemOfUnits import GeV
import sys
from argparse import ArgumentParser


ComponentAccumulator.debugMode="trackCA trackEventAlgo ..."


def GeoModelConversionAlgCfg(flags):

    acc = ComponentAccumulator()

    alg = CompFactory.GeoModelConversion("GeoModelConversionAlg")

    acc.addEventAlgo(alg)

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

    cfg.merge(GeoModelConversionAlgCfg(flags))

    cfg.run(maxEvents=1)

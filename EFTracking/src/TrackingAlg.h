//
// Copyright (C) 2023-2025 CERN for the benefit of the ATLAS collaboration
//
#ifndef EFTRACKING_TRACKINGALG_H
#define EFTRACKING_TRACKINGALG_H

// Athena/Gaudi includes.
#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/IChronoStatSvc.h"

#include "traccc/geometry/silicon_detector_description.hpp"
#include "vecmem/memory/host_memory_resource.hpp"

#include "EFTracking/ITrackingHitInputTool.h"
#include "EFTracking/ITrackingRecoTool.h"

class TrackingAlg : public AthAlgorithm
{
   public:
      TrackingAlg(const std::string &name, ISvcLocator* pSvcLocator);
      virtual ~TrackingAlg() {};

      /// Initialize the algorithm.
      StatusCode initialize() override;
      /// Execute the algorithm.
      StatusCode execute() override;
      /// Finalize the algorithm.
      StatusCode finalize() override;

   private:
      ToolHandle<ITrackingHitInputTool> m_hitInputTool {
         this,
         "HitInputTool",
         "EFTracking/TrackingHitInputTool",
         "Hit Input Tool"};

      ToolHandle<ITrackingRecoTool> m_recoTool {
         this,
         "RecoTool",
         "EFTracking/TrackingRecoTool",
         "Reco Tool"};

      std::unordered_map<uint64_t, Identifier> m_DetrayToAthenaMap;
      std::unordered_map<Identifier, uint64_t> m_AthenaToDetrayMap;
      void initAthenaDetrayConversionMaps(
           std::unordered_map<uint64_t, Identifier> & detray_to_athena_map
         , std::unordered_map<Identifier, uint64_t> & athena_to_detray_map
      );

      vecmem::host_memory_resource host_mr;
      traccc::silicon_detector_description::host m_dd{host_mr};

      Gaudi::Property<bool> m_makeStandaloneFiles{
         this, "makeStandaloneFiles", false,
         "If you also want to make Traccc standalone files."};
      Gaudi::Property<std::string> m_output_dir{
         this, "outputDir", "",
         "full path to existing location for traccc standalone file writing"};

      Gaudi::Property<bool> m_doRecoFromClusters{
         this, "recoFromClusters", false,
         "Create tracks from clusters"};
      Gaudi::Property<bool> m_doRecoFromHits{
         this, "recoFromHits", true,
         "Create tracks from hits"};
      Gaudi::Property<bool> m_doTruth{
         this, "doTruth", true,
         "Create output containers and link to truth"};
      Gaudi::Property<std::string> m_filesDir{
         this, "filesDir", "/eos/project/a/atlas-eftracking/GPU/ITk_data/",
         "full path to location with detector files"};

      int m_n_event;

      ServiceHandle<IChronoStatSvc> m_chrono{this, "ChronoStatService", "ChronoStatSvc"};

}; // class TrackingAlg

#endif // EFTRACKING_TRACKINGALG_H

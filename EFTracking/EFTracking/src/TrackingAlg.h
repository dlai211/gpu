//
// Copyright (C) 2023 CERN for the benefit of the ATLAS collaboration
//
#ifndef EFTRACKING_TRACKINGALG_H
#define EFTRACKING_TRACKINGALG_H

// Athena/Gaudi includes.
#include "AthenaBaseComps/AthAlgorithm.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"
#include "StoreGate/ReadHandleKey.h"
#include "StoreGate/WriteHandleKey.h"

#include "traccc/cuda/fitting/fitting_algorithm.hpp"
#include "traccc/geometry/silicon_detector_description.hpp"

#include "EFTracking/ITrackingHitInputTool.h"
#include "EFTracking/ITrackingHitMapTool.h"
#include "EFTracking/ITrackingRecoTool.h"
#include "EFTracking/ITrackConversionTool.h"


class TrackingAlg : public AthAlgorithm{

   public:
      TrackingAlg(const std::string &name, ISvcLocator* pSvcLocator);
      virtual ~TrackingAlg() {};

      /// Initialize the algorithm.
      virtual StatusCode initialize() override;
      /// Execute the algorithm.
      virtual StatusCode execute() override;
      /// Finalize the algorithm.
      virtual StatusCode finalize() override;

      /// @}

   private:
   
      void getMaps();
      ToolHandle<ITrackingHitInputTool> m_hitInputTool {this, "HitInputTool", "EFTracking/TrackingHitInputTool", "Hit Input Tool"};
      ToolHandle<ITrackingHitMapTool> m_hitMapTool {this, "HitMapTool", "EFTracking/TrackingHitMapTool", "Hit Mapping Tool"};
      ToolHandle<ITrackingRecoTool> m_recoTool {this, "RecoTool", "EFTracking/TrackingRecoTool", "Reco Tool"};
      ToolHandle<ITrackConversionTool> m_cnvTool {this, "TrackConversionTool", "EFTracking/TrackConversionTool", "Track Conversion Tool"};

      std::map<Identifier, Identifier> m_atlasHumanIDtoIdentifier;
      std::map<std::uint64_t, Identifier> m_DetrayToAtlasMap;
      std::map<Identifier, std::uint64_t> m_AtlasToDetrayMap;
      
      vecmem::host_memory_resource host_mr;
      traccc::silicon_detector_description::host m_dd{host_mr};

}; // class TrackingAlg

#endif // EFTRACKING_TRACKINGALG_H

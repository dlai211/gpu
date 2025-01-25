//
// Copyright (C) 2023 CERN for the benefit of the ATLAS collaboration
//

// Local include(s).
#include "TrackingAlg.h"

// Framework include(s).
#include "AthContainers/tools/copyAuxStoreThinned.h"
#include "StoreGate/ReadHandle.h"
#include "StoreGate/WriteHandle.h"
#include "xAODCore/AuxContainerBase.h"

#include <xAODEventInfo/EventInfo.h>

TrackingAlg::TrackingAlg(const std::string &name, ISvcLocator* pSvcLocator) : AthAlgorithm(name, pSvcLocator){}

StatusCode TrackingAlg::initialize(){

  ATH_CHECK(m_hitInputTool.retrieve());
  ATH_MSG_DEBUG("Initializing Hit Input Tool");

  ATH_CHECK(m_hitMapTool.retrieve());
  ATH_MSG_DEBUG("Initializing Hit Map Tool");

  ATH_CHECK(m_recoTool.retrieve());
  ATH_MSG_DEBUG("Initializing Reco Tool");

  ATH_CHECK(m_cnvTool.retrieve());
  ATH_MSG_DEBUG("Initializing Conversion Tool"); 

  // Return gracefully.
  return StatusCode::SUCCESS;
}

void TrackingAlg::getMaps(){
  
  std::tuple<std::map<std::uint64_t, int>,std::map<int, std::uint64_t>,traccc::silicon_detector_description::host,std::map<int, Identifier>> maps = m_hitMapTool->createMaps();

  m_DetrayToAtlasMap = std::get<0>(maps);
  m_AtlasToDetrayMap = std::get<1>(maps);
  m_dd = std::get<2>(maps);
  m_atlasHumanIDtoIdentifier = std::get<3>(maps);
  
  if(m_atlasHumanIDtoIdentifier.empty() || m_DetrayToAtlasMap.empty() || m_AtlasToDetrayMap.empty()){
    ATH_MSG_DEBUG("Mapping GeoModule<->Detray likely went wrong, will fail gracefully.");
  }
}

StatusCode TrackingAlg::execute(){

  auto context = Gaudi::Hive::currentContext();

  ATH_MSG_DEBUG("Reading event hit information");
  std::vector<clusterInfo> clusters = m_hitInputTool->readClusters(context);
  std::vector<hitInfo> hits = m_hitInputTool->readHits(context);

  if(m_atlasHumanIDtoIdentifier.empty() || m_DetrayToAtlasMap.empty() || m_AtlasToDetrayMap.empty()){
    getMaps();
  }  

  ATH_MSG_DEBUG("Mapping GeoModule<->Detray");
  //std::vector<clusterInfo> detray_clusters = m_hitMapTool->mapClusters(clusters,m_AtlasToDetrayMap);
  //std::vector<hitInfo> detray_hits = m_hitMapTool->mapHits(hits,m_AtlasToDetrayMap);

  ATH_MSG_DEBUG("Passing info to traccc");
  traccc::host_container<traccc::fitting_result<detray::cmath<float> >, traccc::track_state<detray::cmath<float> > > traccc_cluster_tracks = m_recoTool->doRecoFromClusters(clusters);
  traccc::host_container<traccc::fitting_result<detray::cmath<float> >, traccc::track_state<detray::cmath<float> > > traccc_hit_tracks = m_recoTool->doRecoFromHits(hits,clusters);

  ATH_MSG_DEBUG("Passing tracks to conversion tool");
  ATH_CHECK(m_cnvTool->convertTracks(context, traccc_hit_tracks, m_atlasHumanIDtoIdentifier, m_DetrayToAtlasMap));

  return StatusCode::SUCCESS;

};

StatusCode TrackingAlg::finalize(){

  // Return gracefully.
  return StatusCode::SUCCESS;
}


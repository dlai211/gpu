// Copyright (C) 2023-2025 CERN for the benefit of the ATLAS collaboration

// Local include(s).
#include "TrackingAlg.h"
#include "AthenaBaseComps/AthMsgStreamMacros.h"
#include <xAODEventInfo/EventInfo.h>
#include <fstream>

TrackingAlg::TrackingAlg(const std::string &name, ISvcLocator* pSvcLocator)
  : AthAlgorithm(name, pSvcLocator)
{}

StatusCode TrackingAlg::initialize()
{
  ATH_MSG_INFO("TrackingAlg::initialize");

  ATH_MSG_INFO("Retrieving Hit Input Tool");
  ATH_CHECK(m_hitInputTool.retrieve());

  ATH_MSG_INFO("Retrieving Reco Tool");
  ATH_CHECK(m_recoTool.retrieve());

  ATH_MSG_INFO("Retrieved timing tool.");
  ATH_CHECK(m_chrono.retrieve());

  initAthenaDetrayConversionMaps(m_DetrayToAthenaMap, m_AthenaToDetrayMap);
  m_hitInputTool->setAthenaDetrayConversionMaps(m_DetrayToAthenaMap, m_AthenaToDetrayMap);
  m_recoTool->setAthenaDetrayConversionMaps(m_DetrayToAthenaMap, m_AthenaToDetrayMap);

  m_n_event = 0;

  ATH_MSG_INFO("TrackingAlg::initialize complete");
  // Return gracefully.
  return StatusCode::SUCCESS;
}

void TrackingAlg::initAthenaDetrayConversionMaps(
    std::unordered_map<uint64_t, Identifier> & detray_to_athena_map
  , std::unordered_map<Identifier, uint64_t> & athena_to_detray_map
  )
{
  ATH_MSG_DEBUG("constructing athena-detray conversion maps from "
    << m_filesDir + "athenaIdentifierToDetrayMap.txt");
  std::ifstream map_file(m_filesDir + "athenaIdentifierToDetrayMap.txt");
  if (!map_file.is_open()) {
    ATH_MSG_FATAL("error opening athena-detray conversion file");
  }
  std::string str;
  while (std::getline(map_file, str)){
    std::stringstream test(str);
    std::vector<std::string> seglist;
    std::string segment;
    while(std::getline(test, segment, ',')){
      seglist.push_back(segment);
    }

    const std::string& a_id = seglist[0];
    std::string d_id = seglist[1];
    uint64_t idD = std::stoull(d_id);
    Identifier idA;
    idA.set(a_id);
    m_AthenaToDetrayMap[idA] = idD;
    m_DetrayToAthenaMap[idD] = idA;
  }
  ATH_MSG_INFO("athena-detray conversion maps size: " << m_AthenaToDetrayMap.size());
}

StatusCode TrackingAlg::execute()
{
  m_chrono->chronoStart("TrackingAlg::execute");

  auto context = Gaudi::Hive::currentContext();

  if(m_doRecoFromClusters){
    ATH_MSG_INFO("Reading event cluster information");
    m_chrono->chronoStart("traccc read clusters");
    std::pair<std::vector<clusterInfo>,std::map<int,int>> clusters = m_hitInputTool->readClusters(context);
    m_chrono->chronoStop("traccc read clusters");

    ATH_MSG_DEBUG("Passing info to traccc");

    ATH_CHECK(m_recoTool->doRecoFromClusters(context,clusters.first,clusters.second));

    if(m_makeStandaloneFiles){
      ATH_CHECK(m_recoTool->makeTracccStandaloneData(m_n_event,clusters.first,m_output_dir));
    }
  }else if(m_doRecoFromHits){
    ATH_MSG_INFO("Reading event cluster information");
    m_chrono->chronoStart("traccc read clusters");
    std::pair<std::vector<clusterInfo>,std::map<int,int>> clusters = m_hitInputTool->readClusters(context);
    m_chrono->chronoStop("traccc read clusters");

    ATH_MSG_DEBUG("Passing info to traccc");

    if(m_doTruth){ATH_CHECK(m_recoTool->doRecoAndWriteFromHits(context));}
    else{ATH_CHECK(m_recoTool->doRecoFromHits(context));}

    if(m_makeStandaloneFiles){
      ATH_CHECK(m_recoTool->makeTracccStandaloneData(m_n_event,context,m_output_dir));
    }
  }else{
    ATH_MSG_FATAL("You did not specify where you want to start the track reconstruction. Will now fail gracefully.");
  }

  m_n_event++;
  ATH_MSG_DEBUG("Done event: " << m_n_event);

  m_chrono->chronoStop("TrackingAlg::execute");
  return StatusCode::SUCCESS;
};

StatusCode TrackingAlg::finalize()
{
  // Return gracefully.
  return StatusCode::SUCCESS;
}


/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#include "StoreGate/DataHandle.h"

#include "PixelReadoutGeometry/PixelModuleDesign.h"

#include "InDetReadoutGeometry/SiDetectorDesign.h"
#include "InDetReadoutGeometry/SiDetectorElement.h"
#include "ReadoutGeometryBase/SiCellId.h"
#include "ReadoutGeometryBase/SiReadoutCellId.h"
#include "InDetRawData/InDetRawDataCollection.h"
#include "InDetRawData/InDetRawDataContainer.h"
#include "InDetRawData/InDetRawDataCLASS_DEF.h"

#include "TrackingHitInputTool.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>


TrackingHitInputTool::TrackingHitInputTool(const std::string& algname, const std::string& name, const IInterface* ifc) :
  base_class(algname, name, ifc) {}

StatusCode TrackingHitInputTool::initialize() {

  ATH_MSG_DEBUG("Initializing hit reading tool");

  ATH_CHECK(detStore()->retrieve(m_PIX_mgr, "ITkPixel"));
  ATH_CHECK(detStore()->retrieve(m_SCT_mgr, "ITkStrip"));

  ATH_CHECK(m_inputPixelClusterContainerKey.initialize() );
  ATH_CHECK(detStore()->retrieve(m_pixelId, "PixelID"));
  ATH_CHECK(m_pixelRDOKey.initialize());

  ATH_CHECK(m_inputStripClusterContainerKey.initialize() );
  ATH_CHECK(detStore()->retrieve(m_stripId, "SCT_ID"));
  ATH_CHECK(m_stripRDOKey.initialize());


  // Setting up maps
  // std::ifstream map_file("/eos/project/a/atlas-eftracking/GPU/ITk_data/ATLAS-P2-RUN4-03-00-00/athenaIdentifierToDetrayMap.txt");
  std::ifstream map_file("/eos/user/j/jlai/itk_data/old/athenaIdentifierToDetrayMap.txt");

  ATH_MSG_INFO("Constructing athena to detray map");
  std::string str; 
  while (std::getline(map_file, str)){
    
    std::stringstream test(str);
    std::vector<std::string> seglist;
    std::string segment;
    while(std::getline(test, segment, ',')){
      seglist.push_back(segment);
    }

    std::string a_id = seglist[0];
    std::string d_id = seglist[1];
    uint64_t idD = std::stoull(d_id);
    // int idA = std::stoi(a_id);
    Identifier idA;
    idA.set(a_id);
    m_AtlasToDetrayMap[idA] = idD;
    // m_DetrayToAtlasMap[idD] = idA;

  }

  ATH_MSG_INFO("Map size: " << m_AtlasToDetrayMap.size());

  ATH_MSG_DEBUG("Initializing complete");
  return StatusCode::SUCCESS;

}

std::vector<clusterInfo> TrackingHitInputTool::readClusters( const EventContext& eventContext){

  std::vector<clusterInfo> strip_clusters;

  std::vector<clusterInfo> pixel_clusters;
  ATH_MSG_INFO("Reading pixel clusters");
  SG::ReadHandle<InDet::PixelClusterContainer> inputPixelClusterContainer(m_inputPixelClusterContainerKey, eventContext);
  int nPix = 0;
  
  for (const auto *const clusterCollection : *inputPixelClusterContainer) {
    if (!clusterCollection) continue;
    for(const auto *const theCluster : *clusterCollection)  {

      nPix++;
      const InDetDD::SiDetectorElement* element=theCluster->detectorElement();
      const Identifier Pixel_ModuleID = element->identify();
      // if (nPix <= 10) {
      //   ATH_MSG_INFO("Pixel Cluster " << nPix << " center global: " << element->center()[0] << "," << element->center()[1] << "," << element->center()[2]);
      // }

      const std::vector<Identifier>& rdoList = theCluster->rdoList();
      
      clusterInfo thiscluster;
      thiscluster.atlas_id = Pixel_ModuleID;
      uint64_t geometry_id = m_AtlasToDetrayMap[Pixel_ModuleID];
      thiscluster.detray_id = geometry_id;
      thiscluster.globalPosition = theCluster->globalPosition();
      thiscluster.localPosition = theCluster->localPosition();
      thiscluster.local_key = 6;
      auto localCov = theCluster->localCovariance();
      std::vector<float> cov;
      cov.push_back(localCov(0,0));
      cov.push_back(localCov(1,1));
      thiscluster.localCov = cov;
      thiscluster.pixel = true;

      // Print only for the first 10 clusters
      if (nPix <= 10) { 
        ATH_MSG_INFO("Pixel Cluster " << nPix << 
                     " - geometry_id: " << geometry_id << 
                     ", globalPosition: (" << thiscluster.globalPosition.x() << ", " << thiscluster.globalPosition.y() << ", " << thiscluster.globalPosition.z() << ")" << 
                     ", localPosition: (" << thiscluster.localPosition.x() << ", " << thiscluster.localPosition.y() << ")" << 
                     " - localCov(0,0): " << localCov(0,0) << ", localCov(1,1): " << localCov(1,1));
      }

      pixel_clusters.push_back(thiscluster);

    }
  }
  ATH_MSG_INFO("Read "<< nPix << " pixel clusters");

  ATH_MSG_INFO("Reading strip clusters");
  SG::ReadHandle<InDet::SCT_ClusterContainer> inputStripClusterContainer(m_inputStripClusterContainerKey, eventContext);
  
  int nStrip = 0;
  int stripPoz = 0;

  for (const auto *const clusterCollection : *inputStripClusterContainer) {
    if (!clusterCollection) continue;
    for(const auto *const theCluster : *clusterCollection)  {

      ++nStrip;
      const InDetDD::SiDetectorElement *element=theCluster->detectorElement();
      const Identifier StripModuleID = m_stripId->module_id(element->identify());

      if(m_stripId->side(element->identify()) != 0) {continue;} // only side 0 
      ++stripPoz;

      // if (nStrip <= 10) {
      //   ATH_MSG_INFO("Strip cluster " << nStrip << " center: " << element->center()[0] << "," << element->center()[1] << "," << element->center()[2]);
      // }
      // ATH_MSG_DEBUG("Strip Cluster center: " << element->center()[0] << "," << element->center()[1] << "," << element->center()[2]);
      clusterInfo thiscluster;
      thiscluster.atlas_id = StripModuleID;
      uint64_t geometry_id = m_AtlasToDetrayMap[StripModuleID];
      thiscluster.detray_id = geometry_id;
      thiscluster.globalPosition = theCluster->globalPosition();
      thiscluster.localPosition = theCluster->localPosition();
      thiscluster.local_key = 2;
      auto localCov = theCluster->localCovariance();
      std::vector<float> cov;
      cov.push_back(localCov(0,0));
      cov.push_back(localCov(1,1));
      thiscluster.localCov = cov;
      thiscluster.pixel = false;

      // Print only for the first 10 clusters
      if (nStrip <= 10) { 
        ATH_MSG_INFO("Strip Cluster " << nStrip << 
                     " - geometry_id: " << geometry_id << 
                     ", globalPosition: (" << thiscluster.globalPosition.x() << ", " << thiscluster.globalPosition.y() << ", " << thiscluster.globalPosition.z() << ")" << 
                     ", localPosition: (" << thiscluster.localPosition.x() << ", " << thiscluster.localPosition.y() << ")" << 
                    //  ", localPosition<1>: (" << (theCluster->localPosition<1>()).x() << ")" <<
                     " - localCov(0,0): " << localCov(0,0) << ", localCov(1,1): " << localCov(1,1));
      }

      strip_clusters.push_back(thiscluster);

    }
  }

  ATH_MSG_INFO("Read "<< nStrip << " strip clusters and " << stripPoz << " non-skipped (outer side) clusters ");
  // return strip_clusters;
  return strip_clusters;

}

std::vector<hitInfo> TrackingHitInputTool::readHits( const EventContext& eventContext){


  IdContext cntx = m_pixelId->wafer_context();  
  std::vector<hitInfo> pixel_hits;
  std::vector<hitInfo> strip_hits;
  std::vector<hitInfo> all_hits;
  
  ATH_MSG_INFO("Reading pixel hits");
  auto pixelRDOHandle = SG::makeHandle(m_pixelRDOKey, eventContext);
  int nPix = 0;
  for (const InDetRawDataCollection<PixelRDORawData>* pixel_rdoCollection : *pixelRDOHandle) {
    if (pixel_rdoCollection == nullptr) { continue; }
    
    auto pix_rdoColl_id = m_pixelId->show_to_string(pixel_rdoCollection->identify(), &cntx);
    
    for (const PixelRDORawData* pixelRawData : *pixel_rdoCollection) {

      nPix++;
      hitInfo thishit;
      Identifier rdoId = pixelRawData->identify();
      
      // get the det element from the det element collection
      const InDetDD::SiDetectorElement* sielement = m_PIX_mgr->getDetectorElement(rdoId); assert(sielement);
      const Identifier Pixel_ModuleID = sielement->identify();

      Amg::Vector2D LocalPos = sielement->rawLocalPositionOfCell(rdoId);
      thishit.atlas_id = Pixel_ModuleID;
      uint64_t geometry_id = m_AtlasToDetrayMap[Pixel_ModuleID];
      thishit.detray_id = geometry_id;
      thishit.globalPosition = sielement->globalPosition(LocalPos);
      InDetDD::SiCellId id = sielement->cellIdFromIdentifier(rdoId);
      
      thishit.channel0 = id.phiIndex();
      thishit.channel1 = id.etaIndex();
      
      thishit.value = pixelRawData->getToT();
      thishit.timestamp = 8;
      pixel_hits.push_back(thishit);
      all_hits.push_back(thishit);

    }
  }


  ATH_MSG_INFO("Reading strip hits");
  int nStrip = 0;
  constexpr int MaxChannelinStripRow = 128;

  // Write out strip cells information in csv file
  std::string cell_csv_filename = "/eos/user/j/jlai/g200/gpu/G-200/traccc-athena/run/strip_hits.csv";
  std::ofstream cell_csv_file(cell_csv_filename);
  cell_csv_file << "geometry_id,athena_id,stripID,chipID,offset,ITkStripID,idPhiIndex,idEtaIndex,barrel_ec,layer_disk,side,phi_module,eta_module,"
   << "local_x,local_y,global_x,global_y,global_z,center0,center1,center2,thickness,width,length,phiPitch,etaPitch,nFiredStrip,waferID\n";


  auto stripRDOHandle = SG::makeHandle(m_stripRDOKey, eventContext);
  for (const InDetRawDataCollection<SCT_RDORawData>* strip_Collection : *stripRDOHandle) {
    if (strip_Collection == nullptr) { continue; }
    for (const SCT_RDORawData* stripRawData : *strip_Collection) {
      const Identifier rdoId = stripRawData->identify(); // firstStripId
      const InDetDD::SiDetectorElement* sielement = m_SCT_mgr->getDetectorElement(rdoId);
      const Identifier Strip_ModuleID = m_stripId->module_id(sielement->identify());
      Amg::Vector2D localPos = sielement->rawLocalPositionOfCell(rdoId);
      InDetDD::SiCellId id = sielement->cellIdFromIdentifier(rdoId);
      int stripID = m_stripId -> strip(rdoId);
      int side_tmp = m_stripId->side(sielement->identify());

      // Group Size Study
      // if (m_maxFiredStrips != 0u) {
      //   unsigned int nFiredStrips = 0u;
      //   nFiredStrips += stripRawData->getGroupSize();
      //   // if (nFiredStrips > m_maxFiredStrips && nStrip < 10) {
      //   if (nStrip < 10) {
      //     ATH_MSG_INFO("m_maxFiredStrips: " << m_maxFiredStrips << " nFiredStrips: " << nFiredStrips);
      //   }
      // }





      // Each ITK ABC chip reads 128 channels in one row, so we just need to divide the current strip with 128 to get the chip index
      // for the Strip ID, it is the remainder left after dividing by 128
      int chipID = stripID / MaxChannelinStripRow;
      int ITkStripID = stripID % MaxChannelinStripRow;
      int side = m_stripId->side(rdoId);
      unsigned int nFiredStrip = stripRawData->getGroupSize();

      // for each ABC chip readout, each reads 256 channels actually. 0-127 corresponds to lower row and then 128-255 corresponds to the 
      // upper. This can be simulated in the code by using the eta module index. Even index are not offest, while odd index, the 
      // strip id is offest by 128
      // One point to not is that for barrel, the eta module index start at 1, and not zero. Hence a shift of 1 is needed
      int offset = m_stripId->eta_module(rdoId) % 2;
      if(m_stripId->barrel_ec(rdoId) == 0) {
        offset = (std::abs(m_stripId->eta_module(rdoId)) - 1) % 2;
      };

      const IdentifierHash Strip_ModuleHash = m_stripId->wafer_hash(Strip_ModuleID);
      const Identifier Strip_ModuleID2 = m_stripId->wafer_id(Strip_ModuleHash + side);


      ITkStripID += offset * MaxChannelinStripRow;
      uint64_t geometry_id_tmp = m_AtlasToDetrayMap[Strip_ModuleID];
      Amg::Vector3D globalPosition_tmp = sielement->globalPosition(localPos);

      // if (m_AtlasToDetrayMap[Strip_ModuleID] == 1536010436653820607) {
      //   ATH_MSG_INFO("geometry_id: " << m_AtlasToDetrayMap[Strip_ModuleID] << "stripID: " << stripID << " chipID: " << chipID << " offset: " << offset << 
      //               " ITkStripID: " << ITkStripID << " id.strip(): " << id.strip() << " side: " << side_tmp << 
      //               " phiIndex(): " << id.phiIndex() << " etaIndex(): " << id.etaIndex()); }

      cell_csv_file << geometry_id_tmp << ","
          << Strip_ModuleID2 << ","
          << stripID << ","
          << chipID << ","
          << offset << ","
          << ITkStripID << ","
          << id.phiIndex() << ","
          << id.etaIndex() << ","
          << m_stripId->barrel_ec(rdoId) << ","
          << m_stripId->layer_disk(rdoId) << ","
          << m_stripId->side(rdoId) << "," // m_strip->side(sielement->identify()) << ","
          << m_stripId->phi_module(rdoId) << ","
          << m_stripId->eta_module(rdoId) << ","
          << localPos.x() << ","
          << localPos.y() << ","
          << globalPosition_tmp.x() << ","
          << globalPosition_tmp.y() << ","
          << globalPosition_tmp.z() << ","
          << sielement->center()[0] << ","
          << sielement->center()[1] << ","
          << sielement->center()[2] << ","
          << sielement->thickness() << ","
          << sielement->width() << ","
          << sielement->length() << ","
          << sielement->phiPitch() << ","
          << sielement->etaPitch() << ","
          << nFiredStrip << "," 
          << Strip_ModuleHash+side << "\n";


      

      if(m_stripId->side(sielement->identify()) == 1) {continue;} // only side 0




      hitInfo thishit;

      // nStrip++;
      // thishit.atlas_id = Strip_ModuleID;
      // uint64_t geometry_id = m_AtlasToDetrayMap[Strip_ModuleID];
      // thishit.detray_id = geometry_id;
      // thishit.globalPosition = sielement->globalPosition(localPos);
      // thishit.channel0 = id.strip();
      // thishit.channel1 = 0;
      // thishit.value = 1;
      // thishit.timestamp = 8;
      // strip_hits.push_back(thishit);
      // all_hits.push_back(thishit);

        
          
      for (int i = 0; i < nFiredStrip; i++) {
          nStrip++;
          thishit.atlas_id = Strip_ModuleID;
          uint64_t geometry_id = m_AtlasToDetrayMap[Strip_ModuleID];
          thishit.detray_id = geometry_id;
          thishit.globalPosition = sielement->globalPosition(localPos);
          thishit.channel0 = id.strip() + i;
          thishit.channel1 = 0;
          thishit.value = 1;
          thishit.timestamp = 8;
          strip_hits.push_back(thishit);
          all_hits.push_back(thishit);
      }

    }
  }

  ATH_MSG_INFO("Read " << nPix << " pixel hits");
  ATH_MSG_INFO("Read " << nStrip << " strip hits");
  cell_csv_file.close();


  // return pixel_hits;
  return strip_hits;
  

}

StatusCode TrackingHitInputTool::finalize() {
  return StatusCode::SUCCESS;
}




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
#include "ReadoutGeometryBase/SiCellId.h"
#include "ReadoutGeometryBase/SiReadoutCellId.h"

#include "TrackingHitInputTool.h"
#include <iostream>
#include <sstream>
#include <fstream>


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

  ATH_MSG_DEBUG("Initializing complete");
  return StatusCode::SUCCESS;

}

std::vector<clusterInfo> TrackingHitInputTool::readClusters( const EventContext& eventContext){

  
  std::vector<clusterInfo> clusters;
  ATH_MSG_INFO("Reading pixel clusters");
  SG::ReadHandle<InDet::PixelClusterContainer> inputPixelClusterContainer(m_inputPixelClusterContainerKey, eventContext);
  int nPix = 0;
  
  for (const auto *const clusterCollection : *inputPixelClusterContainer) {
    if (!clusterCollection) continue;
    for(const auto *const theCluster : *clusterCollection)  {

      nPix++;
      const InDetDD::SiDetectorElement* element=theCluster->detectorElement();
      const Identifier Pixel_ModuleID = element->identify();
      
      int humanReadableID = 1000000*1 /*1 = Pixel*/
      + 100000*m_pixelId->layer_disk(Pixel_ModuleID)
      + 1000*(10+m_pixelId->eta_module(Pixel_ModuleID))
      + m_pixelId->phi_module(Pixel_ModuleID);
      if ( m_pixelId->barrel_ec(Pixel_ModuleID) != 0 ) {
          humanReadableID = m_pixelId->barrel_ec(Pixel_ModuleID)*(humanReadableID + 10000000);
      }

      ATH_MSG_DEBUG("Cluster center global: " << element->center()[0] << "," << element->center()[1] << "," << element->center()[2]);

      const std::vector<Identifier>& rdoList = theCluster->rdoList();
      if(humanReadableID == -22020000){
        
        ATH_MSG_DEBUG("Module ID: " << humanReadableID);
        ATH_MSG_DEBUG("Cluster center: " << theCluster->globalPosition()[0] << "," << theCluster->globalPosition()[1] << "," << theCluster->globalPosition()[2]);
        ATH_MSG_DEBUG("Local coordinates: " << theCluster->localPosition()[0] << "," << theCluster->localPosition()[1]);
        ATH_MSG_DEBUG("Hits associated to this cluster: ");
	      std::vector<Identifier>::const_iterator nextRDO; 
        int p = 0;
        for(nextRDO=rdoList.begin(); nextRDO !=rdoList.end(); ++nextRDO){ 
          Identifier rdoId = (*nextRDO); 
                ATH_MSG_DEBUG("hit n. " << p << " : " << m_pixelId->phi_index(rdoId) << " / " << m_pixelId->eta_index(rdoId));
                p++;
        }
        
      }
      clusterInfo thiscluster;
      thiscluster.atlas_id = humanReadableID;
      thiscluster.globalPosition = theCluster->globalPosition();
      thiscluster.localPosition = theCluster->localPosition();
      thiscluster.local_key = 6;
      auto localCov = theCluster->localCovariance();
      std::vector<float> cov;
      cov.push_back(localCov(0,0));
      cov.push_back(localCov(1,1));
      thiscluster.localCov = cov;
      thiscluster.pixel = true;
      clusters.push_back(thiscluster);

    }
  }
  ATH_MSG_INFO("Read "<< nPix << " pixel clusters");

  SG::ReadHandle<InDet::SCT_ClusterContainer> inputStripClusterContainer(m_inputStripClusterContainerKey, eventContext);
  int nStrip = 0;
  for (const auto *const clusterCollection : *inputStripClusterContainer) {
    if (!clusterCollection) continue;
    for(const auto *const theCluster : *clusterCollection)  {

      ++nStrip;
      const InDetDD::SiDetectorElement *element=theCluster->detectorElement();
      const Identifier StripModuleID = m_stripId->module_id(element->identify());

      if(m_stripId->side(element->identify()) != 0) {continue;} // not inner side case 

      int humanReadableID = 1000000*2 /*2 = Strip*/
      + 100000*m_stripId->layer_disk(StripModuleID)
      + 1000*(10+m_stripId->eta_module(StripModuleID))
      + m_stripId->phi_module(StripModuleID);
      if ( m_stripId->barrel_ec(StripModuleID) != 0 ) {
          humanReadableID = m_stripId->barrel_ec(StripModuleID)*(humanReadableID + 10000000);
      }

      ATH_MSG_DEBUG("Cluster center: " << element->center()[0] << "," << element->center()[1] << "," << element->center()[2]);
      clusterInfo thiscluster;
      thiscluster.atlas_id = humanReadableID;
      thiscluster.globalPosition = theCluster->globalPosition();
      thiscluster.localPosition = theCluster->localPosition();
      thiscluster.local_key = 2;
      auto localCov = theCluster->localCovariance();
      std::vector<float> cov;
      cov.push_back(localCov(0,0));
      cov.push_back(localCov(1,1));
      thiscluster.localCov = cov;
      thiscluster.pixel = false;
      clusters.push_back(thiscluster);

    }
  }

  ATH_MSG_INFO("Read "<< nStrip << " strip clusters");
  return clusters;

}

std::vector<hitInfo> TrackingHitInputTool::readHits( const EventContext& eventContext){

  IdContext cntx = m_pixelId->wafer_context();  
  std::vector<hitInfo> hits;
  
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

      int humanReadableID = 1000000*1 /*1 = Pixel*/
      + 100000*m_pixelId->layer_disk(Pixel_ModuleID)
      + 1000*(10+m_pixelId->eta_module(Pixel_ModuleID))
      + m_pixelId->phi_module(Pixel_ModuleID);
      if ( m_pixelId->barrel_ec(Pixel_ModuleID) != 0 ) {
          humanReadableID = m_pixelId->barrel_ec(Pixel_ModuleID)*(humanReadableID + 10000000);
      }

      Amg::Vector2D LocalPos = sielement->rawLocalPositionOfCell(rdoId);
      thishit.atlas_id = humanReadableID;
      thishit.globalPosition = sielement->globalPosition(LocalPos);
      InDetDD::SiCellId id = sielement->cellIdFromIdentifier(rdoId);
      
      thishit.channel0 = id.phiIndex();
      thishit.channel1 = id.etaIndex();
      
      thishit.value = pixelRawData->getToT();
      thishit.timestamp = 8;
      hits.push_back(thishit);

    }
  }

  return hits;

  int nStrip = 0;
  auto stripRDOHandle = SG::makeHandle(m_stripRDOKey, eventContext);
  for (const InDetRawDataCollection<SCT_RDORawData>* strip_Collection : *stripRDOHandle) {
    if (strip_Collection == nullptr) { continue; }
    for (const SCT_RDORawData* stripRawData : *strip_Collection) {
      const Identifier rdoId = stripRawData->identify();
      const InDetDD::SiDetectorElement* sielement = m_SCT_mgr->getDetectorElement(rdoId);
      const Identifier Strip_ModuleID = m_stripId->module_id(sielement->identify());
      Amg::Vector2D localPos = sielement->rawLocalPositionOfCell(rdoId);
      InDetDD::SiCellId id = sielement->cellIdFromIdentifier(rdoId);

      nStrip++;

      int humanReadableID = 1000000*2 //2 = Strip
      + 100000*m_stripId->layer_disk(Strip_ModuleID)
      + 1000*(10+m_stripId->eta_module(Strip_ModuleID))
      + m_stripId->phi_module(Strip_ModuleID);
      if ( m_stripId->barrel_ec(Strip_ModuleID) != 0 ) {
        humanReadableID = m_stripId->barrel_ec(Strip_ModuleID)*(humanReadableID + 10000000);
      }

      if(humanReadableID == 22464004){ATH_MSG_INFO("I am a strip module!");}

      hitInfo thishit;
      thishit.atlas_id = humanReadableID;
      thishit.globalPosition = sielement->globalPosition(localPos);
      thishit.channel0 = id.phiIndex();
      
      // do i need this?
      //thishit.value = pixelRawData->getToT();
      thishit.timestamp = 8;
      hits.push_back(thishit);

    }
  }

  ATH_MSG_INFO("Read " << nPix << " pixel hits");
  ATH_MSG_INFO("Read " << nStrip << " strip hits");

  //return hits;
  

}

StatusCode TrackingHitInputTool::finalize() {
  return StatusCode::SUCCESS;
}



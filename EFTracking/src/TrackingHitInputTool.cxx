
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

#include "SCT_ReadoutGeometry/SCT_ModuleSideDesign.h"
#include "SCT_ReadoutGeometry/StripStereoAnnulusDesign.h"
#include "SCT_ReadoutGeometry/SCT_BarrelModuleSideDesign.h"
#include "SCT_ReadoutGeometry/SCT_ForwardModuleSideDesign.h"
#include "xAODInDetMeasurement/PixelClusterAuxContainer.h"
#include "xAODInDetMeasurement/StripClusterAuxContainer.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

// Athena StripClusteringTool.cxx
#include <InDetRawData/SCT3_RawData.h>
#include <SCT_ReadoutGeometry/SCT_ModuleSideDesign.h>
#include <SCT_ReadoutGeometry/StripStereoAnnulusDesign.h>
#include <TrkSurfaces/Surface.h>
#include <algorithm>
#include <stdexcept>

// SiLocalPosition.h
#include "ReadoutGeometryBase/SiLocalPosition.h"


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

  ATH_CHECK(m_lorentzAngleTool.retrieve());


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
  cell_csv_file << std::fixed << std::setprecision(10);
  cell_csv_file << "geometry_id,athena_id,stripID,chipID,offset,ITkStripID,idPhiIndex,idEtaIndex,barrel_ec,layer_disk,side,phi_module,eta_module,"
   << "local_x,local_y,global_x,global_y,global_z,center0,center1,center2,thickness,width,length,phiPitch,etaPitch,nFiredStrip,waferID,lorentzShift,angleShift,"
   << "strip,row,xPhi,xEta\n";

  std::string annulus_test_filename = "/eos/user/j/jlai/g200/gpu/G-200/traccc-athena/run/annulus_test.csv";
  std::ofstream annulus_test_file(annulus_test_filename);
  annulus_test_file << std::fixed << std::setprecision(10);
  annulus_test_file << "id,rdoId,geometry_id,athena_id,local_x,local_y,side,strip,row,"
                    << "x_beam,y_beam,rad_beam,phi_beam,x_strip,y_strip,rad_strip,phi_strip,"
                    << "stereo,waferCentreR,centreR,phiPitch,minR,maxR,lengthBF,m_nStrips0,m_pitch0,type\n";

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

      double shift = m_lorentzAngleTool->getLorentzShift(Strip_ModuleHash+side, eventContext);
      double angle_shift = m_lorentzAngleTool->getTanLorentzAngle(Strip_ModuleHash+side, eventContext);

      ITkStripID += offset * MaxChannelinStripRow;
      uint64_t geometry_id_tmp = m_AtlasToDetrayMap[Strip_ModuleID];
      Amg::Vector3D globalPosition_tmp = sielement->globalPosition(localPos);
      
      double localpos_x = localPos.x();
      double localpos_y = localPos.y();

      InDetDD::SiLocalPosition pos_tmp;
      int strip = -1, row = -1;
      // double stereo = -1., waferCentreR = -1., centreR = -1., phiPitch = -1., minR = -1., maxR = -1., lengthBF = -1.;

      if (m_stripId->barrel_ec(rdoId) != 0) { // Annulus
        const InDetDD::StripStereoAnnulusDesign* annulus_design = dynamic_cast<const InDetDD::StripStereoAnnulusDesign*>(&sielement->design());
        if (!annulus_design) {
          ATH_MSG_WARNING("Failed to cast to StripStereoAnnulusDesign");
          continue;
        }

        auto [strip_tmp, row_tmp] = annulus_design->getStripRow(id);
        strip = strip_tmp;
        row = row_tmp;


        pos_tmp = annulus_design->localPositionOfCellPC(id);

        // Save one module information for annulus
        if (nStrip == 0) {
          for (int stripIndex = 0; stripIndex < 1024; ++stripIndex) {
            InDetDD::SiCellId id_tmp(stripIndex, 0);

            Identifier rdoId_tmp = sielement->identifierFromCellId(id_tmp);
            InDetDD::SiDetectorElement* sielement_tmp = m_SCT_mgr->getDetectorElement(rdoId_tmp);
            Amg::Vector2D localPos_tmp = sielement_tmp->rawLocalPositionOfCell(rdoId_tmp); // beam frame {x_beam, y_beam}
            int side_tmp = m_stripId->side(rdoId_tmp);

            auto [strip, row] = annulus_design->getStripRow(id_tmp);
            InDetDD::SiLocalPosition posxy_beam = annulus_design->localPositionOfCell(id_tmp); // beam frame {x_beam, y_beam}
            InDetDD::SiLocalPosition posrphi_beam = annulus_design->localPositionOfCellPC(id_tmp); // beam frame {xEta=>rPrime, xPhi=>phi}
            InDetDD::SiLocalPosition posxy_strip = annulus_design->beamToStrip(posxy_beam); // strip frame {x_strip, y_strip}
            InDetDD::SiLocalPosition posrphi_strip = annulus_design->beamToStripPC(posrphi_beam); // strip frame {rad_strip, phi_strip}

            int m_nStrips0 = annulus_design->diodesInRow(row);
            double m_pitch0 = annulus_design->phiPitchPhi(id_tmp);
            InDetDD::DetectorType type = annulus_design->type();

            double stereo = annulus_design->stereo();
            double waferCentreR = annulus_design->waferCentreR();
            double centreR = annulus_design->centreR();
            double phiPitch = annulus_design->phiPitch(id_tmp);
            double minR = annulus_design->minR();
            double maxR = annulus_design->maxR();  
            double lengthBF = 2. * waferCentreR * std::sin(stereo / 2.);

            annulus_test_file << id_tmp << "," << rdoId_tmp << "," << geometry_id_tmp << "," 
            << Strip_ModuleID2 << "," << localPos_tmp.x() << "," << localPos_tmp.y() << ","
            << side_tmp << "," << strip << "," << row << "," 
            << posxy_beam.xEta() << "," 
            << posxy_beam.xPhi() << ","
            << posrphi_beam.xEta() << ","
            << posrphi_beam.xPhi() << ","
            << posxy_strip.xEta() << ","
            << posxy_strip.xPhi() << ","
            << posrphi_strip.xEta() << ","
            << posrphi_strip.xPhi() << ","
            << stereo << ","
            << waferCentreR << ","
            << centreR << ","
            << phiPitch << ","
            << minR << ","
            << maxR << ","
            << lengthBF << ","
            << m_nStrips0 << ","
            << m_pitch0 << ","
            << type << "\n";
          }
        }
      }

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
          << localpos_x << ","
          << localpos_y << ","
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
          << Strip_ModuleHash+side << "," 
          << shift << ","
          << angle_shift << ","
          << strip << ","
          << row << ","
          << pos_tmp.xPhi() << ","
          << pos_tmp.xEta() << "\n";


      

      if(m_stripId->side(sielement->identify()) == 1) {continue;} // only side 0

      // Annulus 
      


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


  // ATH_MSG_INFO("Start performing Lorentz Shift");
  // std::string measurement_filename = "/eos/user/j/jlai/Traccc/traccc/data/output/measurements.csv";
  // std::ifstream file(measurement_filename);
  // std::string line; 

  // std::getline(file, line); // Skip the header line
  // std::vector<float> local0_values;
  // std::vector<float> local1_values;

  // while (std::getline(file, line)) {
  //   std::stringstream ss(line);
  //   std::string token;
  //   int col_idx = 0;
  //   float local0 = 0.0, local1 = 0.0;

  //   while (std::getline(ss, token, '\t')) {
  //     if (col_idx == 0) local0 = std::stof(token);
  //     if (col_idx == 1) local1 = std::stof(token);
  //     col_idx++;
  //   }

  //   local0_values.push_back(local0);
  //   local1_values.push_back(local1);
  // }

  // // Print the first 5 values of local0 and local1
  // for (size_t i = 0; i < std::min(local0_values.size(), size_t(5)); ++i) {
  //       std::cout << "measurement.csv: "
  //                 << "  local0: " << local0_values[i]
  //                 << ", local1: " << local1_values[i] << std::endl;
  // }

  // Shift
  // double shift = m_lorentzAngleTool->getLorentzAngle(Strip_ModuleHash, eventContext);


  // return pixel_hits;
  return strip_hits;
  

}

StatusCode TrackingHitInputTool::finalize() {
  return StatusCode::SUCCESS;
}



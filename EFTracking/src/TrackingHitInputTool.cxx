
/*
  Copyright (C) 2002-2025 CERN for the benefit of the ATLAS collaboration
*/

#include "TrackingHitInputTool.h"

#include "PixelReadoutGeometry/PixelModuleDesign.h"

#include "InDetReadoutGeometry/SiDetectorDesign.h"
#include "InDetReadoutGeometry/SiDetectorElement.h"
#include "ReadoutGeometryBase/SiCellId.h"
#include "ReadoutGeometryBase/SiCellId.h"

#include "PixelReadoutGeometry/PixelModuleDesign.h"
#include "SCT_ReadoutGeometry/SCT_ModuleSideDesign.h"
#include "SCT_ReadoutGeometry/StripStereoAnnulusDesign.h"
#include "SCT_ReadoutGeometry/SCT_BarrelModuleSideDesign.h"
#include "SCT_ReadoutGeometry/SCT_ForwardModuleSideDesign.h"
#include "xAODInDetMeasurement/PixelClusterAuxContainer.h"
#include "xAODInDetMeasurement/StripClusterAuxContainer.h"
#include <cstdint>
#include "AthenaDetrayConversion.h"
#include <fstream>

TrackingHitInputTool::TrackingHitInputTool(const std::string& algname, const std::string& name, const IInterface* ifc)
  : base_class(algname, name, ifc)
  , m_AthenaToDetrayMap(nullptr)
  , m_DetrayToAthenaMap(nullptr)
{}

StatusCode TrackingHitInputTool::initialize()
{
  ATH_MSG_INFO("TrackingHitInputTool::initialize");

  ATH_CHECK(detStore()->retrieve(m_pixelManager));
  ATH_CHECK(detStore()->retrieve(m_stripManager));

  ATH_CHECK(m_inputPixelClusterContainerKey.initialize() );
  ATH_CHECK(detStore()->retrieve(m_pixelID, "PixelID"));
  ATH_CHECK(m_pixelRDOKey.initialize());

  ATH_CHECK(m_inputStripClusterContainerKey.initialize() );
  ATH_CHECK(detStore()->retrieve(m_stripID, "SCT_ID"));
  ATH_CHECK(m_stripRDOKey.initialize());

  ATH_CHECK(m_xAODPixelClusterFromTracccClusterKey.initialize());
  ATH_CHECK(m_xAODStripClusterFromTracccClusterKey.initialize());
  ATH_CHECK(m_xAODSpacepointFromTracccClusterKey.initialize());

  ATH_CHECK(m_xAODPixelClusterFromInDetClusterKey.initialize());
  ATH_CHECK(m_xAODStripClusterFromInDetClusterKey.initialize());
  ATH_CHECK(m_xAODSpacepointFromInDetClusterKey.initialize());

  ATH_CHECK( m_pixelDetEleCollKey.initialize() );
  ATH_CHECK( m_stripDetEleCollKey.initialize() );

  ATH_CHECK(m_lorentzAngleTool.retrieve());

  std::ifstream acts_map_file(m_filesDir+"actsToAthenaIdentifierMap.txt");
  std::string str;
  while (std::getline(acts_map_file, str)){

    std::stringstream test(str);
    std::vector<std::string> seglist;
    std::string segment;
    while(std::getline(test, segment, ',')){
      seglist.push_back(segment);
    }

    const std::string& acts_id = seglist[0];
    std::string athena_id = seglist[1];
    Identifier idA;
    idA.set(athena_id);
    m_AthenaToActsMap[idA] = acts_id;

  }
  ATH_MSG_INFO("Acts to Athena conversion map: " << m_AthenaToActsMap.size());

  ATH_MSG_INFO("TrackingHitInputTool::initialize complete");
  return StatusCode::SUCCESS;

}

void TrackingHitInputTool::setAthenaDetrayConversionMaps(
  std::unordered_map<uint64_t, Identifier> const & detray_to_athena_map,
  std::unordered_map<Identifier, uint64_t> const & athena_to_detray_map
  )
{
  m_AthenaToDetrayMap = &athena_to_detray_map;
  m_DetrayToAthenaMap = &detray_to_athena_map;
}

std::pair<std::vector<clusterInfo>,std::map<int,int>>
  TrackingHitInputTool::readClusters(const EventContext& eventContext)
{
  std::vector<clusterInfo> clusters;
  int nPix = 0;
  int nStrip = 0;
  m_tracccMeasToxAODClusterMap.clear();

  ATH_MSG_INFO("Converting InDet clusters to xAOD");
  SG::WriteHandle<xAOD::PixelClusterContainer> xAODPixelContainerFromInDetClusters(m_xAODPixelClusterFromInDetClusterKey, eventContext);
  if( (xAODPixelContainerFromInDetClusters.record(std::make_unique<xAOD::PixelClusterContainer>(),
            std::make_unique<xAOD::PixelClusterAuxContainer>())).isFailure() ){
              ATH_MSG_FATAL("Could not record xAOD Pixel container from InDet Clusters");
  }

  SG::WriteHandle<xAOD::StripClusterContainer> xAODStripContainerFromInDetClusters(m_xAODStripClusterFromInDetClusterKey, eventContext);
  if( (xAODStripContainerFromInDetClusters.record(std::make_unique<xAOD::StripClusterContainer>(),
            std::make_unique<xAOD::StripClusterAuxContainer>())).isFailure() ){
              ATH_MSG_FATAL("Could not record xAOD Strip container from InDet Clusters");
  }

  SG::ReadCondHandle<InDetDD::SiDetectorElementCollection> pixelDetEleHandle(m_pixelDetEleCollKey, eventContext);
  const InDetDD::SiDetectorElementCollection* pixElements(*pixelDetEleHandle);

  if (not pixelDetEleHandle.isValid() or pixElements==nullptr) {
    ATH_MSG_FATAL(m_pixelDetEleCollKey.fullKey() << " is not available.");
  }

  SG::ReadCondHandle<InDetDD::SiDetectorElementCollection> stripDetEleHandle(m_stripDetEleCollKey, eventContext);
  const InDetDD::SiDetectorElementCollection* stripElements(*stripDetEleHandle);
  if (not stripDetEleHandle.isValid() or stripElements==nullptr) {
    ATH_MSG_FATAL(m_stripDetEleCollKey.fullKey() << " is not available.");
  }

  SG::WriteHandle<xAOD::SpacePointContainer> xAODSpacepointContainerFromInDetClusters(m_xAODSpacepointFromInDetClusterKey, eventContext);
  if( xAODSpacepointContainerFromInDetClusters.record(std::make_unique<xAOD::SpacePointContainer>(),
						 std::make_unique<xAOD::SpacePointAuxContainer>()).isFailure() ){
              throw std::runtime_error("creation of InDet spacepoint containers failed");
  }

  ATH_MSG_INFO("Reading clusters");
  SG::ReadHandle<InDet::PixelClusterContainer> inputPixelClusterContainer(m_inputPixelClusterContainerKey, eventContext);

  for (const auto *const clusterCollection : *inputPixelClusterContainer) {
    if (!clusterCollection) continue;
    for(const auto *const theCluster : *clusterCollection)  {

      const InDetDD::SiDetectorElement* element=theCluster->detectorElement();
      const Identifier Pixel_ModuleID = element->identify();

      ATH_MSG_DEBUG("Cluster center global: " << element->center()[0] << "," << element->center()[1] << "," << element->center()[2]);

      clusterInfo thiscluster;
      thiscluster.atlas_id = Pixel_ModuleID;
      uint64_t const geometry_id = findDetrayID(Pixel_ModuleID, m_AthenaToDetrayMap);
      thiscluster.detray_id = geometry_id;
      thiscluster.globalPosition = theCluster->globalPosition();
      thiscluster.localPosition = theCluster->localPosition();
      thiscluster.local_key = 6;
      thiscluster.measurement_id = nPix+nStrip;

      thiscluster.pixel = true;
      clusters.push_back(thiscluster);

      xAOD::PixelCluster * pixelCl = new xAOD::PixelCluster();
      xAODPixelContainerFromInDetClusters->push_back(pixelCl);
      if((convertInDetToXaodCluster(*theCluster, *element, *pixelCl)).isFailure()){
        ATH_MSG_FATAL("Could not convert InDet pixel cluster to xAOD");
      }

      xAOD::SpacePoint *xaod_sp = new xAOD::SpacePoint();
      xAODSpacepointContainerFromInDetClusters->push_back(xaod_sp);
      const IdentifierHash Pixel_ModuleHash = m_pixelID->wafer_hash(Pixel_ModuleID);
      Eigen::Matrix<float,2,1> globalVariance(0, 0);
      Eigen::Matrix<float,3,1> globalPosition((theCluster->globalPosition()).x(), (theCluster->globalPosition()).y(), (theCluster->globalPosition()).z());
      xaod_sp->setSpacePoint(Pixel_ModuleHash, globalPosition, globalVariance(0,0), globalVariance(1,0), {pixelCl});

      size_t index = xAODPixelContainerFromInDetClusters->size() - 1;
      m_tracccMeasToxAODClusterMap[nPix+nStrip] = index;

      nPix++;
    }
  }
  ATH_MSG_INFO("Read "<< nPix << " pixel clusters");

  SG::ReadHandle<InDet::SCT_ClusterContainer> inputStripClusterContainer(m_inputStripClusterContainerKey, eventContext);

  // for debugging purposes,
  // printing all clusters and the associated cells on one endcap strip module
  bool printClusterInfo = false;

  // Save the strip clusters output to csv file
  std::string athena_strip_filename = "strip_cluster_athena.csv";
  std::ofstream athena_strip_file(athena_strip_filename);
  athena_strip_file << std::fixed << std::setprecision(10);
  athena_strip_file << "geometry_id,athena_id,barrel_ec,local_x,local_y,global_x,global_y,global_z\n";

  for (const auto *const clusterCollection : *inputStripClusterContainer) {
    if (!clusterCollection) continue;
    for(const auto *const theCluster : *clusterCollection)  {

      const InDetDD::SiDetectorElement *element=theCluster->detectorElement();
      const Identifier strip_moduleID = m_stripID->module_id(element->identify()); //from wafer id to module id
      const IdentifierHash Strip_ModuleHash = m_stripID->wafer_hash(strip_moduleID);

      // Extract the correct Strip_ModuleID
      int side = m_stripID->side(element->identify());
      const Identifier Strip_ModuleID = m_stripID->wafer_id(Strip_ModuleHash + side);

      if(Strip_ModuleID.get_compact() == 0xd4e418000000000){printClusterInfo = true;}

      if(printClusterInfo){
        const std::vector<Identifier>& rdoList = theCluster->rdoList();

        for (auto rdoIter : rdoList) {
          const InDetDD::SiCellId& chargeCellId =  element->cellIdFromIdentifier(rdoIter);
          if(printClusterInfo){ATH_MSG_DEBUG("hit id at this cluster: " << chargeCellId);}
        }

      }

      ATH_MSG_DEBUG("Cluster center: " << element->center()[0] << "," << element->center()[1] << "," << element->center()[2]);
      clusterInfo thiscluster;
      thiscluster.atlas_id = Strip_ModuleID;
      uint64_t const geometry_id = findDetrayID(Strip_ModuleID, m_AthenaToDetrayMap);
      thiscluster.detray_id = geometry_id;
      thiscluster.globalPosition = theCluster->globalPosition();
      thiscluster.localPosition = theCluster->localPosition();
      if(printClusterInfo){ATH_MSG_DEBUG("local pos of this cluster: " << theCluster->localPosition().x() << "," << theCluster->localPosition().y());}
      thiscluster.local_key = 2;
      thiscluster.measurement_id = nPix+nStrip;

      if (m_stripID->barrel_ec(Strip_ModuleID) == 0) {

        auto boundsType = element->bounds().type();
        if(boundsType == Trk::SurfaceBounds::Annulus){
          thiscluster.local_key = 4;
        }
      }else{

        auto boundsType = element->bounds().type();
        if(boundsType == Trk::SurfaceBounds::Annulus){
          thiscluster.local_key = 4;
        }
      }

      thiscluster.pixel = false;
      clusters.push_back(thiscluster);

      xAOD::StripCluster * stripCl = new xAOD::StripCluster();
      xAODStripContainerFromInDetClusters->push_back(stripCl);
      if((convertInDetToXaodCluster(*theCluster, *element, *stripCl)).isFailure()){
        ATH_MSG_FATAL("Could not convert InDet strip cluster to xAOD");
      }
      size_t index = xAODStripContainerFromInDetClusters->size() - 1;
      m_tracccMeasToxAODClusterMap[nPix+nStrip] = index;
      if(printClusterInfo){ATH_MSG_DEBUG("Done cluster: " << index);}

      nStrip++;

      athena_strip_file << geometry_id << ","
        << Strip_ModuleID << ","
        << m_stripID->barrel_ec(Strip_ModuleID) << ","
        << theCluster->localPosition().x() << ","
        << theCluster->localPosition().y() << ","
        << theCluster->globalPosition().x() << ","
        << theCluster->globalPosition().y() << ","
        << theCluster->globalPosition().z() << "\n";
    }
  }

  athena_strip_file.close();

  ATH_MSG_INFO("Read "<< nStrip << " strip clusters");
  ATH_MSG_INFO("Combined we have " << clusters.size() << " clusters.");

  return std::make_pair(clusters,m_tracccMeasToxAODClusterMap);

}

std::map<int,int> TrackingHitInputTool::convertClusters(const EventContext& eventContext,traccc::edm::silicon_cluster_collection::host& traccc_clusters, traccc::measurement_collection_types::host& traccc_measurements, traccc::edm::silicon_cell_collection::host& traccc_cells){

  SG::WriteHandle<xAOD::PixelClusterContainer> PixelContainerFromTracccClusters(m_xAODPixelClusterFromTracccClusterKey, eventContext);
  if( PixelContainerFromTracccClusters.record(std::make_unique<xAOD::PixelClusterContainer>(),
						 std::make_unique<xAOD::PixelClusterAuxContainer>()).isFailure() ){
              throw std::runtime_error("creation of traccc pixel cluster containers failed");
  }

  SG::WriteHandle<xAOD::StripClusterContainer> StripContainerFromTracccClusters(m_xAODStripClusterFromTracccClusterKey, eventContext);
  if( StripContainerFromTracccClusters.record(std::make_unique<xAOD::StripClusterContainer>(),
						 std::make_unique<xAOD::StripClusterAuxContainer>()).isFailure() ){
              throw std::runtime_error("creation of traccc strip cluster containers failed");
  }

  SG::WriteHandle<xAOD::SpacePointContainer> SpacepointContainerFromTracccClusters(m_xAODSpacepointFromTracccClusterKey, eventContext);
  if( SpacepointContainerFromTracccClusters.record(std::make_unique<xAOD::SpacePointContainer>(),
						 std::make_unique<xAOD::SpacePointAuxContainer>()).isFailure() ){
              throw std::runtime_error("creation of traccc spacepoint containers failed");
  }

  int nPix = 0;
  int nStrip = 0;
  m_tracccMeasToxAODClusterMap.clear();

  // for debugging purposes,
  // printing all clusters and the associated cells on one endcap strip module
  bool printClusterInfo = false;

  ATH_MSG_INFO("Converting clusters from Traccc");
  // Save the strip clusters output to csv file
  std::string traccc_strip_filename = "strip_cluster_traccc.csv";
  std::ofstream traccc_strip_file(traccc_strip_filename);
  traccc_strip_file << std::fixed << std::setprecision(10);
  traccc_strip_file << "athena_id,barrel_ec,local_x,local_y,global_x,global_y,global_z\n";

  for(size_t i = 0; i < traccc_clusters.size(); i++){

    const auto& cluster = traccc_clusters[i];
    const auto& meas = traccc_measurements[i];
    uint64_t detray_id = meas.surface_link.value();

    Identifier athena_ID = findAthenaID(detray_id, m_DetrayToAthenaMap);
    if(athena_ID.get_compact() == 0xd4e418000000000){printClusterInfo = true;}

    std::vector<Identifier> rdoList;

    int colMin = 0;
    int colMax = 0;
    int rowMin = 0;
    int rowMax = 0;



    if (meas.meas_dim == 2u) {

      xAOD::PixelCluster *xaod_pcl = new xAOD::PixelCluster();
      PixelContainerFromTracccClusters->push_back(xaod_pcl);
      size_t index = PixelContainerFromTracccClusters->size() - 1;

      xAOD::SpacePoint *xaod_sp = new xAOD::SpacePoint();
      SpacepointContainerFromTracccClusters->push_back(xaod_sp);

      m_tracccMeasToxAODClusterMap[meas.measurement_id] = index;
      nPix++;

      for(const unsigned int cell_idx : cluster.cell_indices()){
        const auto cell = traccc_cells.at(cell_idx);
        int phiIndex = cell.channel0();
        int etaIndex = cell.channel1();
        if(phiIndex < rowMin){rowMin = phiIndex;}
        if(phiIndex > rowMax){rowMax = phiIndex;}
        if(etaIndex < colMin){colMin = etaIndex;}
        if(etaIndex > colMax){colMax = etaIndex;}
        Identifier hit_id = m_pixelID->pixel_id(athena_ID, phiIndex, etaIndex);
        rdoList.push_back(hit_id);
      }

      const IdentifierHash Pixel_ModuleHash = m_pixelID->wafer_hash(athena_ID);
      const InDetDD::SiDetectorElement* pDE = m_pixelManager->getDetectorElement(Pixel_ModuleHash);
      const InDetDD::PixelModuleDesign* design (dynamic_cast<const InDetDD::PixelModuleDesign*>(&pDE->design()));

      Amg::Vector2D localPos(meas.local[0],meas.local[1]);

      Amg::Vector3D globalPos = pDE->globalPosition(localPos);

      double const etaWidth = static_cast<double>(colMax - colMin) / 2;
      double const phiWidth = static_cast<double>(rowMax - rowMin) / 2;
      double const etaW = design->widthFromColumnRange(colMin, colMax-1);
      double const phiW = design->widthFromRowRange(rowMin, rowMax-1);

      InDet::SiWidth siWidth(Amg::Vector2D(phiWidth,etaWidth),Amg::Vector2D(phiW,etaW));

      Amg::MatrixX cov(2,2);
      cov.setZero();

      cov(0,0) = siWidth.phiR()*siWidth.phiR()/12;
      cov(1,1) = siWidth.z()*siWidth.z()/12;
      bool split = false;
      float splitProb1 = 0;
      float splitProb2 = 0;

      ATH_MSG_DEBUG("Measurement covariance (calculated): " << cov(0,0) << " , " << cov(1,1));
      ATH_MSG_DEBUG("Measurement variance (detray): " << meas.variance[0] << " , " << meas.variance[1]);

      Eigen::Matrix<float,2,1> localPosition(localPos.x(), localPos.y());
      Eigen::Matrix<float,3,1> globalPosition(globalPos.x(), globalPos.y(), globalPos.z());
      Eigen::Matrix<float,2,2> localCovariance;
      localCovariance.setZero();
      localCovariance(0, 0) = cov(0, 0);
      localCovariance(1, 1) = cov(1, 1);

      xaod_pcl->setMeasurement<2>(Pixel_ModuleHash, localPosition, localCovariance);
      xaod_pcl->setIdentifier( rdoList.front().get_compact() );
      xaod_pcl->setRDOlist(rdoList);
      xaod_pcl->globalPosition() = globalPosition;
      xaod_pcl->setChannelsInPhiEta(siWidth.colRow()[0], siWidth.colRow()[1]);
      xaod_pcl->setWidthInEta(static_cast<float>(siWidth.widthPhiRZ()[1]));
      xaod_pcl->setIsSplit(split);
      xaod_pcl->setSplitProbabilities(splitProb1, splitProb2);

      Eigen::Matrix<double,2,1> globalVariance(0., 0.);
      xaod_sp->setSpacePoint(Pixel_ModuleHash, globalPosition, globalVariance(0,0), globalVariance(1,0), {xaod_pcl});

    }
    else {
      xAOD::StripCluster *xaod_scl = new xAOD::StripCluster();
      StripContainerFromTracccClusters->push_back(xaod_scl);
      size_t index = StripContainerFromTracccClusters->size() - 1;

      m_tracccMeasToxAODClusterMap[meas.measurement_id] = index;
      nStrip++;

      if(printClusterInfo){ATH_MSG_DEBUG("Strip cluster no. " << nStrip);}

      std::string acts_geom_id = m_AthenaToActsMap[athena_ID];
      std::istringstream stream(acts_geom_id);

      std::vector<std::string> parts;
      std::string part;
      while (std::getline(stream, part, '|')) {
          parts.push_back(part);
      }

      int volume = std::stoi(parts[0].substr(4));    // Remove "vol="
      int sensor = std::stoi(parts[2].substr(4));    // Remove "sen="

      // Make the correct athena_id
      if (volume == 22 || volume == 23 || volume == 24) {
        if (sensor % 2 == 0) { // even
          athena_ID = Identifier(athena_ID.get_compact() | 0x0008000000000);
        }
      }

      for(const unsigned int cell_idx : cluster.cell_indices()){
        const auto cell = traccc_cells.at(cell_idx);
        int phiIndex = cell.channel0();

        if(phiIndex < colMin){colMin = phiIndex;}
        if(phiIndex > colMax){colMax = phiIndex;}
        if(printClusterInfo){ATH_MSG_DEBUG("cell at this cluster: " << phiIndex);}
        Identifier hit_id = m_stripID->strip_id(athena_ID, int(phiIndex));
        rdoList.push_back(hit_id);
      }

      const IdentifierHash Strip_ModuleHash = m_stripID->wafer_hash(athena_ID);
      const InDetDD::SiDetectorElement* pDE = m_stripManager->getDetectorElement(Strip_ModuleHash);

      const InDetDD::SCT_ModuleSideDesign* design;
      if (pDE->isBarrel()){
        design = (static_cast<const InDetDD::SCT_ModuleSideDesign*>(&pDE->design()));
      } else{
        design = (static_cast<const InDetDD::StripStereoAnnulusDesign*>(&pDE->design()));
      }

      const int firstStrip = m_stripID->strip(rdoList.front());
      const int lastStrip = m_stripID->strip(rdoList.back());

      const int row = m_stripID->row(rdoList.front());
      const int firstStrip1D = design->strip1Dim (firstStrip, row );
      const int lastStrip1D = design->strip1Dim( lastStrip, row );
      const InDetDD::SiCellId cell1(firstStrip1D);
      const InDetDD::SiCellId cell2(lastStrip1D);
      const InDetDD::SiLocalPosition firstStripPos( pDE->rawLocalPositionOfCell(cell1 ));
      const InDetDD::SiLocalPosition lastStripPos( pDE->rawLocalPositionOfCell(cell2) );
      const InDetDD::SiLocalPosition centre( (firstStripPos+lastStripPos) * 0.5 );
      const double width = design->stripPitch() * ( lastStrip - firstStrip + 1 );

      const std::pair<InDetDD::SiLocalPosition, InDetDD::SiLocalPosition> ends( design->endsOfStrip(centre) );
      const double stripLength( std::abs(ends.first.xEta() - ends.second.xEta()) );

      double const phiWidth = static_cast<double>(rowMax - rowMin) / 2;
      double shift =  m_lorentzAngleTool->getLorentzShift(Strip_ModuleHash, eventContext);

      InDet::SiWidth siWidth(Amg::Vector2D(phiWidth,1), Amg::Vector2D(width,stripLength) );
      Amg::Vector2D localPos(centre.xPhi() + shift,  centre.xEta());
      if(printClusterInfo){ATH_MSG_DEBUG("local pos of this cluster: " << localPos.x() << "," << localPos.y());}

      Amg::Vector3D globalPos = pDE->globalPosition(localPos);
      Eigen::Matrix<float,1,1> localPosition;
      Eigen::Matrix<float,1,1> localCovariance;
      localCovariance.setZero();

      if (pDE->isBarrel()) {
        localPosition(0, 0) = localPos.x();
        localCovariance(0, 0) = pDE->phiPitch() * pDE->phiPitch() * (1./12.);
      } else {
        //InDetDD::SiCellId cellId = pDE->cellIdOfPosition(localPos);
        const InDetDD::StripStereoAnnulusDesign *annulusDesign = dynamic_cast<const InDetDD::StripStereoAnnulusDesign *>(&pDE->design());
        if ( annulusDesign == nullptr ){
          ATH_MSG_ERROR("Could not get strip stereo design!!");
        }

        // InDetDD::SiLocalPosition localInPolar = annulusDesign->localPositionOfCellPC(cellId);
        localPosition(0, 0) = localPos.x(); // localInPolar.xPhi();
        if(printClusterInfo){ATH_MSG_DEBUG("local pos of this cluster: " << localPosition);}
        localCovariance(0, 0) = annulusDesign->phiPitchPhi() * annulusDesign->phiPitchPhi() * (1./12.);
      }

      Eigen::Matrix<float,3,1> globalPosition(globalPos[0],globalPos[1],globalPos[2]);

      xaod_scl->setMeasurement<1>(Strip_ModuleHash, localPosition, localCovariance);
      xaod_scl->setIdentifier( rdoList.front().get_compact() );
      xaod_scl->setRDOlist(rdoList);
      xaod_scl->globalPosition() = globalPosition;
      xaod_scl->setChannelsInPhi(siWidth.colRow()[0]);

      traccc_strip_file << athena_ID << ","
        << m_stripID->barrel_ec(athena_ID) << ","
        << localPos.x() << ","
        << localPos.y() << ","
        << globalPos[0] << ","
        << globalPos[1] << ","
        << globalPos[2] << "\n";
    }

  }
  traccc_strip_file.close();

  ATH_MSG_INFO("xAOD::PixelClusterContainer with size: " << PixelContainerFromTracccClusters->size());
  ATH_MSG_INFO("xAOD::StripClusterContainer with size: " << StripContainerFromTracccClusters->size());
  return m_tracccMeasToxAODClusterMap;
}

StatusCode TrackingHitInputTool::convertToxAODClusters(const EventContext& eventContext){

  SG::WriteHandle<xAOD::PixelClusterContainer> xAODPixelContainerFromInDetClusters(m_xAODPixelClusterFromInDetClusterKey, eventContext);
  ATH_CHECK( xAODPixelContainerFromInDetClusters.record(std::make_unique<xAOD::PixelClusterContainer>(),
						 std::make_unique<xAOD::PixelClusterAuxContainer>()) );

  SG::WriteHandle<xAOD::StripClusterContainer> xAODStripContainerFromInDetClusters(m_xAODStripClusterFromInDetClusterKey, eventContext);
  ATH_CHECK( xAODStripContainerFromInDetClusters.record(std::make_unique<xAOD::StripClusterContainer>(),
						 std::make_unique<xAOD::StripClusterAuxContainer>()) );

  SG::ReadHandle<InDet::PixelClusterContainer> inputPixelClusterContainer(m_inputPixelClusterContainerKey, eventContext);
  SG::ReadHandle<InDet::SCT_ClusterContainer> inputStripClusterContainer(m_inputStripClusterContainerKey, eventContext);

  SG::ReadCondHandle<InDetDD::SiDetectorElementCollection> pixelDetEleHandle(m_pixelDetEleCollKey, eventContext);
  const InDetDD::SiDetectorElementCollection* pixElements(*pixelDetEleHandle);

  if (not pixelDetEleHandle.isValid() or pixElements==nullptr) {
    ATH_MSG_FATAL(m_pixelDetEleCollKey.fullKey() << " is not available.");
  }

  SG::ReadCondHandle<InDetDD::SiDetectorElementCollection> stripDetEleHandle(m_stripDetEleCollKey, eventContext);
  const InDetDD::SiDetectorElementCollection* stripElements(*stripDetEleHandle);
  if (not stripDetEleHandle.isValid() or stripElements==nullptr) {
    ATH_MSG_FATAL(m_stripDetEleCollKey.fullKey() << " is not available.");
  }

  for (const auto *const clusterCollection : *inputPixelClusterContainer) {
    if (!clusterCollection) continue;
    for(const auto *const theCluster : *clusterCollection)  {

      const InDetDD::SiDetectorElement* element=theCluster->detectorElement();

      xAOD::PixelCluster * pixelCl = new xAOD::PixelCluster();
      xAODPixelContainerFromInDetClusters->push_back(pixelCl);
      ATH_CHECK(convertInDetToXaodCluster(*theCluster, *element, *pixelCl));

    }
  }

  for (const auto *const clusterCollection : *inputStripClusterContainer) {
    if (!clusterCollection) continue;
    for(const auto *const theCluster : *clusterCollection)  {

      const InDetDD::SiDetectorElement *element=theCluster->detectorElement();

      xAOD::StripCluster * stripCl = new xAOD::StripCluster();
      xAODStripContainerFromInDetClusters->push_back(stripCl);
      ATH_CHECK(convertInDetToXaodCluster(*theCluster, *element, *stripCl));

    }

  }

  ATH_MSG_DEBUG("xAOD::PixelClusterContainer with size: " << xAODPixelContainerFromInDetClusters->size());
  ATH_MSG_DEBUG("xAOD::StripClusterContainer with size: " << xAODStripContainerFromInDetClusters->size());
  return StatusCode::SUCCESS;
}

StatusCode TrackingHitInputTool::convertInDetToXaodCluster(const InDet::PixelCluster& indetCluster,
				       const InDetDD::SiDetectorElement& element,
				       xAOD::PixelCluster& xaodCluster){
  IdentifierHash idHash = element.identifyHash();

  auto localPos = indetCluster.localPosition();
  auto localCov = indetCluster.localCovariance();

  Eigen::Matrix<float,2,1> localPosition(localPos.x(), localPos.y());

  Eigen::Matrix<float,2,2> localCovariance;
  localCovariance.setZero();
  localCovariance(0, 0) = localCov(0, 0);
  localCovariance(1, 1) = localCov(1, 1);

  auto globalPos = indetCluster.globalPosition();
  Eigen::Matrix<float, 3, 1> globalPosition(globalPos.x(), globalPos.y(), globalPos.z());

  const auto& RDOs = indetCluster.rdoList();
  const auto& ToTs = indetCluster.totList();
  const auto& charges = indetCluster.chargeList();
  const auto& width = indetCluster.width();
  auto isSplit = indetCluster.isSplit();
  auto splitProbability1 = indetCluster.splitProbability1();
  auto splitProbability2 = indetCluster.splitProbability2();

  xaodCluster.setMeasurement<2>(idHash, localPosition, localCovariance);
  xaodCluster.setIdentifier( indetCluster.identify().get_compact() );
  xaodCluster.setRDOlist(RDOs);
  xaodCluster.globalPosition() = globalPosition;
  xaodCluster.setToTlist(ToTs);
  //xaodCluster.setTotalToT( xAOD::xAODInDetMeasurement::Utilities::computeTotalToT(ToTs) );
  xaodCluster.setChargelist(charges);
  //xaodCluster.setTotalCharge( xAOD::xAODInDetMeasurement::Utilities::computeTotalCharge(charges) );
  xaodCluster.setLVL1A(indetCluster.LVL1A());
  xaodCluster.setChannelsInPhiEta(width.colRow()[0], width.colRow()[1]);
  xaodCluster.setWidthInEta(static_cast<float>(width.widthPhiRZ()[1]));
  xaodCluster.setIsSplit(isSplit);
  xaodCluster.setSplitProbabilities(splitProbability1, splitProbability2);

  return StatusCode::SUCCESS;
}

StatusCode TrackingHitInputTool::convertInDetToXaodCluster(const InDet::SCT_Cluster& indetCluster,
				       const InDetDD::SiDetectorElement& element,
				       xAOD::StripCluster& xaodCluster){
  static const double one_over_twelve = 1. / 12.;
  IdentifierHash idHash = element.identifyHash();

  auto localPos = indetCluster.localPosition();

  Eigen::Matrix<float,1,1> localPosition;
  Eigen::Matrix<float,1,1> localCovariance;
  localCovariance.setZero();

  if (element.isBarrel()) {
    localPosition(0, 0) = localPos.x();
    localCovariance(0, 0) = element.phiPitch() * element.phiPitch() * one_over_twelve;
  } else {
    InDetDD::SiCellId cellId = element.cellIdOfPosition(localPos);
    const InDetDD::StripStereoAnnulusDesign *design = dynamic_cast<const InDetDD::StripStereoAnnulusDesign *>(&element.design());
    if ( design == nullptr ) {
      return StatusCode::FAILURE;
    }
    InDetDD::SiLocalPosition localInPolar = design->localPositionOfCellPC(cellId);
    localPosition(0, 0) = localInPolar.xPhi();
    localCovariance(0, 0) = design->phiPitchPhi() * design->phiPitchPhi() * one_over_twelve;
  }

  auto globalPos = indetCluster.globalPosition();
  Eigen::Matrix<float, 3, 1> globalPosition(globalPos.x(), globalPos.y(), globalPos.z());

  const auto& RDOs = indetCluster.rdoList();
  const auto& width = indetCluster.width();

  xaodCluster.setMeasurement<1>(idHash, localPosition, localCovariance);
  xaodCluster.setIdentifier( indetCluster.identify().get_compact() );
  xaodCluster.setRDOlist(RDOs);
  xaodCluster.globalPosition() = globalPosition;
  xaodCluster.setChannelsInPhi(width.colRow()[0]);

  return StatusCode::SUCCESS;
}

StatusCode TrackingHitInputTool::finalize() {
  return StatusCode::SUCCESS;
}

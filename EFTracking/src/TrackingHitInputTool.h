/*
  Copyright (C) 2002-2025 CERN for the benefit of the ATLAS collaboration
*/

#ifndef EFTRACKING_TRACKINGHITINPUTTOOL_H
#define EFTRACKING_TRACKINGHITINPUTTOOL_H

// Athena includes
#include "AthenaBaseComps/AthAlgTool.h"
#include "StoreGate/ReadHandleKey.h"

#include "InDetIdentifier/PixelID.h"
#include "InDetPrepRawData/PixelClusterContainer.h"
#include "InDetIdentifier/SCT_ID.h"
#include "InDetPrepRawData/SCT_ClusterContainer.h"
#include "xAODInDetMeasurement/SpacePointContainer.h"
#include "xAODInDetMeasurement/SpacePointAuxContainer.h"

#include "InDetIdentifier/PixelID.h"
#include "InDetIdentifier/SCT_ID.h"
#include "xAODInDetMeasurement/PixelClusterContainer.h"
#include "xAODInDetMeasurement/StripClusterContainer.h"
#include "PixelReadoutGeometry/PixelDetectorManager.h"
#include "SCT_ReadoutGeometry/SCT_DetectorManager.h"
#include "InDetCondTools/ISiLorentzAngleTool.h"

#include "InDetRawData/PixelRDO_Container.h"
#include "InDetRawData/SCT_RDO_Container.h"
#include "InDetReadoutGeometry/SiDetectorDesign.h"
#include "InDetReadoutGeometry/SiDetectorManager.h"

#include "EFTracking/ITrackingHitInputTool.h"
#include <unordered_map>

class PixelID;
class SCT_ID;

class TrackingHitInputTool : public extends<AthAlgTool, ITrackingHitInputTool>
{
  public:
        TrackingHitInputTool(const std::string&, const std::string&, const IInterface*);
        virtual ~TrackingHitInputTool() { ; }

        StatusCode initialize() override;
        StatusCode finalize()   override;
        std::pair<std::vector<clusterInfo>,std::map<int,int>> readClusters( const EventContext& eventContext)  override;
        std::map<int,int> convertClusters(const EventContext& eventContext,traccc::edm::silicon_cluster_collection::host& traccc_clusters, traccc::measurement_collection_types::host& traccc_measurements, traccc::edm::silicon_cell_collection::host& traccc_cells) override;

        void setAthenaDetrayConversionMaps(
          std::unordered_map<uint64_t, Identifier> const & detray_to_athena_map,
          std::unordered_map<Identifier, uint64_t> const & athena_to_detray_map
          ) override;
  private:

        StatusCode convertInDetToXaodCluster(const InDet::PixelCluster& indetCluster,
				       const InDetDD::SiDetectorElement& element,
				       xAOD::PixelCluster& xaodCluster);
        StatusCode convertInDetToXaodCluster(const InDet::SCT_Cluster& indetCluster,
				       const InDetDD::SiDetectorElement& element,
				       xAOD::StripCluster& xaodCluster);

        StatusCode convertToxAODClusters(const EventContext& eventContext);

        SG::ReadHandleKey<InDet::PixelClusterContainer> m_inputPixelClusterContainerKey {this, "InputPixelClustersName", "ITkPixelClusters", "name of the input InDet pixel cluster container"};
        SG::ReadHandleKey<InDet::SCT_ClusterContainer>  m_inputStripClusterContainerKey {this, "InputStripClustersName", "ITkStripClusters", "name of the input InDet strip cluster container"};

        SG::WriteHandleKey<xAOD::PixelClusterContainer> m_xAODPixelClusterFromTracccClusterKey{this, "xAODPixelClusterFromTracccClusterKey","xAODPixelClustersFromTracccCluster","Traccc cluster->xAOD PixelClusters Container"};
        SG::WriteHandleKey<xAOD::StripClusterContainer> m_xAODStripClusterFromTracccClusterKey{this, "xAODStripClusterFromTracccClusterKey","xAODStripClustersFromTracccCluster","Traccc cluster ->xAOD StripClusters Container"};
        SG::WriteHandleKey<xAOD::SpacePointContainer> m_xAODSpacepointFromTracccClusterKey{this, "xAODSpacepointFromTracccClusterKey","xAODSpacepointFromTracccCluster","Traccc cluster->xAOD Spacepoint Container"};

        SG::WriteHandleKey<xAOD::PixelClusterContainer> m_xAODPixelClusterFromInDetClusterKey{this, "xAODPixelClusterFromInDetClusterKey","xAODPixelClustersFromInDetCluster","InDet cluster->xAOD PixelClusters Container"};
        SG::WriteHandleKey<xAOD::StripClusterContainer> m_xAODStripClusterFromInDetClusterKey{this, "xAODStripClusterFromInDetClusterKey","xAODStripClustersFromInDetCluster","InDet cluster ->xAOD StripClusters Container"};
        SG::WriteHandleKey<xAOD::SpacePointContainer> m_xAODSpacepointFromInDetClusterKey{this, "xAODSpacepointFromInDetClusterKey","xAODSpacepointFromInDetCluster","InDet cluster->xAOD Spacepoint Container"};

        Gaudi::Property<std::string> m_filesDir{this, "filesDir", "/eos/project/a/atlas-eftracking/GPU/ITk_data/", "full path to location with detector files"};

        ToolHandle<ISiLorentzAngleTool> m_lorentzAngleTool {this, "LorentzAngleTool", "SiLorentzAngleTool/SCTLorentzAngleTool", "Tool to retrieve Lorentz angle of SCT"};
        // Gaudi::Property<bool> m_doShift{this, "applyLorentsShift", false, "Apply Lorentz shift to cluster local position"};

        SG::ReadCondHandleKey<InDetDD::SiDetectorElementCollection> m_pixelDetEleCollKey{this, "PixelDetEleCollKey", "ITkPixelDetectorElementCollection", "Key of SiDetectorElementCollection for Pixel"};
        SG::ReadCondHandleKey<InDetDD::SiDetectorElementCollection> m_stripDetEleCollKey{this, "StripDetEleCollKey", "ITkStripDetectorElementCollection", "Key of SiDetectorElementCollection for Strip"};

        SG::ReadHandleKey<PixelRDO_Container> m_pixelRDOKey  { this, "PixelRDO", "ITkPixelRDOs" };
        SG::ReadHandleKey<SCT_RDO_Container> m_stripRDOKey  { this, "StripRDO", "ITkStripRDOs" };

        const PixelID* m_pixelID;
        const SCT_ID*  m_stripID;
        const InDetDD::PixelDetectorManager* m_pixelManager{nullptr};
        const InDetDD::SCT_DetectorManager* m_stripManager{nullptr};

        std::unordered_map<Identifier, std::uint64_t> const * m_AthenaToDetrayMap;
        std::unordered_map<std::uint64_t, Identifier> const * m_DetrayToAthenaMap;
        std::map<Identifier, std::string> m_AthenaToActsMap;
        std::map<int, int> m_tracccMeasToxAODClusterMap;
};

#endif

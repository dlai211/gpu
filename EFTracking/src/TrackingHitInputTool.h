/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#ifndef EFTRACKING_TRACKINGHITINPUTTOOL_H
#define EFTRACKING_TRACKINGHITINPUTTOOL_H

// Athena includes
#include "AthenaBaseComps/AthAlgTool.h"
#include "AthenaBaseComps/AthService.h"
#include "GaudiKernel/ServiceHandle.h"
#include "StoreGate/ReadHandleKey.h"

#include "InDetIdentifier/PixelID.h"
#include "InDetPrepRawData/PixelClusterCollection.h"
#include "InDetPrepRawData/PixelClusterContainer.h"

#include "InDetIdentifier/SCT_ID.h"
#include "InDetPrepRawData/SCT_ClusterContainer.h"


#include "InDetRawData/PixelRDO_Container.h"
#include "InDetReadoutGeometry/SiDetectorDesign.h"
#include "InDetRawData/InDetRawDataCollection.h"
#include "InDetRawData/InDetRawDataContainer.h"
#include "InDetRawData/InDetRawDataCLASS_DEF.h"
#include "InDetReadoutGeometry/SiDetectorManager.h"
#include "ReadoutGeometryBase/SiCellId.h"
#include "ReadoutGeometryBase/SiReadoutCellId.h"

#include "EFTracking/ITrackingHitInputTool.h"

class PixelID;
class SCT_ID;

class TrackingHitInputTool : public extends<AthAlgTool, ITrackingHitInputTool> {

  public:
        TrackingHitInputTool(const std::string&, const std::string&, const IInterface*);
        virtual ~TrackingHitInputTool() { ; }

        virtual StatusCode initialize() override;
        virtual StatusCode finalize()   override;
        virtual std::vector<hitInfo> readHits( const EventContext& eventContext)  override;
        virtual std::vector<clusterInfo> readClusters( const EventContext& eventContext)  override;

  private:

         SG::ReadHandleKey<InDet::PixelClusterContainer> m_inputPixelClusterContainerKey {this, "InputPixelClustersName", "ITkPixelClusters", "name of the input InDet pixel cluster container"};
         SG::ReadHandleKey<InDet::SCT_ClusterContainer>  m_inputStripClusterContainerKey {this, "InputStripClustersName", "ITkStripClusters", "name of the input InDet strip cluster container"};

         SG::ReadHandleKey<PixelRDO_Container> m_pixelRDOKey  { this, "PixelRDO", "ITkPixelRDOs" };
         SG::ReadHandleKey<SCT_RDO_Container> m_stripRDOKey  { this, "StripRDO", "ITkStripRDOs" };
         const InDetDD::SiDetectorManager* m_PIX_mgr{nullptr};
         const InDetDD::SiDetectorManager* m_SCT_mgr = nullptr;
         const PixelID* m_pixelId;
         const SCT_ID*  m_stripId;

         std::map<Identifier, std::uint64_t> m_AtlasToDetrayMap;
         std::map<std::uint64_t, Identifier> m_DetrayToAtlasMap;
};

#endif

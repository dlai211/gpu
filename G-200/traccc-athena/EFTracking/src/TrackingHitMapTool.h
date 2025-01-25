/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#ifndef EFTRACKING_TRACKINGHITMAPTOOL_H
#define EFTRACKING_TRACKINGHITMAPTOOL_H

// Athena includes
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

#include "InDetReadoutGeometry/SiDetectorDesign.h"
#include "InDetIdentifier/PixelID.h"

#include "InDetRawData/InDetRawDataCollection.h"
#include "InDetRawData/InDetRawDataContainer.h"
#include "InDetRawData/InDetRawDataCLASS_DEF.h"
#include "InDetReadoutGeometry/SiDetectorManager.h"
#include "ReadoutGeometryBase/SiCellId.h"
#include "ReadoutGeometryBase/SiReadoutCellId.h"

// Detray includes
#include "detray/io/frontend/detector_reader_config.hpp"
#include "detray/core/detector.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/geometry/tracking_surface.hpp"
#include "detray/navigation/volume_graph.hpp"

// vecmem includes
#include <vecmem/memory/host_memory_resource.hpp>

// Traccc includes
#include "traccc/geometry/silicon_detector_description.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/write.hpp"

#include "EFTracking/ITrackingHitMapTool.h"
#include <vector>

class PixelID;
class SCT_ID;


using Amg::Vector3D;
using detector_t = detray::detector<>;
using point3_t = typename detector_t::point3_type;

class TrackingHitMapTool : public extends<AthAlgTool, ITrackingHitMapTool> {

  public:
        TrackingHitMapTool(const std::string&, const std::string&, const IInterface*);
        virtual ~TrackingHitMapTool() { ; }

        virtual StatusCode initialize() override;
        virtual StatusCode finalize()   override;
        virtual std::vector<hitInfo> mapHits(std::vector<hitInfo>& hits, std::map<int, std::uint64_t>& AtlasToDetrayMap)  override;
        virtual std::vector<clusterInfo> mapClusters(std::vector<clusterInfo>& clusters, std::map<int, std::uint64_t>& AtlasToDetrayMap)  override;
        std::tuple<std::map<std::uint64_t, int>,std::map<int, std::uint64_t>,traccc::silicon_detector_description::host,std::map<int, Identifier>> createMaps();
         
  private:

         Gaudi::Property<bool> m_writeDigi{this,"writeDigi",false,"If you also want to write out the digitization cfg"};

         vecmem::host_memory_resource host_mr;
         using host_detector_type = detray::detector<>;
         host_detector_type m_detector{host_mr};

         void createAthenaMap();
         std::map<int, Amg::Vector3D >  m_atlasModuleMap;
         std::map<int, moduleInfo> m_atlasModuleInfo;
         std::map<int, Identifier> m_atlasHumanIDToIdentifier;

         void createAthenaToDetrayMap();
         std::map<std::uint64_t, int> m_DetrayToAtlasMap;
         std::map<int, std::uint64_t> m_AtlasToDetrayMap;
         std::map<std::uint64_t, Acts::GeometryIdentifier> m_DetrayToActsMap;

         std::map<std::uint64_t, Amg::Vector3D >  m_detrayModuleMap;
         traccc::silicon_detector_description::host m_dd{host_mr};
         void fill_digi_info();
         void dump_digi_cfg(traccc::silicon_detector_description::host& dd);
         void createDetrayMap();
         

         int m_AthenaHashedID;
         int m_HumanReadableID;
         std::map<Identifier, HepGeom::Point3D<double> >  m_ModuleList;

         SG::ReadCondHandleKey<InDetDD::SiDetectorElementCollection> m_pixelDetEleCollKey{this, "PixelDetEleCollKey", "ITkPixelDetectorElementCollection", "Key of SiDetectorElementCollection for Pixel"};
         SG::ReadCondHandleKey<InDetDD::SiDetectorElementCollection> m_stripDetEleCollKey{this, "SCTDetEleCollKey", "ITkStripDetectorElementCollection", "Key of SiDetectorElementCollection for SCT"};
         
         const PixelID* m_pixelId;
         const SCT_ID*  m_stripId;   
         
};

#endif

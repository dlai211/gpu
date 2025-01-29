
/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#include "StoreGate/DataHandle.h"

#include "InDetReadoutGeometry/SiDetectorDesign.h"
#include "PixelReadoutGeometry/PixelModuleDesign.h"
#include "SCT_ReadoutGeometry/SCT_BarrelModuleSideDesign.h"
#include "SCT_ReadoutGeometry/SCT_ForwardModuleSideDesign.h"

#include "InDetReadoutGeometry/SiDetectorElement.h"
#include "ReadoutGeometryBase/SiCellId.h"
#include "ReadoutGeometryBase/SiReadoutCellId.h"
#include "InDetRawData/InDetRawDataCollection.h"
#include "InDetRawData/InDetRawDataContainer.h"
#include "InDetRawData/InDetRawDataCLASS_DEF.h"
#include "InDetIdentifier/SCT_ID.h"

#include "TrkSurfaces/AnnulusBounds.h"
#include "TrkSurfaces/RectangleBounds.h"
#include "TrkSurfaces/Surface.h"
#include "TrkSurfaces/SurfaceBounds.h"
#include "TrkSurfaces/TrapezoidBounds.h"

#include "TrackingHitMapTool.h"

#include <iostream>
#include <sstream>
#include <fstream>

using point3 = detray::dpoint3D<traccc::default_algebra>;

TrackingHitMapTool::TrackingHitMapTool(const std::string& algname, const std::string& name, const IInterface* ifc) :
  base_class(algname, name, ifc) {}

StatusCode TrackingHitMapTool::initialize() {

  ATH_MSG_DEBUG("Initializing hit mapping  tool");

  ATH_CHECK(detStore()->retrieve(m_pixelId, "PixelID"));
  ATH_CHECK(detStore()->retrieve(m_stripId, "SCT_ID"));
  ATH_CHECK(m_pixelDetEleCollKey.initialize());
  ATH_CHECK(m_stripDetEleCollKey.initialize());

  // Read the detector.
  detray::io::detector_reader_config reader_cfg{};
  reader_cfg.add_file("/eos/project/a/atlas-eftracking/GPU/ITk_data/ITk_DetectorBuilder_geometry.json");
  auto [host_det, _] = detray::io::read_detector<host_detector_type>(host_mr, reader_cfg);
  m_detector = std::move(host_det);

  ATH_MSG_DEBUG("Initializing complete");

  return StatusCode::SUCCESS;

}

std::tuple<std::map<std::uint64_t, int>,std::map<int, std::uint64_t>,traccc::silicon_detector_description::host,std::map<int, Identifier>> TrackingHitMapTool::createMaps(){

  createAthenaMap();
  createDetrayMap();
  createAthenaToDetrayMap();
  fill_digi_info();
  
  return std::make_tuple(m_DetrayToAtlasMap,m_AtlasToDetrayMap,m_dd,m_atlasHumanIDToIdentifier);

}

std::vector<clusterInfo> TrackingHitMapTool::mapClusters(std::vector<clusterInfo>& clusters, std::map<int, std::uint64_t>& AtlasToDetrayMap){

  typename detector_t::geometry_context context{};

  int not_matched = 0;
  int not_filled = 0;
  ATH_MSG_INFO("Number of clusters to match: " << clusters.size());
  for(std::vector<clusterInfo>::size_type i = 0; i < clusters.size(); i++){

    clusterInfo &thiscluster = clusters[i];
    auto atlas_id = thiscluster.atlas_id;
    if(!atlas_id){
      not_filled++;
      ATH_MSG_DEBUG("Could not match this cluster: " << i);
    }
    auto detray_id = AtlasToDetrayMap[atlas_id];
    if (!detray_id){
      ATH_MSG_DEBUG("Could not match this id: " << atlas_id);
      not_matched++;
      continue;
    }
    thiscluster.detray_id = detray_id;
    const point3 gl_pos{thiscluster.globalPosition[0],thiscluster.globalPosition[1],thiscluster.globalPosition[2]};
    const auto& sf = detray::geometry::barcode{detray_id};
    const detray::tracking_surface<detray::detector<> > surface{m_detector, sf};
    const auto loc_pos = surface.transform(context).point_to_local(gl_pos);
    thiscluster.localPosition[0] = loc_pos[0];
    thiscluster.localPosition[1] = loc_pos[1];

  }

  ATH_MSG_INFO("Number of unmatched clusters: " << not_matched);
  ATH_MSG_INFO("Number of unfilled clusters : " << not_filled);
  return clusters;
}

std::vector<hitInfo> TrackingHitMapTool::mapHits(std::vector<hitInfo>& hits, std::map<int, std::uint64_t>& AtlasToDetrayMap){

  int not_matched = 0;
  int not_filled = 0;

  ATH_MSG_INFO("Number of hits to match: " << hits.size());
  for(std::vector<hitInfo>::size_type i = 0; i < hits.size(); i++){ 

    hitInfo &thishit = hits[i];
    auto atlas_id = thishit.atlas_id;
    if(!atlas_id){
      not_filled++;
      ATH_MSG_DEBUG("Could not match this id: " << atlas_id);
    }
    if (AtlasToDetrayMap.find(atlas_id) != AtlasToDetrayMap.end()) {

      auto detray_id = AtlasToDetrayMap.at(atlas_id);
      thishit.detray_id = detray_id;

    } else {
      ATH_MSG_INFO("Could not match this id: " << atlas_id);
      not_matched++;
      thishit.detray_id = 0;
    }

  }

  ATH_MSG_INFO("Number of unmatched hits: " << not_matched);
  ATH_MSG_INFO("Number of unfilled hits : " << not_filled);
  return hits;
}

StatusCode TrackingHitMapTool::finalize() {
  return StatusCode::SUCCESS;
}

void TrackingHitMapTool::createAthenaToDetrayMap(){

  std::ifstream map_file("/eos/project/a/atlas-eftracking/GPU/ITk_data/athenaToDetrayMap.txt");

  ATH_MSG_INFO("Constructing athena to detray map");
  std::string str; 
  while (std::getline(map_file, str)){
    
    ATH_MSG_DEBUG(str);
    std::stringstream test(str);
    std::vector<std::string> seglist;
    std::string segment;
    while(std::getline(test, segment, ',')){
      seglist.push_back(segment);
    }

    std::string a_id = seglist[0];
    std::string d_id = seglist[1];
    ATH_MSG_DEBUG( a_id << " mapps to " << d_id);
    uint64_t idD = std::stoull(d_id);
    int idA = std::stoi(a_id);
    m_AtlasToDetrayMap[idA] = idD;
    m_DetrayToAtlasMap[idD] = idA;

  }

}

void TrackingHitMapTool::createAthenaMap() {

  int nStrip = 0;
  int nPixel = 0;

  SG::ReadCondHandle<InDetDD::SiDetectorElementCollection> stripDetEleHandle(m_stripDetEleCollKey);
  const InDetDD::SiDetectorElementCollection* strip_elements(*stripDetEleHandle);
  if (not stripDetEleHandle.isValid() or strip_elements==nullptr) {
     ATH_MSG_FATAL(m_stripDetEleCollKey.fullKey() << " is not available.");
  }
  for (const InDetDD::SiDetectorElement *element: *strip_elements) {
    const Identifier Strip_ModuleID = m_stripId->module_id(element->identify()); //from wafer id to module id
    const IdentifierHash Strip_ModuleHash = m_stripId->wafer_hash(Strip_ModuleID);
    if (m_ModuleList.find(Strip_ModuleID) == m_ModuleList.end()){

      if (m_stripId->side(element->identify()) != 0) { continue; } // not inner side case

      // Write out Visualization Lookup Tree
      m_AthenaHashedID = Strip_ModuleID.get_identifier32().get_compact();
      m_HumanReadableID = 1000000*2 //2 = Strip
      + 100000*m_stripId->layer_disk(Strip_ModuleID)
      + 1000*(10+m_stripId->eta_module(Strip_ModuleID))
      + m_stripId->phi_module(Strip_ModuleID);
      if ( m_stripId->barrel_ec(Strip_ModuleID) != 0 ) {
        m_HumanReadableID = m_stripId->barrel_ec(Strip_ModuleID)*(m_HumanReadableID + 10000000);
      }

      m_atlasHumanIDToIdentifier[m_HumanReadableID] = Strip_ModuleID;
      
      if(m_atlasModuleInfo.find(m_HumanReadableID) != m_atlasModuleInfo.end()){continue;}

      ++nStrip;
      ATH_MSG_VERBOSE( "Strip module " << nStrip );

      moduleInfo thismod;
      //const InDetDD::SCT_ModuleSideDesign* s_design;
      if (m_stripId->barrel_ec(Strip_ModuleID) == 0) {
        const InDetDD::SCT_BarrelModuleSideDesign* s_design = (static_cast<const InDetDD::SCT_BarrelModuleSideDesign*>(&element->design()));
        auto boundsType = element->bounds().type();
        thismod.thickness = element->thickness();
        thismod.pixel = false;

        if (boundsType == Trk::SurfaceBounds::Rectangle) {

            ATH_MSG_DEBUG("Strip module in barrel with ID " << m_HumanReadableID << " is designed like a rectangle");
            
            thismod.module_width = s_design->width();
            thismod.module_length = s_design->length();
            thismod.rows = (s_design->length())/(s_design->stripPitch());

        }
        if (boundsType == Trk::SurfaceBounds::Trapezoid) {
            
            thismod.module_width = s_design->width();
            thismod.module_length = s_design->length();
            thismod.rows = (s_design->length())/(s_design->stripPitch());
            
            ATH_MSG_DEBUG("Strip module in barrel with ID " << m_HumanReadableID << " is designed like a trapezoid");
            ATH_MSG_DEBUG("width " << s_design->width());
            ATH_MSG_DEBUG("lenght / pitch " << s_design->length() << " / " << s_design->stripPitch());
        }
        if (boundsType == Trk::SurfaceBounds::Annulus)   {
            
            thismod.module_width = s_design->width();
            thismod.module_length = s_design->length();
            thismod.rows = (s_design->length())/(s_design->stripPitch());
            
            ATH_MSG_DEBUG("Strip module in barrel with ID " << m_HumanReadableID << " is designed like a annulus");
            ATH_MSG_DEBUG("width " << s_design->width());
            ATH_MSG_DEBUG("lenght / pitch " << s_design->length() << " / " << s_design->stripPitch());
        }
      } else {
        const InDetDD::SCT_ForwardModuleSideDesign* s_design = (static_cast<const InDetDD::SCT_ForwardModuleSideDesign*>(&element->design()));
        auto boundsType = element->bounds().type();
        thismod.thickness = element->thickness();
        thismod.pixel = false;

        if (boundsType == Trk::SurfaceBounds::Rectangle) {
            ATH_MSG_DEBUG("Strip module FWD with ID " << m_HumanReadableID << " is designed like a rectangle");
            thismod.module_width = s_design->width();
            thismod.module_length = s_design->length();
            thismod.rows = (s_design->length())/(s_design->stripPitch());

        }
        if (boundsType == Trk::SurfaceBounds::Trapezoid) {
            
            thismod.module_width = s_design->width();
            thismod.module_length = s_design->length();
            thismod.rows = (s_design->length())/(s_design->stripPitch());
            
            ATH_MSG_DEBUG("Strip module FWD with ID " << m_HumanReadableID << " is designed like a trapezoid");
            ATH_MSG_DEBUG("width " << s_design->width());
            ATH_MSG_DEBUG("lenght / pitch " << s_design->length() << " / " << s_design->stripPitch());
        }
        if (boundsType == Trk::SurfaceBounds::Annulus)   {
            
            thismod.module_width = s_design->width();
            thismod.module_length = s_design->length();
            thismod.rows = (s_design->length())/(s_design->stripPitch());
            
            ATH_MSG_DEBUG("Strip module FWD with ID " << m_HumanReadableID << " is designed like a annulus");
            ATH_MSG_DEBUG("width " << s_design->width());
            ATH_MSG_DEBUG("lenght / pitch " << s_design->length() << " / " << s_design->stripPitch());
        }
          
      }
      m_atlasModuleInfo[m_HumanReadableID] = thismod;
     
      const InDetDD::SiDetectorElement *module = strip_elements->getDetectorElement(Strip_ModuleHash);
      std::vector<const Trk::Surface*> module_surfaces = module->surfaces();

      m_atlasModuleMap[m_HumanReadableID][0] = module->center()[0];
      m_atlasModuleMap[m_HumanReadableID][1] = module->center()[1];
      m_atlasModuleMap[m_HumanReadableID][2] = module->center()[2];
      m_atlasHumanIDToIdentifier[m_HumanReadableID] = Strip_ModuleID;
      
      ATH_MSG_VERBOSE( "Human Readable ID: " << m_HumanReadableID << module->center()[0] << "," << module->center()[1] << "," << module->center()[2]);
             
    }
  }

  SG::ReadCondHandle<InDetDD::SiDetectorElementCollection> pixelDetEleHandle(m_pixelDetEleCollKey);
  const InDetDD::SiDetectorElementCollection* pixel_elements(*pixelDetEleHandle);
  if (not pixelDetEleHandle.isValid() or pixel_elements==nullptr) {
    ATH_MSG_FATAL(m_pixelDetEleCollKey.fullKey() << " is not available.");
  }
  for (const InDetDD::SiDetectorElement *element: *pixel_elements) {
    // get the ID
    const Identifier Pixel_ModuleID = element->identify();
    const IdentifierHash Pixel_ModuleHash = m_pixelId->wafer_hash(Pixel_ModuleID);
    if (Pixel_ModuleID.is_valid()) {
      if (m_ModuleList.find(Pixel_ModuleID) == m_ModuleList.end()) {

        const InDetDD::SiDetectorElement *module = pixel_elements->getDetectorElement(Pixel_ModuleHash);
        std::vector<const Trk::Surface*> module_surfaces = module->surfaces();
        
        ++nPixel;
        ATH_MSG_DEBUG( "Pixel module " << nPixel );

        // Write out Visualization Lookup Tree
        m_AthenaHashedID = Pixel_ModuleID.get_identifier32().get_compact();
        m_HumanReadableID = 1000000*1 /*1 = Pixel*/
        + 100000*m_pixelId->layer_disk(Pixel_ModuleID)
        + 1000*(10+m_pixelId->eta_module(Pixel_ModuleID))
        + m_pixelId->phi_module(Pixel_ModuleID);
        if ( m_pixelId->barrel_ec(Pixel_ModuleID) != 0 ) {
          m_HumanReadableID = m_pixelId->barrel_ec(Pixel_ModuleID)*(m_HumanReadableID + 10000000);
        }
        
        m_atlasHumanIDToIdentifier[m_HumanReadableID] = Pixel_ModuleID;
        
        if(m_atlasModuleInfo.find(m_HumanReadableID) != m_atlasModuleInfo.end()){continue;}

        moduleInfo thismod;
        const InDetDD::PixelModuleDesign *p_design = static_cast<const InDetDD::PixelModuleDesign*>(&element->design());
        auto boundsType = module->bounds().type();
        thismod.thickness = module->thickness();
        thismod.pixel = true;

        if (boundsType == Trk::SurfaceBounds::Rectangle) {
          
          thismod.module_width = p_design->width();
          thismod.module_length = p_design->length();
          thismod.rows = p_design->rows();
          thismod.columns = p_design->columns();

        }
        if (boundsType == Trk::SurfaceBounds::Trapezoid) {
          
          thismod.module_width = p_design->width();
          thismod.module_length = p_design->length();
          thismod.rows = p_design->rows();
          thismod.columns = p_design->columns();
          
          ATH_MSG_DEBUG("Pixel module with ID " << m_HumanReadableID << " is designed like a trapezoid");
          ATH_MSG_DEBUG("pitch x / columns " << p_design->width() << " / " << p_design->columns());
          ATH_MSG_DEBUG("pitch y / rows " << p_design->length() << " / " << p_design->rows());
        }
        if (boundsType == Trk::SurfaceBounds::Annulus)   {
          
          thismod.module_width = p_design->width();
          thismod.module_length = p_design->length();
          thismod.rows = p_design->rows();
          thismod.columns = p_design->columns();
          
          ATH_MSG_DEBUG("Pixel module with ID " << m_HumanReadableID << " is designed like a annulus");
          ATH_MSG_DEBUG("pitch x / columns " << p_design->width() << " / " << p_design->columns());
          ATH_MSG_DEBUG("pitch y / rows " << p_design->length() << " / " << p_design->rows());
        }
        m_atlasModuleInfo[m_HumanReadableID] = thismod;
      } 
    } else {
       ATH_MSG_INFO( "Not a valid PIXEL Module ID (setup)" );
    }
  }

  ATH_MSG_INFO(  nPixel << " Atlas Pixel modules found." );
  ATH_MSG_INFO(  nStrip   << " Atlas Strip  modules found." );
  ATH_MSG_INFO( "Wrote segmentation info for " << m_atlasModuleInfo.size() << " modules");

}

void TrackingHitMapTool::dump_digi_cfg(traccc::silicon_detector_description::host& dd){

    std::vector<std::pair<Acts::GeometryIdentifier, traccc::module_digitization_config>> module_digi;

    for(unsigned int i = 0; i < dd.acts_geometry_id().size(); i++){

        traccc::module_digitization_config config;
        Acts::BinUtility bUtility;
        bUtility += Acts::BinUtility(
        (-2*dd.reference_x().at(i))/((dd.pitch_x().at(i))), (dd.reference_x().at(i)), (-dd.reference_x().at(i)), Acts::open, Acts::BinningValue::binX // X binning
        );
        bUtility += Acts::BinUtility(
        1, -1.0, 1.0, Acts::open, Acts::BinningValue::binY
        // (-2*dd.reference_y().at(i))/((dd.pitch_y().at(i))), (dd.reference_y().at(i)), (-dd.reference_y().at(i)), Acts::open, Acts::BinningValue::binY // X binning
        );
        config.segmentation = bUtility;
        module_digi.emplace_back(dd.acts_geometry_id().at(i),config);

    }

    const traccc::digitization_config digi_conf(module_digi);
    traccc::io::write("ITk_digitization_config_with_strips.json",traccc::data_format::json,digi_conf);

}

void TrackingHitMapTool::fill_digi_info(){

    ATH_MSG_INFO("Filling detector description");
    typename detector_t::geometry_context context{};

    m_dd.reserve(m_detector.surfaces().size());
    ATH_MSG_INFO("Detector has " << m_detector.surfaces().size() << " surfaces");

    for (const auto& surface : m_detector.surfaces()) {
      auto sf = detray::tracking_surface{m_detector, surface};
      if (!surface.is_sensitive()){continue;}
      
      auto id = surface.barcode().value();

      const auto geo_id = surface.source;
      Acts::GeometryIdentifier acts_geom_id{geo_id};

      m_dd.resize(m_dd.size() + 1);
      m_dd.geometry_id().back() = detray::geometry::barcode{surface.barcode()};
      m_dd.placement().back() = sf.transform(context);
      m_dd.acts_geometry_id().back() = id;//surface.source;

      int a_mod = m_DetrayToAtlasMap[id];

      moduleInfo thismod = m_atlasModuleInfo[a_mod];

      // Set a hard-coded threshold for which cells should be considered for
      // clusterization on this module / surface.
      m_dd.threshold().back() = 0.f;
      m_dd.reference_x().back() = (-thismod.module_width*0.5);
      m_dd.reference_y().back() = (-thismod.module_length*0.5);
      // m_dd.pitch_y().back() = thismod.module_length/(thismod.columns);
      m_dd.pitch_y().back() = 1.f;
      m_dd.pitch_x().back() = thismod.module_width/(thismod.rows);

      m_dd.dimensions().back() = 2;

      // dump_digi_cfg(m_dd);
    }

    std::vector<std::pair<Acts::GeometryIdentifier, traccc::module_digitization_config>> module_digi;

    for(unsigned int i = 0; i < m_dd.acts_geometry_id().size(); i++){

        traccc::module_digitization_config config;
        Acts::BinUtility bUtility;
        bUtility += Acts::BinUtility(
        (-2*m_dd.reference_x().at(i))/((m_dd.pitch_x().at(i))), (m_dd.reference_x().at(i)), (-m_dd.reference_x().at(i)), Acts::open, Acts::BinningValue::binX // X binning
        );
        bUtility += Acts::BinUtility(
        1, -1.0, 1.0, Acts::open, Acts::BinningValue::binY
        // (-2*m_dd.reference_y().at(i))/((m_dd.pitch_y().at(i))), (m_dd.reference_y().at(i)), (-m_dd.reference_y().at(i)), Acts::open, Acts::BinningValue::binY // X binning
        );
        config.segmentation = bUtility;
        module_digi.emplace_back(m_dd.acts_geometry_id().at(i),config);

     }

    const traccc::digitization_config digi_conf(module_digi);

    traccc::io::write("/eos/user/j/jlai/g200/gpu/G-200/run/ITk_digitization_config_with_strips.json",traccc::data_format::json,digi_conf);


    //   if(thismod.pixel){m_dd.dimensions().back() = 2;}
    //   else{m_dd.dimensions().back() = 1;}

    // }
    // if(m_writeDigi){
    //   dump_digi_cfg(m_dd);
    // }
}

void TrackingHitMapTool::createDetrayMap(){

  ATH_MSG_INFO("Filling detray map");
  std::stringstream ss;

  using scalar_t = typename detector_t::scalar_type;
  typename detector_t::geometry_context context{};

  int n_surfaces = 0;
  
  for (const auto& sf_desc : m_detector.surfaces()) {
    auto sf = detray::tracking_surface{m_detector, sf_desc};
    if (!sf_desc.is_sensitive()){continue;}

    const auto geo_id = sf_desc.source;
    Acts::GeometryIdentifier acts_geom_id{geo_id};
    
    const auto position = sf.transform(context).point_to_global(sf.centroid());
    auto id = sf_desc.barcode().value();
    if(std::isnan(position[0]) || std::isnan(position[1]) || std::isnan(position[2])){
      ATH_MSG_INFO("Issue translating module center to global coordinates");
      continue;
    }
    
    m_DetrayToActsMap[id] = acts_geom_id;
    n_surfaces++;
  }
  ATH_MSG_INFO(  n_surfaces << " Detray surfaces found." );

}

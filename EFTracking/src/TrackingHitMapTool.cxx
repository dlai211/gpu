
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
#include <nlohmann/json.hpp>  // JSON library


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

std::tuple<std::map<std::uint64_t, Identifier>,std::map<Identifier, std::uint64_t>,traccc::silicon_detector_description::host,std::map<Identifier, Identifier>> TrackingHitMapTool::createMaps(){
 
  createAthenaMap();
  createDetrayMap();
  createAthenaToDetrayMap();
  fill_digi_info();
  
  return std::make_tuple(m_DetrayToAtlasMap,m_AtlasToDetrayMap,m_dd,m_atlasHumanIDToIdentifier);

}

std::vector<clusterInfo> TrackingHitMapTool::mapClusters(std::vector<clusterInfo>& clusters, std::map<Identifier, std::uint64_t>& AtlasToDetrayMap){

  typename detector_t::geometry_context context{};

  int not_matched = 0;
  int not_filled = 0;
  ATH_MSG_INFO("Number of clusters to match: " << clusters.size());
  for(std::vector<clusterInfo>::size_type i = 0; i < clusters.size(); i++){

    clusterInfo &thiscluster = clusters[i];
    auto atlas_id = thiscluster.atlas_id;
    if(atlas_id == Identifier()){
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

std::vector<hitInfo> TrackingHitMapTool::mapHits(std::vector<hitInfo>& hits, std::map<Identifier, std::uint64_t>& AtlasToDetrayMap){

  int not_matched = 0;
  int not_filled = 0;

  ATH_MSG_INFO("Number of hits to match: " << hits.size());
  for(std::vector<hitInfo>::size_type i = 0; i < hits.size(); i++){ 

    hitInfo &thishit = hits[i];
    auto atlas_id = thishit.atlas_id;
    if(atlas_id == Identifier()){
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

  std::ifstream map_file("/eos/project/a/atlas-eftracking/GPU/ITk_data/athenaIdentifierToDetrayMap.txt");

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
    Identifier idA;
    idA.set(a_id);
    m_AtlasToDetrayMap[idA] = idD;
    m_DetrayToAtlasMap[idD] = idA;

  }
}


void TrackingHitMapTool::createAthenaMap() {

  int nStrip = 0;
  int nPixel = 0;
  int nRect1 = 0;
  int nTrap1 = 0;
  int nAnnu1 = 0;
  int nnonbar = 0;
  int nbar = 0;
  int nRect2 = 0;
  int nTrap2 = 0;
  int nAnnu2 = 0;

  // Open CSV file for writing
  std::string csv_filename = "/eos/user/j/jlai/g200/gpu/G-200/traccc-athena/run/ITk_athena_map.csv";
  std::ofstream csv_file(csv_filename);
  if (!csv_file.is_open()) {
    ATH_MSG_FATAL("Failed to open CSV file for writing: " << csv_filename);
    return;
  }

  // Write CSV Header
  csv_file << "Module_ID,Barrel_ec,Layer_Disk,Phi,Eta,Side,StripWafer_Hash\n";
  

  SG::ReadCondHandle<InDetDD::SiDetectorElementCollection> stripDetEleHandle(m_stripDetEleCollKey);
  const InDetDD::SiDetectorElementCollection* strip_elements(*stripDetEleHandle);
  if (not stripDetEleHandle.isValid() or strip_elements==nullptr) {
     ATH_MSG_FATAL(m_stripDetEleCollKey.fullKey() << " is not available.");
  }
  for (const InDetDD::SiDetectorElement *element: *strip_elements) {
    const Identifier Strip_ModuleID = m_stripId->module_id(element->identify()); //from wafer id to module id
    const IdentifierHash Strip_ModuleHash = m_stripId->wafer_hash(Strip_ModuleID);

    // Extract side information (0 = inner, 1 = outer)
    int layer_disk = m_stripId->layer_disk(element->identify());
    int eta_module = m_stripId->eta_module(element->identify());
    int phi_module = m_stripId->phi_module(element->identify());
    int barrel_ec = m_stripId->barrel_ec(element->identify());
    int side = m_stripId->side(element->identify());
    int strip_tmp = m_stripId->strip(element->identify());


    if (m_ModuleList.find(Strip_ModuleID) == m_ModuleList.end()){

      // if (side != 0) { continue; }
      // if (m_stripId->side(element->identify()) != 0) { continue; } // not inner side case

      csv_file << Strip_ModuleID << ","
               << barrel_ec << ","
               << layer_disk << ","
               << phi_module << ","
               << eta_module << ","
               << side << ","
               << strip_tmp << ","
               << Strip_ModuleHash << "\n";

      // m_atlasHumanIDToIdentifier[m_HumanReadableID] = Strip_ModuleID;
      m_atlasHumanIDToIdentifier[Strip_ModuleID] = Strip_ModuleID;

      // if(m_atlasModuleInfo.find(m_HumanReadableID) != m_atlasModuleInfo.end()){
        // ATH_MSG_INFO("repeating sideAwareID, Weird");
        // continue;}

      ++nStrip;
      ATH_MSG_VERBOSE( "Strip module " << nStrip );

      moduleInfo thismod;
      // thismod.side = side; // store side information

      //const InDetDD::SCT_ModuleSideDesign* s_design;
      if (m_stripId->barrel_ec(Strip_ModuleID) == 0) {
        ++nbar;
        const InDetDD::SCT_BarrelModuleSideDesign* s_design = (static_cast<const InDetDD::SCT_BarrelModuleSideDesign*>(&element->design()));
        auto boundsType = element->bounds().type();
        thismod.thickness = element->thickness();
        thismod.pixel = false; 
        thismod.side = side; // store side information

        if (nbar <= 15) {
          ATH_MSG_INFO("barrel: " << m_stripId->barrel_ec(Strip_ModuleID) << " Module ID: " << Strip_ModuleID << 
          ", Side: " << side << ", Side: " << thismod.side << " Thickness: " << thismod.thickness);
        }  


        if (boundsType == Trk::SurfaceBounds::Rectangle) {
            ++nRect1;
            ATH_MSG_DEBUG("Strip module in barrel with ID " << Strip_ModuleID << " is designed like a rectangle");
            
            thismod.module_width = s_design->width();
            thismod.module_length = s_design->length();
            thismod.rows = s_design->cells(); // same as s_design->width() / s_design_stripPitch();

            if (nRect1 < 10) {
              ATH_MSG_INFO(" nRect1, s_design->width(): " << s_design->width() << 
                           " nRect1, s_design->length(): " << s_design->length() <<
                           " nRect1, s_design->cells(): " << s_design->cells() <<
                           " nRect1, s_design->stripPitch(): " << s_design->stripPitch() << 
                           " nRect1, s_design->stripLength(): " << s_design->stripLength() );
            }

        }
        if (boundsType == Trk::SurfaceBounds::Trapezoid) {
            ++nTrap1;
            thismod.module_width = s_design->width();
            thismod.module_length = s_design->length();
            thismod.rows = s_design->cells();
            
            ATH_MSG_DEBUG("Strip module in barrel with ID " << Strip_ModuleID << " is designed like a trapezoid");
            ATH_MSG_DEBUG("width " << s_design->width());
            ATH_MSG_DEBUG("lenght / pitch " << s_design->length() << " / " << s_design->stripPitch());
        }
        if (boundsType == Trk::SurfaceBounds::Annulus)   {
            ++nAnnu1;
            thismod.module_width = s_design->width();
            thismod.module_length = s_design->length();
            thismod.rows = s_design->cells();
            
            ATH_MSG_DEBUG("Strip module in barrel with ID " << Strip_ModuleID << " is designed like a annulus");
            ATH_MSG_DEBUG("width " << s_design->width());
            ATH_MSG_DEBUG("lenght / pitch " << s_design->length() << " / " << s_design->stripPitch());
        }
      } else {
        ++nnonbar;
        const InDetDD::SCT_ForwardModuleSideDesign* s_design = (static_cast<const InDetDD::SCT_ForwardModuleSideDesign*>(&element->design()));
        auto boundsType = element->bounds().type();
        thismod.thickness = element->thickness();
        thismod.pixel = false;
        thismod.side = side; 

        if (nnonbar <= 15) {
          ATH_MSG_INFO("barrel: " << m_stripId->barrel_ec(Strip_ModuleID) << " Module ID: " << Strip_ModuleID << 
          ", Side: " << side << ", Side: " << thismod.side << " Thickness: " << thismod.thickness);
        }  

        if (boundsType == Trk::SurfaceBounds::Rectangle) {
            ++nRect2;
            ATH_MSG_DEBUG("Strip module FWD with ID " << Strip_ModuleID << " is designed like a rectangle");
            thismod.module_width = s_design->width();
            thismod.module_length = s_design->length();
            thismod.rows = s_design->cells();

        }
        if (boundsType == Trk::SurfaceBounds::Trapezoid) {
            ++nTrap2;
            thismod.module_width = s_design->width();
            thismod.module_length = s_design->length();
            thismod.rows = s_design->cells();
            
            ATH_MSG_DEBUG("Strip module FWD with ID " << Strip_ModuleID << " is designed like a trapezoid");
            ATH_MSG_DEBUG("width " << s_design->width());
            ATH_MSG_DEBUG("lenght / pitch " << s_design->length() << " / " << s_design->stripPitch());
        }
        if (boundsType == Trk::SurfaceBounds::Annulus)   {
            ++nAnnu2;
            thismod.module_width = s_design->width();
            thismod.module_length = s_design->length();
            thismod.rows = s_design->cells();

            if (nAnnu2 < 10) {
              ATH_MSG_INFO(" nAnnu2, s_design->width(): " << s_design->width() << 
                           " nAnnu2, s_design->length(): " << s_design->length() <<
                           " nAnnu2, s_design->cells(): " << s_design->cells() <<
                           " nAnnu2, s_design->stripPitch(): " << s_design->stripPitch());
            }
        }
          
      }
      m_atlasModuleInfo[Strip_ModuleID] = thismod;
     
      const InDetDD::SiDetectorElement *module = strip_elements->getDetectorElement(Strip_ModuleHash);
      std::vector<const Trk::Surface*> module_surfaces = module->surfaces();

      m_atlasModuleMap[Strip_ModuleID][0] = module->center()[0];
      m_atlasModuleMap[Strip_ModuleID][1] = module->center()[1];
      m_atlasModuleMap[Strip_ModuleID][2] = module->center()[2];
      m_atlasHumanIDToIdentifier[Strip_ModuleID] = Strip_ModuleID;
      
      ATH_MSG_VERBOSE( "Human Readable ID: " << Strip_ModuleID << module->center()[0] << "," << module->center()[1] << "," << module->center()[2]);
             
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
        
        m_atlasHumanIDToIdentifier[Pixel_ModuleID] = Pixel_ModuleID;
        
        if(m_atlasModuleInfo.find(Pixel_ModuleID) != m_atlasModuleInfo.end()){continue;}

        moduleInfo thismod;
        const InDetDD::PixelModuleDesign *p_design = static_cast<const InDetDD::PixelModuleDesign*>(&element->design());
        auto boundsType = module->bounds().type();
        thismod.thickness = module->thickness();
        thismod.pixel = true;
        thismod.side = 2; // for pixel

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
          
          ATH_MSG_DEBUG("Pixel module with ID " << Pixel_ModuleID << " is designed like a trapezoid");
          ATH_MSG_DEBUG("pitch x / columns " << p_design->width() << " / " << p_design->columns());
          ATH_MSG_DEBUG("pitch y / rows " << p_design->length() << " / " << p_design->rows());
        }
        if (boundsType == Trk::SurfaceBounds::Annulus)   {
          
          thismod.module_width = p_design->width();
          thismod.module_length = p_design->length();
          thismod.rows = p_design->rows();
          thismod.columns = p_design->columns();
          
          ATH_MSG_DEBUG("Pixel module with ID " << Pixel_ModuleID << " is designed like a annulus");
          ATH_MSG_DEBUG("pitch x / columns " << p_design->width() << " / " << p_design->columns());
          ATH_MSG_DEBUG("pitch y / rows " << p_design->length() << " / " << p_design->rows());
        }
        m_atlasModuleInfo[Pixel_ModuleID] = thismod;
      } 
    } else {
       ATH_MSG_INFO( "Not a valid PIXEL Module ID (setup)" );
    }
  }

  ATH_MSG_INFO(  nPixel << " Atlas Pixel modules found." );
  ATH_MSG_INFO(  nStrip   << " Atlas Strip  modules found." );
  ATH_MSG_INFO(" Barrel_ec != 0 " << nnonbar << " Barrel_ec == 0 " << nbar);

  ATH_MSG_INFO("Surface Bounds Type Counts: ");
  ATH_MSG_INFO(" Rectangle1: " << nRect1 << " Rectangle2: " << nRect2);
  ATH_MSG_INFO(" Trapezoid1: " << nTrap1 << " Trapezoid2: " << nTrap2);
  ATH_MSG_INFO(" Annulus1: " << nAnnu1 << " Annulus2: " << nAnnu2);
  ATH_MSG_INFO( "Wrote segmentation info for " << m_atlasModuleInfo.size() << " modules");

  // Close CSV file
  csv_file.close();
  ATH_MSG_INFO("Successfully saved Athena map data to CSV: " << csv_filename);
  

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
        (-2*dd.reference_y().at(i))/((dd.pitch_y().at(i))), (dd.reference_y().at(i)), (-dd.reference_y().at(i)), Acts::open, Acts::BinningValue::binY // X binning
        );
        config.segmentation = bUtility;
        module_digi.emplace_back(dd.acts_geometry_id().at(i),config);

    }

    const traccc::digitization_config digi_conf(module_digi);
    traccc::io::write("ITk_digitization_config_with_strips.json",traccc::data_format::json,digi_conf);

}

void TrackingHitMapTool::fill_digi_info(){

    ATH_MSG_DEBUG("Filling detector description");
    typename detector_t::geometry_context context{};

    m_dd.reserve(m_detector.surfaces().size());
    int pcount = 0;
    int scount = 0;
    std::set<Identifier> processed_a_mods;  // set to track processed Athena Geometry IDs
    std::map<Identifier, int> a_mod_count;

    // Open the ActsToAthena map
    std::ifstream map_file("/eos/project/a/atlas-eftracking/GPU/ITk_data/actsToAthenaIdentifierMap.txt");
    if (!map_file.is_open()) {
      ATH_MSG_ERROR("Failed to open actsToAthenaIdentifierMap.txt");
      return;
    }

    // Read the ActsToAthena mapping and store it in a temporary map
    std::map<Identifier, std::tuple<int, int, int>> AthenaGeometryMap;  // Key: a_mod, Value: (volume, layer, sensor)
    std::string line;

    while (std::getline(map_file, line)) {

        std::stringstream ss(line);
        std::vector<std::string> segments;
        std::string temp;

        // Split line by '|' first, then extract the last value after ','
        while (std::getline(ss, temp, '|')) {
            segments.push_back(temp);
        }

        // Ensure correct format
        if (segments.size() < 3) {
            ATH_MSG_WARNING("Invalid format in actsToAthenaIdentifierMap.txt: " << line);
            continue;
        }

        // Extract volume, layer, sensor ID
        int volume = std::stoi(segments[0].substr(4));    // Remove "vol="
        int layer = std::stoi(segments[1].substr(4));     // Remove "lay="
        int sensor = std::stoi(segments[2].substr(4, segments[2].find(',')));  // Extract "sen="

        // Extract Athena geometry ID (a_mod)
        std::string a_mod_str = segments[2].substr(segments[2].find(',') + 1);
        Identifier a_mod = Identifier(std::stoull(a_mod_str, nullptr, 16));  // Convert to int

        // Track count of a_mod occurrences
        a_mod_count[a_mod]++;

        // Store in the map
        AthenaGeometryMap[a_mod] = std::make_tuple(volume, layer, sensor);
    }
    map_file.close();

    ATH_MSG_INFO("Loaded Acts-to-Athena mapping with " << AthenaGeometryMap.size() << " entries");

    // Saving digi config file
    nlohmann::json json_output;
    json_output["acts-geometry-hierarchy-map"] = {
        {"format-version", 0},
        {"value-identifier", "digitization-configuration"}
    };
    json_output["entries"] = nlohmann::json::array();

    // Also saving all strip information to a csv file
    std::string csv_filename = "/eos/user/j/jlai/g200/gpu/G-200/traccc-athena/run/ITk_strip_info.csv";
    std::ofstream csv_file(csv_filename);
    
    if (!csv_file.is_open()) {
        ATH_MSG_ERROR("Failed to open CSV file for writing: " << csv_filename);
        return;
    }

    // Write CSV Header
    csv_file << "Detray_Geometry_ID,Athena_Geometry_ID,Volume,Layer,Sensor,Side,"
                "Placement_X,Placement_Y,Placement_Z,"
                "Rotation_00,Rotation_01,Rotation_02,Rotation_10,Rotation_11,Rotation_12,Rotation_20,Rotation_21,Rotation_22,"
                "Reference_X,Reference_Y,Pitch_X,Pitch_Y,Module_Width,Module_Length,Rows,Columns,Dimensions\n";    


    for (const auto& surface : m_detector.surfaces()) {
      auto sf = detray::tracking_surface{m_detector, surface};
      if (!surface.is_sensitive()){continue;}
      
      auto id = surface.barcode().value();
      const auto geo_id = surface.source;
      Acts::GeometryIdentifier acts_geom_id{geo_id};

      Identifier a_mod = Identifier(m_DetrayToAtlasMap[id]);

      // If a_mod is repeated more than once, skip and print
      if (a_mod_count[a_mod] > 2) {
        ATH_MSG_WARNING("Skipping duplicate Athena ID: " << a_mod << " which appeared " << a_mod_count[a_mod] << " times.");
        continue; }

      // If it appears exactly once, increment sensor count by 1
      if (processed_a_mods.find(a_mod) != processed_a_mods.end()) {
          AthenaGeometryMap[a_mod] = std::make_tuple(std::get<0>(AthenaGeometryMap[a_mod]), 
                                                      std::get<1>(AthenaGeometryMap[a_mod]), 
                                                      std::get<2>(AthenaGeometryMap[a_mod]) - 1); }
      processed_a_mods.insert(a_mod);

      // if (m_atlasModuleInfo.find(a_mod) == m_atlasModuleInfo.end()) {
      //  ATH_MSG_WARNING("Cannot find this module " << a_mod )
      // continue;
      moduleInfo thismod = m_atlasModuleInfo[a_mod];

      // // Check if side-aware ID exists in m_atlasModuleInfo
      // if (m_atlasModuleInfo.find(sideAwareID) != m_atlasModuleInfo.end()) {
      //     thismod = m_atlasModuleInfo[sideAwareID]; 
      //     if (scount <= 10) {
      //     ATH_MSG_INFO("Retrieving module " << a_mod << 
      //                  " sideAwareID " << sideAwareID <<
      //                 " with Side: " << thismod.side);}
      // }

      // if (scount <= 50) {
      //   ATH_MSG_INFO("Retrieving module " << a_mod << " with Side: " << thismod.side);
      // }

      // Fetch volume, layer, sensor from AthenaGeometryMap
      // Declare local variables for volume, layer, sensor
      int volume = -1, layer = -1, sensor = -1;
      if (AthenaGeometryMap.find(a_mod) != AthenaGeometryMap.end()) {
        std::tie(volume, layer, sensor) = AthenaGeometryMap[a_mod];
      } else {
          ATH_MSG_WARNING("No mapping found for Athena Geometry ID: " << a_mod);
      }

      int side = thismod.side;


      if (thismod.pixel) {         
        if (pcount <= 10) { // Print the Pixel data but do NOT save in m_dd
        ATH_MSG_INFO("Pixel - Detray geometry_id: " << id << " Athena geometry_id: " << a_mod << " vol: " << volume << " lay: " << layer << " sen: " << sensor 
                    << " side: " << side
                    << ", placement: (" 
                    << sf.transform(context).translation()[0] << ", "
                    << sf.transform(context).translation()[1] << ", "
                    << sf.transform(context).translation()[2] << ")"
                    << ", rotation: (" 
                    << sf.transform(context).rotation()[0] << ", "
                    << sf.transform(context).rotation()[1] << ", "
                    << sf.transform(context).rotation()[2] << ")"
                    << ", threshold: 0"
                    << ", reference_x: " << (-thismod.module_width * 0.5)
                    << ", reference_y: " << (-thismod.module_length * 0.5)
                    << ", pitch_x: " << (thismod.module_width / thismod.rows)
                    << ", pitch_y: " << (thismod.module_length / thismod.columns)
                    << ", module width: " << (thismod.module_width)
                    << ", module length: " << (thismod.module_length)
                    << ", module rows: " << (thismod.rows)
                    << ", module columns: " << (thismod.columns)
                    << ", dimensions: 2"); // Pixels are always dim = 2
        pcount++;}
        // continue; 
        json_output["entries"].push_back({
          {"layer", layer},
          {"sensitive", sensor},
          {"volume", volume},
          {"value", {
              {"geometric", {
                  {"segmentation", {
                      {"binningdata", {
                          {
                              {"bins", thismod.rows},
                              {"max", thismod.module_width * 0.5},
                              {"min", -thismod.module_width * 0.5},
                              {"option", "open"}, // change back to open
                              {"type", "equidistant"},
                              {"value", "binX"}
                          },
                          {
                              {"bins", thismod.columns}, // for pixel, this must be column 
                              {"max", thismod.module_length * 0.5},
                              {"min", -thismod.module_length * 0.5},
                              {"option", "open"},
                              {"type", "equidistant"},
                              {"value", "binY"}
                          }
                      }}
                    }}
                  }}
                }}
              });

      }
      else {
        if (scount <= 10) {
          ATH_MSG_INFO("Strip - Detray geometry_id: " << id << " Athena geometry_id: " << a_mod << " vol: " << volume << " lay: " << layer << " sen: " << sensor 
                      << " side: " << side
                      << ", placement: (" 
                      << sf.transform(context).translation()[0] << ", "
                      << sf.transform(context).translation()[1] << ", "
                      << sf.transform(context).translation()[2] << ")"
                      << ", rotation: (" 
                      << sf.transform(context).rotation()[0] << ", "
                      << sf.transform(context).rotation()[1] << ", "
                      << sf.transform(context).rotation()[2] << ")"
                      << ", threshold: 0"
                      << ", reference_x: " << (-thismod.module_width * 0.5)
                      << ", reference_y: " << (-thismod.module_length * 0.5)
                      << ", pitch_x: " << (thismod.module_width / thismod.rows)
                      << ", pitch_y: " << (thismod.module_length / thismod.columns)
                      << ", module width: " << (thismod.module_width)
                      << ", module length: " << (thismod.module_length)
                      << ", module rows: " << (thismod.rows)
                      << ", module columns: " << (thismod.columns)
                      << ", dimensions: 1"); 
          scount++;
        }

        json_output["entries"].push_back({
          {"layer", layer},
          {"sensitive", sensor},
          {"volume", volume},
          {"value", {
              {"geometric", {
                  {"segmentation", {
                      {"binningdata", {
                          {
                              {"bins", thismod.rows},
                              {"max", thismod.module_width * 0.5},
                              {"min", -thismod.module_width * 0.5},
                              {"option", "open"}, // change back to open
                              {"type", "equidistant"},
                              {"value", "binX"}
                          },
                          {
                              {"bins", 1}, // for strip, bin must be 1
                              {"max", thismod.module_length * 0.5},
                              {"min", -thismod.module_length * 0.5},
                              {"option", "open"},
                              {"type", "equidistant"},
                              {"value", "binY"}
                          }
                      }}
                    }}
                  }}
                }}
              });

        csv_file << id << "," << a_mod << "," << volume << "," << layer << "," << sensor << ","
                << side << ","
                << sf.transform(context).translation()[0] << "," 
                << sf.transform(context).translation()[1] << "," 
                << sf.transform(context).translation()[2] << ","
                << sf.transform(context).rotation()[0][0] << "," 
                << sf.transform(context).rotation()[0][1] << "," 
                << sf.transform(context).rotation()[0][2] << ","
                << sf.transform(context).rotation()[1][0] << "," 
                << sf.transform(context).rotation()[1][1] << "," 
                << sf.transform(context).rotation()[1][2] << ","
                << sf.transform(context).rotation()[2][0] << "," 
                << sf.transform(context).rotation()[2][1] << "," 
                << sf.transform(context).rotation()[2][2] << ","
                << (-thismod.module_width * 0.5) << "," 
                << (-thismod.module_length * 0.5) << ","
                << (thismod.module_width / thismod.rows) << "," 
                << (thismod.module_length / thismod.columns) << ","
                << (thismod.module_width) << "," 
                << (thismod.module_length) << ","
                << (thismod.rows) << "," 
                << (thismod.columns) << ","
                << 1 << "\n";

        m_dd.resize(m_dd.size() + 1);
        m_dd.geometry_id().back() = detray::geometry::barcode{surface.barcode()};
        m_dd.placement().back() = sf.transform(context);
        m_dd.acts_geometry_id().back() = id;//surface.source;


        // Set a hard-coded threshold for which cells should be considered for
        // clusterization on this module / surface.
        m_dd.threshold().back() = 0.f;
        m_dd.reference_x().back() = (-thismod.module_width*0.5);
        m_dd.reference_y().back() = (-thismod.module_length*0.5);
        m_dd.pitch_y().back() = thismod.module_length/(thismod.columns);
        m_dd.pitch_x().back() = thismod.module_width/(thismod.rows);
        m_dd.dimensions().back() = 1;
        
        // if(thismod.pixel){m_dd.dimensions().back() = 2;}
        // else{m_dd.dimensions().back() = 1;}
      }
    }

    // Save JSON to file
    std::string output_filename = "/eos/user/j/jlai/g200/gpu/G-200/traccc-athena/run/ITk_digitization_config_with_strips.json";
    try {
        std::ofstream json_file(output_filename);
        json_file << json_output.dump(4);  // Pretty print with indentation = 4 spaces
        json_file.close();
        ATH_MSG_INFO("Successfully saved digitization config to " << output_filename);
    } catch (const std::exception &e) {
        ATH_MSG_ERROR("Error saving digitization JSON file: " << e.what());
    }

    // Close CSV file
    csv_file.close();
    ATH_MSG_INFO("Successfully saved strip data to CSV: " << csv_filename);

    // ATH_MSG_INFO("Printing first 25 entries of strip:");

    // int count = 0;
    // for (size_t i = 0; i < m_dd.size(); ++i) {
    //     if (count >= 25) break;  // Print only the first 10 entries

    //         // Retrieve moduleInfo again using the geometry ID
    //     int atlas_id = m_DetrayToAtlasMap[m_dd.geometry_id()[i].value()];
    //     moduleInfo thismod = m_atlasModuleInfo[atlas_id]; // Retrieve from the map
    
    //     ATH_MSG_INFO("Strip " << count + 1 
    //                  << " - geometry_id: " << m_dd.geometry_id()[i].value()
    //                  << ", placement: (" 
    //                  << m_dd.placement()[i].translation()[0] << ", "  // Using array indexing
    //                  << m_dd.placement()[i].translation()[1] << ", "
    //                  << m_dd.placement()[i].translation()[2] << ")"
    //                  << ", rotation: (" 
    //                  << m_dd.placement()[i].rotation()[0] << ", "  // If rotation exists, use indexing
    //                  << m_dd.placement()[i].rotation()[1] << ", "
    //                  << m_dd.placement()[i].rotation()[2] << ")"
    //                  << ", acts_geometry_id: " << m_dd.acts_geometry_id()[i]
    //                  << ", threshold: " << m_dd.threshold()[i]
    //                  << ", reference_x: " << m_dd.reference_x()[i]
    //                  << ", reference_y: " << m_dd.reference_y()[i]
    //                  << ", pitch_x: " << m_dd.pitch_x()[i]
    //                  << ", pitch_y: " << m_dd.pitch_y()[i]
    //                  << ", module width: " << m_dd.pitch_x()[i] * thismod.rows
    //                  << ", module length: " << m_dd.pitch_y()[i] * thismod.columns
    //                  << ", module rows: " << (thismod.rows)
    //                  << ", module columns: " << (thismod.columns)
    //                  << ", dimensions: " << static_cast<int>(m_dd.dimensions()[i])
                    //  << ", vol, lay, sen " << (thismod.volume) << (thismod.layer) << (thismod.sensitive) );
    //     count++;
    // }
    
    // dump_digi_cfg(m_dd);

    // std::vector<std::pair<Acts::GeometryIdentifier, traccc::module_digitization_config>> module_digi;

    // for(unsigned int i = 0; i < m_dd.acts_geometry_id().size(); i++){

    //     traccc::module_digitization_config config;
    //     Acts::BinUtility bUtility;
    //     bUtility += Acts::BinUtility(
    //     (-2*m_dd.reference_x().at(i))/((m_dd.pitch_x().at(i))), (m_dd.reference_x().at(i)), (-m_dd.reference_x().at(i)), Acts::open, Acts::BinningValue::binX // X binning
    //     );
    //     bUtility += Acts::BinUtility(
    //     (-2*m_dd.reference_y().at(i))/((m_dd.pitch_y().at(i))), (m_dd.reference_y().at(i)), (-m_dd.reference_y().at(i)), Acts::open, Acts::BinningValue::binY // X binning
    //     );
    //     config.segmentation = bUtility;
    //     module_digi.emplace_back(m_dd.acts_geometry_id().at(i),config);

    // }

    // const traccc::digitization_config digi_conf(module_digi);
    // for (size_t i = 0; i < module_digi.size(); ++i) {
    //   if (std::isnan(m_dd.pitch_x().at(i)) || std::isnan(m_dd.pitch_y().at(i)) ||
    //       std::isnan(m_dd.reference_x().at(i)) || std::isnan(m_dd.reference_y().at(i))) {
    //       ATH_MSG_ERROR("NaN detected in module " << module_digi[i].first);
    //       continue; } 
    //   if (m_dd.pitch_x().at(i) <= 0 || m_dd.pitch_y().at(i) <= 0) {
    //     ATH_MSG_ERROR("Invalid pitch_x or pitch_y for module " << module_digi[i].first
    //                   << " - pitch_x: " << m_dd.pitch_x().at(i) 
    //                   << ", pitch_y: " << m_dd.pitch_y().at(i));
    //     continue; } }
  
    // ATH_MSG_INFO("Writing digitization config with " << module_digi.size() << " modules");

    // try {
    //   traccc::io::write("/eos/user/j/jlai/g200/gpu/G-200/traccc-athena/run/ITk_digitization_config_with_strips.json", 
    //                      traccc::data_format::json, digi_conf);
    //     ATH_MSG_INFO("Successfully wrote ITk_digitization_config_with_strips.json");
    // } catch (const std::exception &e) {
    //     ATH_MSG_ERROR("Failed to write JSON: " << e.what());
    // } catch (...) {
    //     ATH_MSG_ERROR("Unknown error occurred while writing JSON");
    // }

    ////////////////////////////////////////////////////////////////////////
  

    // traccc::io::write("ITk_digitization_config_with_strips.json",traccc::data_format::json,digi_conf);


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

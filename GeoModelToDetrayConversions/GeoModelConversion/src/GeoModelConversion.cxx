#include "GeoModelConversion/GeoModelConversion.h"
#include "ActsGeometry/ActsDetectorElement.h"
#include "Acts/Plugins/GeoModel/GeoModelReader.hpp"
#include "Acts/Plugins/GeoModel/GeoModelBlueprintCreater.hpp"
#include "Acts/Definitions/Units.hpp"

#include "InDetReadoutGeometry/SiDetectorElement.h"
#include "InDetIdentifier/SCT_ID.h"

#include "InDetReadoutGeometry/SiDetectorDesign.h"
#include "PixelReadoutGeometry/PixelModuleDesign.h"
#include "SCT_ReadoutGeometry/SCT_BarrelModuleSideDesign.h"
#include "SCT_ReadoutGeometry/SCT_ForwardModuleSideDesign.h"
#include "SCT_ReadoutGeometry/StripStereoAnnulusDesign.h"

#include "Acts/Visualization/ObjVisualization3D.hpp"

#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"
#include "Acts/Detector/interface/IDetectorBuilder.hpp"

#include <fstream>

GeoModelConversion::GeoModelConversion ( const std::string &n, ISvcLocator *pSvcLoc )
    : AthAlgorithm ( n, pSvcLoc )
{
}


StatusCode GeoModelConversion::initialize(){

	ATH_CHECK(detStore()->retrieve(m_pixelId, "PixelID"));
	ATH_CHECK(detStore()->retrieve(m_stripId, "SCT_ID"));
	ATH_CHECK(m_pixelDetEleCollKey.initialize());
	ATH_CHECK(m_stripDetEleCollKey.initialize());

	if (!m_trackingGeometryTool.empty()) {
		ATH_CHECK(m_trackingGeometryTool.retrieve());
		m_trackingGeometry = m_trackingGeometryTool->trackingGeometry();

    	m_trackingGeometry->visitSurfaces([&](const Acts::Surface *surface) {
			// find acts surface with associated detector element ID
			if (!surface) return;
			m_allsurfaces.push_back(std::const_pointer_cast<Acts::Surface>(surface->getSharedPtr()));
			const auto *actsElement = dynamic_cast<const ActsDetectorElement *>(surface->associatedDetectorElement());
			if (!actsElement) return;
			const auto *geoElement = actsElement->upstreamDetectorElement();
			const auto *detElem = dynamic_cast<const InDetDD::SiDetectorElement*>(geoElement);
			if (!geoElement || !detElem) return;
			m_surfaces.push_back(std::const_pointer_cast<Acts::Surface>(surface->getSharedPtr()));

   		});
	}

	return StatusCode::SUCCESS;
}

void GeoModelConversion::load_detray_detector(){

	// Read the detector.
	detray::io::detector_reader_config reader_cfg{};
	reader_cfg.add_file("./ITk_DetectorBuilder_geometry.json");
	auto [host_det, _] = detray::io::read_detector<detector_t>(host_mr, reader_cfg);
	m_detector = std::move(host_det);
}

void GeoModelConversion::build_detray_detector(){

	std:: string db;
	auto context = Gaudi::Hive::currentContext();
	Acts::GeometryContext geoContext = m_trackingGeometryTool->getGeometryContext(context).context();

	std::string sqliteDbName = "/eos/user/n/nribaric/ITk_data/ITk_geometry_data/ITk_SQLite_acts.db";
	auto gmTree = Acts::GeoModelReader::readFromDb(sqliteDbName);

	auto gmBlueprintConfig = Acts::GeoModelBlueprintCreater::Config();
	gmBlueprintConfig.detectorSurfaces = m_surfaces;
	gmBlueprintConfig.kdtBinning = {Acts::AxisDirection::AxisZ, Acts::AxisDirection::AxisR};

	auto gmBlueprintOptions = Acts::GeoModelBlueprintCreater::Options();
	gmBlueprintOptions.table = "ActsBlueprint";
	gmBlueprintOptions.topEntry = "ITk";

	auto logger = Acts::getDefaultLogger("GeoModelBlueprintCreater", Acts::Logging::INFO);
	auto gmBlueprintCreater = Acts::GeoModelBlueprintCreater(gmBlueprintConfig, std::move(logger));
	auto gmBlueprint = gmBlueprintCreater.create(geoContext, gmTree, gmBlueprintOptions);

	auto gmCylindricalBuilder = std::make_shared<Acts::Experimental::CylindricalContainerBuilder>(*gmBlueprint.topNode,Acts::Logging::VERBOSE);

	std::ofstream fout;
	fout.open("Itk.dot");
	Acts::Experimental::detail::BlueprintDrawer::dotStream(fout, *gmBlueprint.topNode);
	fout.close();

	auto gmGeoIdConfig = Acts::Experimental::GeometryIdGenerator::Config();
	auto logger_3 = Acts::getDefaultLogger("GeoModelBlueprintCreater3", Acts::Logging::INFO);
	auto gmGeoIdGenerator = std::make_shared<Acts::Experimental::GeometryIdGenerator>(gmGeoIdConfig, std::move(logger_3));

	auto gmDetectorConfig = Acts::Experimental::DetectorBuilder::Config();
	gmDetectorConfig.name = "ITk_DetectorBuilder";
	gmDetectorConfig.builder = gmCylindricalBuilder;
	gmDetectorConfig.geoIdGenerator = gmGeoIdGenerator;
	gmDetectorConfig.materialDecorator = nullptr; // WE WILL FILL THIS LATER
	gmDetectorConfig.auxiliary = "GeoModel based Acts::Detector from Db";

	auto logger_2 = Acts::getDefaultLogger("GeoModelBlueprintCreater2", Acts::Logging::INFO);
	auto gmDetectorBuilder = Acts::Experimental::DetectorBuilder(gmDetectorConfig, std::move(logger_2));

	auto detector = gmDetectorBuilder.construct(geoContext);
	auto opts = Acts::DetrayConverter::Options();
	opts.writeToJson = true;
	auto detrayConverter = Acts::DetrayConverter();
	m_detector = detrayConverter.convert(geoContext, *detector, host_mr, opts);


}

void GeoModelConversion::build_geometry_maps(){

	int nPix = 0;
	int nStrip = 0;
	int nF = 0;
	int nE = 0;
	int nM = 0;

	ATH_MSG_INFO("Making ACTS to Athena map");

	for(const auto& surface : m_surfaces){
		const auto *actsElement = dynamic_cast<const ActsDetectorElement *>(surface->associatedDetectorElement());
		if (!actsElement){continue;}
		auto acts_id = actsElement->identify();
		const auto geo_id = (surface.get())->geometryId();

		const auto *geoElement = actsElement->upstreamDetectorElement();
		const auto *detElem = dynamic_cast<const InDetDD::SiDetectorElement*>(geoElement);

		if (!geoElement){
			nE++;
			ATH_MSG_DEBUG("We have this ACTS detector element, but it does not have a corresponding Athena detector element: " << geo_id);
			continue;
		}
		if(!detElem){
			nM++;
			ATH_MSG_DEBUG("We have this Athena detector element, but it is not a SiDetector element: " << geo_id);
			continue;
		}

		Identifier athenaID;

		if(detElem->isPixel()){
			nPix++;
			athenaID = detElem->identify();

		} else if(detElem->isSCT()){
			nStrip++;
			athenaID = m_stripId->module_id(detElem->identify());
			const IdentifierHash Strip_ModuleHash = m_stripId->wafer_hash(athenaID);

			// Extract the correct Strip_ModuleID
			int side = m_stripId->side(detElem->identify());
			athenaID = m_stripId->wafer_id(Strip_ModuleHash + side);
		}else{
			ATH_MSG_INFO("Not a SiDetElement " << geo_id);
		}

		if(m_actsToGeoModel.find(geo_id) != m_actsToGeoModel.end()){
			ATH_MSG_INFO("We already have this ACTS ID in map: " << geo_id << " , " << athenaID);
		}
		m_actsToGeoModel[geo_id] = athenaID;


	}

	ATH_MSG_INFO("Not done " << nF << " modules, out of " << m_surfaces.size() << " surfaces used");
	ATH_MSG_INFO("Not found " << nE << " SiElements, out of those " << nM << " modules not found");

	typename detector_t::geometry_context context{};
	ATH_MSG_INFO("Making Detray to ACTS map");
	nF = 0;
	int nFs = 0;

	ATH_MSG_INFO("Made " << m_detector.surfaces().size() << " surfaces in detray");

	std::vector<Acts::GeometryIdentifier> detray_surfaces;

	std::ofstream AtoD_map_file("./athenaIdentifierToDetrayMap.txt",std::ofstream::out);
	std::ofstream AtoA_map_file("./actsToAthenaIdentifierMap.txt",std::ofstream::out);

	for (const auto& surface : m_detector.surfaces()) {

		// Acts geometry identifier(s) for the surface.
        const auto geo_id = surface.source;
        const Acts::GeometryIdentifier acts_geom_id{geo_id};
		auto sf = detray::tracking_surface{m_detector, surface};
		auto detray_id = sf.barcode().value();

		detray_surfaces.push_back(acts_geom_id);

		if(m_actsToGeoModel.contains(acts_geom_id)){
			auto athena_id = m_actsToGeoModel.at(acts_geom_id);
			AtoD_map_file << athena_id << "," << detray_id << "\n";
			m_DetrayToAtlasMap[detray_id] = athena_id;
			AtoA_map_file << acts_geom_id << "," << athena_id << "\n";

		}
		else{
			nF++;
			ATH_MSG_DEBUG("we did not save key " << acts_geom_id);
			if (surface.is_sensitive()){continue;}
			nFs++;
			ATH_MSG_DEBUG("found this passive surface in detray: " << acts_geom_id);

		}

	}

	AtoA_map_file.close();
	AtoD_map_file.close();


	int n_un = 0;

	for(const auto& surface : m_surfaces){

		const auto geo_id = (surface.get())->geometryId();

		if(std::find(detray_surfaces.begin(), detray_surfaces.end(), geo_id) == detray_surfaces.end()){

			ATH_MSG_INFO("This ACTS surface ID was not found in Detray: " << geo_id);
			n_un++;
		}
	}

	ATH_MSG_INFO("Not matched " << nF << " detray surfaces, " << nFs << " passives.");

}

void GeoModelConversion::createAthenaMap() {

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

	int scount = 0;
	// Saving digi config file
	nlohmann::json json_output;
	json_output["acts-geometry-hierarchy-map"] = {
		{"format-version", 0},
		{"value-identifier", "digitization-configuration"}
	};
	json_output["entries"] = nlohmann::json::array();

	SG::ReadCondHandle<InDetDD::SiDetectorElementCollection> stripDetEleHandle(m_stripDetEleCollKey);
	const InDetDD::SiDetectorElementCollection* strip_elements(*stripDetEleHandle);
	if (not stripDetEleHandle.isValid() or strip_elements==nullptr) {
		ATH_MSG_FATAL(m_stripDetEleCollKey.fullKey() << " is not available.");
	}

	for (const InDetDD::SiDetectorElement *element: *strip_elements) {
		const Identifier strip_moduleID = m_stripId->module_id(element->identify()); //from wafer id to module id
		const IdentifierHash Strip_ModuleHash = m_stripId->wafer_hash(strip_moduleID);

		// Extract the correct Strip_ModuleID
		int side = m_stripId->side(element->identify());
		const Identifier Strip_ModuleID = m_stripId->wafer_id(Strip_ModuleHash + side);

		// Extract side information (0 = inner, 1 = outer)
		int layer_disk = m_stripId->layer_disk(element->identify());
		float eta_module = m_stripId->eta_module(element->identify());
		float phi_module = m_stripId->phi_module(element->identify());
		int barrel_ec = m_stripId->barrel_ec(element->identify());
		double stripPitch = 0.0;
		int cell = 0.0;
		double length = 0.0;
		double width = 0.0;

		if(scount < 10){
			ATH_MSG_INFO("Strip side: " << side);
			ATH_MSG_INFO("Strip module ID: " << Strip_ModuleID << ",  " << strip_moduleID);
			ATH_MSG_INFO("Strip B/EC: " << barrel_ec << ", layer: " << layer_disk << ", eta: " << eta_module << ", phi: " << phi_module);
			scount++;
		}

		++nStrip;
		ATH_MSG_VERBOSE( "Strip module " << nStrip );
		moduleInfo thismod;

		thismod.center0 = 0.0;
		thismod.center1 = 0.0;
		thismod.center2 = 0.0;

		if (m_stripId->barrel_ec(Strip_ModuleID) == 0) {
			++nbar;
			const InDetDD::SCT_BarrelModuleSideDesign* s_design = (static_cast<const InDetDD::SCT_BarrelModuleSideDesign*>(&element->design()));
			auto boundsType = element->bounds().type();
			thismod.pixel = false;
			thismod.side = side; // store side information

			if (nbar <= 15) {
				ATH_MSG_INFO("barrel: " << m_stripId->barrel_ec(Strip_ModuleID) << " Module ID: " << Strip_ModuleID <<
				", Side: " << side << ", Side: " << thismod.side);
			}


			if (boundsType == Trk::SurfaceBounds::Rectangle) {
				++nRect1;
				ATH_MSG_DEBUG("Strip module in barrel with ID " << Strip_ModuleID << " is designed like a rectangle");

				thismod.module_width = s_design->width();
				thismod.module_length = s_design->length();
				thismod.rows = s_design->cells(); // same as s_design->width() / s_design->stripPitch();

				// To save in df_map csv file
				width = s_design->width();
				length = s_design->length();
				cell = s_design->cells();
				stripPitch = s_design->stripPitch();

				if (nRect1 < 10) {
					ATH_MSG_INFO(" nRect1, s_design->width(): " << s_design->width() <<
								" length(): " << s_design->length() <<
								" cells(): " << s_design->cells() <<
								" stripPitch(): " << s_design->stripPitch() << "\n" <<
								" stripLength(): " << s_design->stripLength() <<
								" phiStripPatternCentre(): " << s_design->phiStripPatternCentre() <<
								" etaStripPatternCentre(): " << s_design->etaStripPatternCentre() <<
								" phiPitch(): " << s_design->phiPitch() <<
								" etaPitch(): " << s_design->etaPitch() );
				}

				// Extract side information (0 = inner, 1 = outer)
				int layer_disk = m_stripId->layer_disk(element->identify());
				float eta_module = m_stripId->eta_module(element->identify());
				float phi_module = m_stripId->phi_module(element->identify());
				int barrel_ec = m_stripId->barrel_ec(element->identify());

				if (s_design->cells() != s_design->width() / s_design->stripPitch()) {
					ATH_MSG_INFO(" nRect1, problem, width(): " << s_design->width() <<
					" length(): " << s_design->length() <<
					" cells(): " << s_design->cells() <<
					" stripPitch(): " << s_design->stripPitch() <<
					" stripLength(): " << s_design->stripLength() );
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

		} else { // EC modules

			++nnonbar;
			const InDetDD::StripStereoAnnulusDesign* annulus_design = (static_cast<const InDetDD::StripStereoAnnulusDesign*>(&element->design()));
			double m_stereo = annulus_design->stereo();
			double m_waferCentreR = annulus_design->waferCentreR();
			double m_lengthBF = 2. * m_waferCentreR * std::sin(m_stereo / 2.);

			const InDetDD::SiCellId annulus_id = element->cellIdFromIdentifier(Strip_ModuleID);
			double phiPitch = annulus_design->phiPitch(annulus_id);
			double m_pitch_row = annulus_design->phiPitchPhi(annulus_id);
			double radius = annulus_design->centreR();

			int m_nStrips0 = annulus_design->diodesInRow(0.);
			double max_phi = (m_nStrips0) * m_pitch_row;
			thismod.module_width = max_phi;
			thismod.module_length = radius;
			thismod.rows = m_nStrips0;

			const InDetDD::SCT_ForwardModuleSideDesign* s_design = (static_cast<const InDetDD::SCT_ForwardModuleSideDesign*>(&element->design()));
			auto boundsType = element->bounds().type();
			thismod.pixel = false;
			thismod.side = side;

			// Testing Forward Module
			if (nnonbar <= 15) {
				ATH_MSG_INFO("barrel: " << m_stripId->barrel_ec(Strip_ModuleID) << " Module ID: " << Strip_ModuleID <<
				", Side: " << side << ", Side: " << thismod.side);
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
				stripPitch = s_design->stripPitch();

			}

		}

		const InDetDD::SiDetectorElement *module = strip_elements->getDetectorElement(Strip_ModuleHash+side);
		std::vector<const Trk::Surface*> module_surfaces = module->surfaces();

		thismod.center0 = module->center()[0];
		thismod.center1 = module->center()[1];
		thismod.center2 = module->center()[2];

		m_atlasModuleInfo[Strip_ModuleID] = thismod;

		if (nStrip < 10) {
			ATH_MSG_INFO( "Human Readable ID: " << Strip_ModuleID << " center2: " << thismod.center2 << " surface: " << module->surfaces() << " wafer_id: " << m_stripId->wafer_id(Strip_ModuleHash+side));
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


			const InDetDD::SiDetectorElement *module = pixel_elements->getDetectorElement(Pixel_ModuleHash);
			std::vector<const Trk::Surface*> module_surfaces = module->surfaces();

			++nPixel;
			ATH_MSG_DEBUG( "Pixel module " << nPixel );

			// Write out Visualization Lookup Tree
			//m_AthenaHashedID = Pixel_ModuleID.get_identifier32().get_compact();

			moduleInfo thismod;
			const InDetDD::PixelModuleDesign *p_design = static_cast<const InDetDD::PixelModuleDesign*>(&element->design());
			auto boundsType = module->bounds().type();
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

}

void GeoModelConversion::fill_digi_info(){

    ATH_MSG_DEBUG("Filling detector description");
    typename detector_t::geometry_context context{};


    int pcount = 0;
    int scount = 0;
    std::set<Identifier> processed_athena_ids;  // set to track processed Athena Geometry IDs
    std::map<Identifier, int> athena_id_count;

    // Saving digi config file
    nlohmann::json json_output;
    json_output["acts-geometry-hierarchy-map"] = {
        {"format-version", 0},
        {"value-identifier", "digitization-configuration"}
    };
    json_output["entries"] = nlohmann::json::array();

	int nFoundMod = 0;
	int nFound = 0;
	int duplicate = 0;
    for (const auto& surface : m_detector.surfaces()) {
      auto sf = detray::tracking_surface{m_detector, surface};
      if (!surface.is_sensitive()){continue;}

      auto id = surface.barcode().value();
      const auto geo_id = surface.source;
      Acts::GeometryIdentifier acts_geom_id{geo_id};


      if(!m_DetrayToAtlasMap.contains(id)){
        if (pcount + scount <= 10) {ATH_MSG_ERROR("Could not find this id: " << id << " in detray to atlas map.");}
        nFound++;
      }

      Identifier athena_id = Identifier(m_DetrayToAtlasMap[id]);
	  if(processed_athena_ids.contains(athena_id)){
		duplicate++;
		athena_id_count[athena_id] += 1;
	  }

      if (pcount + scount <= 10) {
        ATH_MSG_INFO("geo_id: " << geo_id << " acts_id: " << acts_geom_id << " id: " << id << " athena_id: " << athena_id);
      }

      std::ostringstream oss;
      oss << acts_geom_id;
      std::stringstream ss(oss.str());

      std::vector<std::string> parts;
      std::string part;
      while (std::getline(ss, part, '|')) {
          parts.push_back(part);
      }

      int volume = std::stoi(parts[0].substr(4));    // Remove "vol="
      int layer = std::stoi(parts[1].substr(4));     // Remove "lay="
      int sensor = std::stoi(parts[2].substr(4));    // Remove "sen="

      // Match the correct athena_id
      double placement_x = sf.transform(context).translation()[0];
      double placement_y = sf.transform(context).translation()[1];
      double placement_z = sf.transform(context).translation()[2];


      // Make the correct athena_id
      if (volume == 22 || volume == 23 || volume == 24) {
        if (sensor % 2 == 0) { // even
          athena_id = Identifier(athena_id.get_compact() | 0x0008000000000);
          if(scount < 10){ATH_MSG_INFO(" athena_id even side: " << athena_id);}
        }
      }


      if(!m_atlasModuleInfo.contains(athena_id)){
        if (pcount + scount <= 10) {ATH_MSG_ERROR("Could not find module for this id: " << athena_id << " in module map.");}
        nFoundMod++;
      }

      moduleInfo thismod = m_atlasModuleInfo[athena_id];

      // Check if the athena_id is correct
      if (!thismod.pixel) {
        if (volume == 22 || volume == 24) { // annulus
          if (abs(placement_z - thismod.center2) > 0.1) {
            ATH_MSG_INFO("Annulus wrong!" << " placement_z: " << placement_z << " center2: " << thismod.center2);
          }
        } else if (volume == 23) { // barrel
          if ((abs(placement_x - thismod.center0) > 0.1) || (abs(placement_y - thismod.center1) > 0.1)) {
            ATH_MSG_INFO("Barrel wrong!" << " placement_x: " << placement_x << " center0: " << thismod.center0);
          }
        }
      }

      int side = thismod.side;

      if (thismod.pixel) {
        pcount++;
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
        scount++;
        if (volume == 22 || volume == 24) { // annulus

          // Convert x, y to radius, angle
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
                              {"max", thismod.module_width * 0.5}, // max theta
                              {"min", -thismod.module_width * 0.5}, // min theta
                              {"option", "open"},
                              {"type", "equidistant"},
                              {"value", "binX"}
                          },
                          {
                              {"bins", 1}, // for strip, bin must be 1
                              {"max", thismod.module_length+0.1}, // radius
                              {"min", thismod.module_length-0.1}, // radius
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

        if (volume == 23) { // barrel
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
        }
      }
    }

	ATH_MSG_INFO("Could not find module info for " << nFoundMod << " modules.");
	ATH_MSG_INFO("Could not find maps for " << nFound << " modules.");
	ATH_MSG_INFO("Made duplicates for " << duplicate << " modules.");

    // Save JSON to file
    std::string output_filename = "ITk_digitization_config_with_strips.json";
    try {
        std::ofstream json_file(output_filename);
        json_file << json_output.dump(4);  // Pretty print with indentation = 4 spaces
        json_file.close();
        ATH_MSG_INFO("Successfully saved digitization config to " << output_filename);
    } catch (const std::exception &e) {
        ATH_MSG_ERROR("Error saving digitization JSON file: " << e.what());
    }

}

StatusCode GeoModelConversion::execute(){


	// first lets build and dump the detray detector into a file
	build_detray_detector();

	// or you can just load it if you only want the mapping to be made
	//load_detray_detector();

	// now let's create acts - detray - geomodel maps
	// and the digitisation config
	build_geometry_maps();
	createAthenaMap();
	fill_digi_info();

	return StatusCode::SUCCESS;
}

StatusCode GeoModelConversion::finalize(){

        return StatusCode::SUCCESS;
}

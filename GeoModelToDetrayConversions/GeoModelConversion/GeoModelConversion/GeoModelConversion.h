#include "AthenaBaseComps/AthAlgorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include "StoreGate/DataHandle.h"

// ACTS includes
#include "ActsGeometryInterfaces/IActsTrackingGeometryTool.h"
//#include "ActsEventCnv/IActsToTrkConverterTool.h"
#include "ActsEvent/TrackContainerHandlesHelper.h"
#include "Acts/Detector/Blueprint.hpp"
#include "Acts/Detector/DetectorBuilder.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/GeoModel/GeoModelReader.hpp"
#include "Acts/Plugins/GeoModel/GeoModelReader.hpp"
#include "Acts/Plugins/GeoModel/GeoModelBlueprintCreater.hpp"
#include "Acts/Plugins/Detray/DetrayConverter.hpp"
#include "Acts/Detector/CylindricalContainerBuilder.hpp"
#include "Acts/Detector/detail/BlueprintDrawer.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"

#include "detray/io/frontend/detector_reader.hpp"
#include "detray/core/detector.hpp"
#include "detray/geometry/tracking_surface.hpp"

#include "InDetReadoutGeometry/SiDetectorDesign.h"
#include "PixelReadoutGeometry/PixelModuleDesign.h"
#include "SCT_ReadoutGeometry/SCT_BarrelModuleSideDesign.h"
#include "SCT_ReadoutGeometry/SCT_ForwardModuleSideDesign.h"
#include "SCT_ReadoutGeometry/StripStereoAnnulusDesign.h"


#include "InDetIdentifier/PixelID.h"
#include "InDetReadoutGeometry/SiDetectorManager.h"

#include "InDetReadoutGeometry/SiDetectorElement.h"
#include "InDetIdentifier/SCT_ID.h"

namespace Acts{
    class Surface;
    class DetectorBuilder;
    class GeometryIdGenerator;
}

struct moduleInfo {
  bool pixel;
  float module_width;
  float module_length;
  float center0;
  float center1;
  float center2;
  int columns;
  int rows;
  float thickness;
  int side;
};

class PixelID;
class SCT_ID;

using detector_t = detray::detector<detray::default_metadata<detray::array<double>>>;

class GeoModelConversion : public AthAlgorithm
{
	public:
		GeoModelConversion(const std::string& n, ISvcLocator* pSvcLocator);
		virtual ~GeoModelConversion() = default;

		virtual StatusCode initialize() override;
		virtual StatusCode execute() override;
		virtual StatusCode finalize() override;

	private:

		void build_geometry_maps();
		void build_detray_detector();
		void load_detray_detector();
		void createAthenaMap();
		void fill_digi_info();

		std::map<Identifier, moduleInfo> m_atlasModuleInfo;
		std::map<std::uint64_t, Identifier> m_DetrayToAtlasMap;

		const PixelID* m_pixelId;
        const SCT_ID*  m_stripId;
		vecmem::host_memory_resource host_mr;
		//using host_detector_type = detray::detector<detray::default_metadata<detray::array<double>>>;
        detector_t m_detray{host_mr};

		std::map<Acts::GeometryIdentifier, Identifier> m_actsToGeoModel;
		std::map<int,std::shared_ptr<Acts::Surface>> m_mapSurfaceToAtlasID;

		ToolHandle<IActsTrackingGeometryTool> m_trackingGeometryTool{this, "TrackingGeometryTool","ActsTrackingGeometryTool"};
        std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
        std::vector<std::shared_ptr<Acts::Surface>> m_surfaces;
		std::vector<std::shared_ptr<Acts::Surface>> m_allsurfaces;

		SG::ReadCondHandleKey<InDetDD::SiDetectorElementCollection> m_pixelDetEleCollKey{this, "PixelDetEleCollKey", "ITkPixelDetectorElementCollection", "Key of SiDetectorElementCollection for Pixel"};
        SG::ReadCondHandleKey<InDetDD::SiDetectorElementCollection> m_stripDetEleCollKey{this, "SCTDetEleCollKey", "ITkStripDetectorElementCollection", "Key of SiDetectorElementCollection for SCT"};

        detector_t m_detector{host_mr};


};

/*
  Copyright (C) 2002-2025 CERN for the benefit of the ATLAS collaboration
*/

#ifndef EFTRACKING_TRACKCONVERSIONTOOL_H
#define EFTRACKING_TRACKCONVERSIONTOOL_H

// Athena includes
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ServiceHandle.h"
#include "StoreGate/ReadHandleKey.h"
#include "GaudiKernel/IChronoStatSvc.h"

#include "ActsGeometryInterfaces/ITrackingGeometryTool.h"
#include "ActsToolInterfaces/IActsToTrkConverterTool.h"
#include "ActsEvent/TrackContainerHandlesHelper.h"
#include "ActsEvent/Seed.h"
#include "ActsEvent/SeedContainer.h"


#include "InDetIdentifier/PixelID.h"
#include "InDetIdentifier/SCT_ID.h"
#include "InDetReadoutGeometry/SiDetectorElement.h"

#include "TrkEventPrimitives/PdgToParticleHypothesis.h"
#include "xAODInDetMeasurement/PixelClusterContainer.h"
#include "xAODInDetMeasurement/StripClusterContainer.h"
#include "xAODInDetMeasurement/SpacePointContainer.h"

#include "traccc/edm/track_state.hpp"
#include "detray/tracks/bound_track_parameters.hpp"

#include "EFTracking/ITrackConversionTool.h"
#include "EFTracking/ITrackingHitInputTool.h"
#include "EFTracking/ITrackingRecoTool.h"
#include "ActsToolInterfaces/ISeedingTool.h"

class TrackConversionTool : public extends<AthAlgTool, ITrackConversionTool> {

  public:
        TrackConversionTool(const std::string&, const std::string&, const IInterface*);
        virtual ~TrackConversionTool() { ; }

        StatusCode initialize() override;
        StatusCode finalize()   override;
        StatusCode convertTracks(
            EventContext const & eventContext
          , traccc::track_state_container_types::host const & resolved_tracks
          , std::map<int,int> const & cluster_map
          , unsigned & nb_output_tracks
          ) override;
        StatusCode convertSeeds(
            EventContext const & eventContext
          , traccc::edm::seed_collection::host const & seed_tracks
          , traccc::edm::spacepoint_collection::host const & traccc_spacepoints
          , traccc::measurement_collection_types::host const & traccc_measurements
          , std::map<int,int> const & cluster_map) override;
        void setAthenaDetrayConversionMaps(
          std::unordered_map<uint64_t, Identifier> const & detray_to_athena_map,
          std::unordered_map<Identifier, uint64_t> const & athena_to_detray_map
          ) override;

  private:

        std::optional<Acts::BoundTrackParameters> convertToActsParameters(const traccc::track_state<traccc::default_algebra>& state);
        const Acts::Surface& actsIdToActsSurface(Acts::GeometryIdentifier& actsID);

        PublicToolHandle<ActsTrk::ITrackingGeometryTool> m_trackingGeometryTool{this, "TrackingGeometryTool", "ActsTrackingGeometryTool"};

        std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
        std::map<Acts::GeometryIdentifier, const Acts::Surface*> m_actsSurfaceMap;
        std::map<Identifier, Acts::GeometryIdentifier> m_AtlasToActs;
        std::unordered_map<Identifier, std::uint64_t> const * m_AthenaToDetrayMap;
        std::unordered_map<std::uint64_t, Identifier> const * m_DetrayToAthenaMap;

        // Create output
        SG::WriteHandleKey<ActsTrk::TrackContainer> m_ActsTracccTrackContainerKey{this, "ActsTracccTracks", "ActsTracccTracks", "Output track collection (ActsTrk variant)"};
        SG::WriteHandleKey<ActsTrk::SeedContainer> m_ActsTracccSeedTrackContainerKey{this, "ActsTracccSeedTracks", "ActsTracccSeedTracks", "Output track seed collection (ActsTrk variant)"};

        SG::ReadHandleKey<xAOD::PixelClusterContainer> m_xAODPixelClusterFromTracccClusterKey{this, "xAODPixelClusterFromTracccClusterKey","xAODPixelClustersFromTracccCluster","Traccc cluster->xAOD PixelClusters Container"};
        SG::ReadHandleKey<xAOD::StripClusterContainer> m_xAODStripClusterFromTracccClusterKey{this, "xAODStripClusterFromTracccClusterKey","xAODStripClustersFromTracccCluster","Traccc cluster ->xAOD StripClusters Container"};
        SG::ReadHandleKey<xAOD::SpacePointContainer> m_xAODSpacepointFromTracccClusterKey{this, "xAODSpacepointFromTracccClusterKey","xAODSpacepointFromTracccCluster","Traccc cluster->xAOD Spacepoint Container"};

        SG::ReadHandleKey<xAOD::PixelClusterContainer> m_xAODPixelClusterFromInDetClusterKey{this, "xAODPixelClusterFromInDetClusterKey","xAODPixelClustersFromInDetCluster","InDet cluster->xAOD PixelClusters Container"};
        SG::ReadHandleKey<xAOD::StripClusterContainer> m_xAODStripClusterFromInDetClusterKey{this, "xAODStripClusterFromInDetClusterKey","xAODStripClustersFromInDetCluster","InDet cluster ->xAOD StripClusters Container"};
        SG::ReadHandleKey<xAOD::SpacePointContainer> m_xAODSpacepointFromInDetClusterKey{this, "xAODSpacepointFromInDetClusterKey","xAODSpacepointFromInDetCluster","InDet cluster->xAOD Spacepoint Container"};

        int getClusterIndex(const traccc::measurement& meas);

        const PixelID* m_pixelID;
        const SCT_ID*  m_stripID;

        // acts helper for the output
        ActsTrk::MutableTrackContainerHandlesHelper m_tracksBackendHandlesHelper{this};
        ActsTrk::MutableTrackContainerHandlesHelper m_seedTracksBackendHandlesHelper{this};
        Acts::BoundTrackParameters translateTrackParametersToActsParameters(detray::bound_track_parameters<traccc::default_algebra>& fit_params);
        Trk::PdgToParticleHypothesis m_pdgToParticleHypothesis;

        ServiceHandle<IChronoStatSvc>   m_chrono{this, "ChronoStatService", "ChronoStatSvc"};

        Gaudi::Property<bool> m_doTruth{this, "doTruth", true, "Create truth links between measurements on track and xAOD clusters"};
        Gaudi::Property<bool> m_recoFromHits{this, "recoFromHits", true, "Match measurements on track to right container"};
        Gaudi::Property<std::string> m_filesDir{this, "filesDir", "/eos/project/a/atlas-eftracking/GPU/ITk_data/", "full path to location with detector files"};

        float chi2_min;
        float chi2_max;
        float chi2_sum;
        float ndf_min;
        float ndf_max;
        float ndf_sum;
        unsigned meas_min;
        unsigned meas_max;
        unsigned meas_sum;
        unsigned excluded_ndf;
        unsigned excluded_weird_sp;
        unsigned excluded_no_sp;
};

#endif

/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#ifndef EFTRACKING_TRACKCONVERSIONTOOL_H
#define EFTRACKING_TRACKCONVERSIONTOOL_H

// Athena includes
#include "AthenaBaseComps/AthAlgTool.h"
#include "AthenaBaseComps/AthService.h"
#include "GaudiKernel/ServiceHandle.h"
#include "StoreGate/ReadHandleKey.h"

#include "ActsGeometryInterfaces/IActsTrackingGeometryTool.h"
#include "ActsEventCnv/IActsToTrkConverterTool.h"
#include "ActsEvent/TrackContainerHandlesHelper.h"

#include "TrkEventPrimitives/PdgToParticleHypothesis.h"

#include "traccc/edm/track_state.hpp"
#include "detray/tracks/bound_track_parameters.hpp"

#include "EFTracking/ITrackConversionTool.h"
#include "EFTracking/ITrackingRecoTool.h"
#include "EFTracking/ITrackingHitMapTool.h"


class TrackConversionTool : public extends<AthAlgTool, ITrackConversionTool> {

  public:
        TrackConversionTool(const std::string&, const std::string&, const IInterface*);
        virtual ~TrackConversionTool() { ; }

        virtual StatusCode initialize() override;
        virtual StatusCode finalize()   override;
        virtual StatusCode convertTracks(const EventContext& eventContext, traccc::track_state_container_types::host& resolved_tracks, std::map<int, Identifier>& atlasHumanIDtoIdentifier, std::map<std::uint64_t, int>& detrayToAtlasMap) override;


  private:
        const Acts::Surface& atlasIdToActsSurface(const Identifier &atlasID);
        const Acts::BoundTrackParameters convertToActsParameters(const traccc::track_state<detray::cmath<float> >& state, std::map<int, Identifier>& atlasHumanIDtoIdentifier, std::map<std::uint64_t, int>& detrayToAtlasMap);
        const Acts::Surface& atlasIdToActsSurface(Identifier& atlasID);

        ToolHandle<IActsTrackingGeometryTool> m_trackingGeometryTool{this, "TrackingGeometryTool", "ActsTrackingGeometryTool"};
        ToolHandle<ITrackingHitMapTool> m_hitMapTool {this, "HitMapTool", "TrackingHitMapTool", "Hit Mapping Tool"};
        std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
        std::map<Identifier, const Acts::Surface*> m_actsSurfaceMap;
        
        ToolHandle<ActsTrk::IActsToTrkConverterTool> m_converterTool{this, "converterTool", "ActsToTrkConverterTool", ""};

        // Create output
        SG::WriteHandleKey<ActsTrk::TrackContainer> m_trackContainerKey{this, "ActsTracccTracks", "ActsTracccTracks", "Output track collection (ActsTrk variant)"};
        // acts helper for the output
        ActsTrk::MutableTrackContainerHandlesHelper m_tracksBackendHandlesHelper;
        Acts::BoundTrackParameters translateTrackParametersToActsParameters(detray::bound_track_parameters<detray::cmath<float> >& fit_params);
        Trk::PdgToParticleHypothesis m_pdgToParticleHypothesis;
         
};

#endif

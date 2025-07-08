/*
  Copyright (C) 2002-2025 CERN for the benefit of the ATLAS collaboration
*/

#ifndef EFTRACKING_ITRACKCONVERSIONTOOL_H
#define EFTRACKING_ITRACKCONVERSIONTOOL_H

#include "Identifier/Identifier.h"
#include "traccc/edm/track_state.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/edm/seed_collection.hpp"
#include <GaudiKernel/EventContext.h>
#include <GaudiKernel/IAlgTool.h>

class ITrackConversionTool : virtual public ::IAlgTool
{
  public:
    DeclareInterfaceID(ITrackConversionTool, 1, 0);
    virtual ~ITrackConversionTool() = default;
    virtual StatusCode convertTracks(EventContext const & eventContext
      , traccc::track_state_container_types::host const & resolved_tracks
      , std::map<int,int> const & clusterMap
      , unsigned & nb_output_tracks) = 0;
    virtual StatusCode convertSeeds(EventContext const & eventContext
      , traccc::edm::seed_collection::host const & seed_tracks
      , traccc::edm::spacepoint_collection::host const & traccc_spacepoints
      , traccc::measurement_collection_types::host const & traccc_measurements
      , std::map<int,int> const & clusterMap) = 0;
    virtual void setAthenaDetrayConversionMaps(
      std::unordered_map<uint64_t, Identifier> const & detray_to_athena_map,
      std::unordered_map<Identifier, uint64_t> const & athena_to_detray_map
      ) = 0;
};

#endif

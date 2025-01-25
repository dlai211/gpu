/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/


#ifndef EFTRACKING_ITRACKCONVERSIONTOOL_H
#define EFTRACKING_ITRACKCONVERSIONTOOL_H

#include "GeoPrimitives/GeoPrimitives.h"

class ITrackConversionTool : virtual public ::IAlgTool {

  public:
        DeclareInterfaceID(ITrackConversionTool, 1, 0);
        virtual ~ITrackConversionTool() = default;
        virtual StatusCode convertTracks(const EventContext& eventContext, traccc::track_state_container_types::host& resolved_tracks, std::map<int, Identifier>& atlasHumanIDtoIdentifier, std::map<std::uint64_t, int>& detrayToAtlasMap) = 0;

};

#endif 

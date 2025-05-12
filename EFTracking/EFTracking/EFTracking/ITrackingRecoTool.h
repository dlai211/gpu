/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/


#ifndef EFTRACKING_ITRACKINGRECOTOOL_H
#define EFTRACKING_ITRACKINGRECOTOOL_H

#include "GeoPrimitives/GeoPrimitives.h"
#include "EFTracking/ITrackingHitInputTool.h"
#include "EFTracking/ITrackingHitMapTool.h"

class ITrackingRecoTool : virtual public ::IAlgTool {

  public:
        DeclareInterfaceID(ITrackingRecoTool, 1, 0);
        virtual ~ITrackingRecoTool() = default;

        // virtual traccc::track_state_container_types::host doRecoFromClusters(std::vector<clusterInfo>& detray_clusters) = 0;
        // virtual traccc::track_state_container_types::host doRecoFromHits(std::vector<hitInfo>& detray_hits) = 0;
        virtual traccc::track_state_container_types::host doRecoFromClusters(std::vector<hitInfo>& detray_hits, std::vector<clusterInfo>& detray_clusters) = 0;


};

#endif 

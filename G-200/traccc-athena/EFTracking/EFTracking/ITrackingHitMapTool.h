
/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/


#ifndef EFTRACKING_ITRACKINGHITMAPTOOL_H
#define EFTRACKING_ITRACKINGHITMAPTOOL_H

#include "GeoPrimitives/GeoPrimitives.h"
#include "EFTracking/ITrackingHitInputTool.h"
#include "traccc/geometry/silicon_detector_description.hpp"

struct moduleInfo {
  bool pixel;
  float module_width;
  float module_length;
  int columns;
  int rows;
  float thickness;
};

class ITrackingHitMapTool : virtual public ::IAlgTool {

  public:
        DeclareInterfaceID(ITrackingHitMapTool, 1, 0);
        virtual ~ITrackingHitMapTool() = default;

        virtual std::vector<hitInfo> mapHits(std::vector<hitInfo>& hits, std::map<int, std::uint64_t>& AtlasToDetrayMap) = 0;
        virtual std::vector<clusterInfo> mapClusters(std::vector<clusterInfo>& clusters, std::map<int, std::uint64_t>& AtlasToDetrayMap) = 0;
        virtual std::tuple<std::map<std::uint64_t, int>,std::map<int, std::uint64_t>,traccc::silicon_detector_description::host,std::map<int, Identifier>> createMaps() = 0;

};

#endif 

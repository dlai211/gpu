
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
  int side;
  int shift;
  float center0;
  float center1;
  float center2;
};

class ITrackingHitMapTool : virtual public ::IAlgTool {

  public:
        DeclareInterfaceID(ITrackingHitMapTool, 1, 0);
        virtual ~ITrackingHitMapTool() = default;

        virtual std::vector<hitInfo> mapHits(std::vector<hitInfo>& hits, std::map<Identifier, std::uint64_t>& AtlasToDetrayMap) = 0;
        // virtual std::vector<hitInfo> mapHits(std::vector<hitInfo>& hits, std::map<Identifier, std::vector<uint64_t>>& AtlasToDetrayMap) = 0;
        virtual std::vector<clusterInfo> mapClusters(std::vector<clusterInfo>& clusters, std::map<Identifier, std::uint64_t>& AtlasToDetrayMap) = 0;
        // virtual std::vector<clusterInfo> mapClusters(std::vector<clusterInfo>& clusters, std::map<Identifier, std::vector<uint64_t>>& AtlasToDetrayMap) = 0;
        virtual std::tuple<std::map<std::uint64_t, Identifier>,std::map<Identifier, std::uint64_t>,traccc::silicon_detector_description::host,std::map<Identifier, Identifier>> createMaps() = 0;

};

#endif 

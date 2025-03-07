/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/


#ifndef EFTRACKING_ITRACKINGHITINPUTTOOL_H
#define EFTRACKING_ITRACKINGHITINPUTTOOL_H

#include "GeoPrimitives/GeoPrimitives.h"
#include "Identifier/Identifier.h"
#include "Acts/Geometry/GeometryIdentifier.hpp"

struct clusterInfo {
  Identifier atlas_id;
  Acts::GeometryIdentifier acts_id;
  std::uint64_t detray_id;
  unsigned int local_key;
  Amg::Vector3D globalPosition;
  Amg::Vector2D localPosition;
  std::vector<float> localCov;
  bool pixel;
};

struct hitInfo {
  Identifier atlas_id;
  Acts::GeometryIdentifier acts_id;
  std::uint64_t detray_id;
  Amg::Vector3D globalPosition;
  int channel0;
  int channel1;
  float value;
  int timestamp;
};

class ITrackingHitInputTool : virtual public ::IAlgTool {

  public:
        DeclareInterfaceID(ITrackingHitInputTool, 1, 0);
        virtual ~ITrackingHitInputTool() = default;

        virtual std::vector<hitInfo> readHits( const EventContext& eventContext) = 0;
        virtual std::vector<clusterInfo> readClusters( const EventContext& eventContext) = 0;

};

#endif 

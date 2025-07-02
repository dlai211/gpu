/*
  Copyright (C) 2002-2025 CERN for the benefit of the ATLAS collaboration
*/

#ifndef EFTRACKING_ITRACKINGHITINPUTTOOL_H
#define EFTRACKING_ITRACKINGHITINPUTTOOL_H

#include "GeoPrimitives/GeoPrimitives.h"
#include "Identifier/Identifier.h"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include <GaudiKernel/IAlgTool.h>
#include <GaudiKernel/EventContext.h>
#include "traccc/edm/silicon_cluster_collection.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/edm/measurement.hpp"

struct clusterInfo {
  Identifier atlas_id;
  Acts::GeometryIdentifier acts_id;
  std::uint64_t detray_id;
  int measurement_id;
  int index;
  unsigned int local_key;
  Amg::Vector3D globalPosition;
  Amg::Vector2D localPosition;
  std::vector<float> localCov;
  bool pixel;

  clusterInfo()
    : atlas_id()
    , acts_id()
    , detray_id()
    , measurement_id(0)
    , index(0)
    , local_key(0)
    , globalPosition()
    , localPosition()
    , localCov()
    , pixel(true)
  {}
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


  hitInfo()
    : atlas_id(0)
    , acts_id(0)
    , detray_id(0)
    , globalPosition({0., 0., 0.})
    , channel0(0)
    , channel1(0)
    , value(0.)
    , timestamp(0)
    {}
};

class ITrackingHitInputTool : virtual public ::IAlgTool
{
  public:
    DeclareInterfaceID(ITrackingHitInputTool, 1, 0);
    virtual ~ITrackingHitInputTool() = default;

    virtual std::pair<std::vector<clusterInfo>, std::map<int,int>> readClusters(
      const EventContext& eventContext
      ) = 0;
    virtual std::map<int,int> convertClusters(
      const EventContext& eventContext,
      traccc::edm::silicon_cluster_collection::host& traccc_clusters,
      traccc::measurement_collection_types::host& traccc_measurements,
      traccc::edm::silicon_cell_collection::host& traccc_cells
      ) = 0;
    virtual void setAthenaDetrayConversionMaps(
      std::unordered_map<uint64_t, Identifier> const & detray_to_athena_map,
      std::unordered_map<Identifier, uint64_t> const & athena_to_detray_map
      ) = 0;

};

#endif

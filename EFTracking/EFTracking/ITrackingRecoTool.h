/*
  Copyright (C) 2002-2025 CERN for the benefit of the ATLAS collaboration
*/

#ifndef EFTRACKING_ITRACKINGRECOTOOL_H
#define EFTRACKING_ITRACKINGRECOTOOL_H

#include "EFTracking/ITrackingHitInputTool.h"

class ITrackingRecoTool : virtual public ::IAlgTool
{
  public:
    DeclareInterfaceID(ITrackingRecoTool, 1, 0);
    virtual ~ITrackingRecoTool() = default;

    virtual StatusCode doRecoFromClusters(const EventContext& evtcontext,std::vector<clusterInfo>& detray_clusters, std::map<int,int>& clusterMap) = 0;
    virtual StatusCode doRecoFromHits(const EventContext& evtcontext) = 0;
    virtual StatusCode doRecoAndWriteFromHits(const EventContext& evtcontext) = 0;
    virtual StatusCode makeTracccStandaloneData(int m_n_event,std::vector<clusterInfo>& detray_clusters, std::string const& outDir) = 0;
    virtual StatusCode makeTracccStandaloneData(int m_n_event,const EventContext& evtcontext, std::string const& outDir) = 0;
    virtual void setAthenaDetrayConversionMaps(
      std::unordered_map<uint64_t, Identifier> const & detray_to_athena_map,
      std::unordered_map<Identifier, uint64_t> const & athena_to_detray_map
      ) = 0;
};

#endif

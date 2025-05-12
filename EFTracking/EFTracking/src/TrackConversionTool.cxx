
/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#include "StoreGate/DataHandle.h"
#include "TrackConversionTool.h"
#include <fstream>

#include "ActsGeometry/ActsDetectorElement.h"
#include "traccc/edm/track_state.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "traccc/definitions/track_parametrization.hpp"
#include "traccc/edm/track_state.hpp"

TrackConversionTool::TrackConversionTool(const std::string& algname, const std::string& name, const IInterface* ifc) :
  base_class(algname, name, ifc) {}

StatusCode TrackConversionTool::initialize() {

  ATH_MSG_DEBUG("Initializing conversion tool");
  ATH_CHECK(m_trackContainerKey.initialize());
  
  if (!m_trackingGeometryTool.empty()) {
    ATH_CHECK(m_trackingGeometryTool.retrieve());
    m_trackingGeometry = m_trackingGeometryTool->trackingGeometry();

    m_trackingGeometry->visitSurfaces([&](const Acts::Surface *surface) {
      // find acts surface with the same detector element ID
      if (!surface)
        return;
      const auto *actsElement = dynamic_cast<const ActsDetectorElement *>(
          surface->associatedDetectorElement());
      if (!actsElement)  return;
      // Conversion from Acts to ATLAS surface impossible for the TRT so the TRT
      // surfaces are not stored in this map
      bool isTRT = (dynamic_cast<const InDetDD::TRT_BaseElement *>(
                        actsElement->upstreamDetectorElement()) != nullptr);
      if (isTRT)  return;

      auto [it, ok] =  m_actsSurfaceMap.insert({actsElement->identify(), surface});
      if (!ok) {
        ATH_MSG_WARNING("ATLAS ID " << actsElement->identify()
                                    << " has two ACTS surfaces: "
                                    << it->second->geometryId() << " and "
                                    << surface->geometryId());
      }
    });
  }

  ATH_CHECK(m_converterTool.retrieve());
  ATH_CHECK(m_hitMapTool.retrieve());

  ATH_CHECK(m_tracksBackendHandlesHelper.initialize(ActsTrk::prefixFromTrackContainerName(m_trackContainerKey.key())));
  ATH_MSG_DEBUG("Creating track container");

  ATH_MSG_DEBUG("Initializing complete");
  return StatusCode::SUCCESS;

}

const Acts::Surface& TrackConversionTool::atlasIdToActsSurface(Identifier& atlasID){

  auto it = m_actsSurfaceMap.find(atlasID);
  if (it != m_actsSurfaceMap.end()) {
    return *it->second;
  }
  ATH_MSG_ERROR("No Acts surface corresponding to this ATLAS id:");
  ATH_MSG_ERROR(atlasID);
  throw std::domain_error("No Acts surface corresponding to the ATLAS one");
}


const Acts::BoundTrackParameters TrackConversionTool::convertToActsParameters(const traccc::track_state<detray::cmath<float> >& state, std::map<Identifier, Identifier>& atlasHumanIDtoIdentifier,std::map<std::uint64_t, Identifier>& detrayToAtlasMap){

  using namespace Acts::UnitLiterals;
  std::shared_ptr<const Acts::Surface> actsSurface;
  Acts::BoundVector params; // = Acts::BoundVector::Random();

  // get the associated surface
  const traccc::measurement& measurement = state.get_measurement();
  const std::uint_least64_t detray_id = measurement.surface_link.value();

  if(detrayToAtlasMap.find(detray_id) != detrayToAtlasMap.end()){
    Identifier atlas_humanID = detrayToAtlasMap[detray_id];
    Identifier atlas_id = atlasHumanIDtoIdentifier[atlas_humanID];
    // get the associated surface
    try {
      actsSurface = atlasIdToActsSurface(atlas_id).getSharedPtr();
    } catch (const std::exception &e) {
      ATH_MSG_ERROR("Could not find ACTS detector surface for this atlas ID:");
      ATH_MSG_ERROR(atlas_id);
      throw;  // Nothing we can do, so just pass exception on...
    }
  }
  // no associated surface create a perigee one
  else {
    ATH_MSG_INFO(
        "traccc track parameters to Acts parameters:: No associated surface found." 
        << " Creating a free surface. Trk parameters:");
    
    /*switch (static_cast<int>(detray_id.shape_id())){
      case Trk::SurfaceType::Plane:
        actsSurface = Acts::Surface::makeShared<const Acts::PlaneSurface>( atlasParameter.associatedSurface().transform());
        break;
      case Trk::SurfaceType::Perigee:
        actsSurface = Acts::Surface::makeShared<const Acts::PerigeeSurface>(  atlasParameter.associatedSurface().transform());
        break;
      // TODO - implement the missing types?
      default:
        ATH_MSG_WARNING("No surface type found for this Trk::Surface. Creating a perigee surface.");
        actsSurface = Acts::Surface::makeShared<const Acts::PerigeeSurface>(atlasParameter.associatedSurface().center());
    }*/
  }
  

  // Construct track parameters
  ATH_MSG_VERBOSE("Constructing track parameters for this state");
  auto atlasParam = state.smoothed(); // for d0,z0 use bound_local
  params << atlasParam.bound_local()[0], atlasParam.bound_local()[1], atlasParam.phi(), atlasParam.theta(), atlasParam.qop(), atlasParam.time();
  //auto position = measurement.local(); shoudl there be a transform here??
  /*
  if (actsSurface->bounds().type() == Acts::SurfaceBounds::BoundsType::eAnnulus) {
    // Annulus surfaces are constructed differently in Acts/Trk so we need to
    // convert local coordinates
    auto result =  actsSurface->globalToLocal(gctx, position, atlasParameter.momentum());
    if (result.ok()) {
      params << atlasParam.l0, atlasParam.l1, atlasParam.phi0, atlasParam.theta, atlasParam.qop, atlasParam.t;
    } else {
      ATH_MSG_WARNING("Unable to convert annulus surface - globalToLocal failed");
    }
  } else {
    params << atlasParam[Trk::locX], atlasParam[Trk::locY],
    atlasParam[Trk::phi0], atlasParam[Trk::theta],
    atlasParameter.charge() / (atlasParameter.momentum().mag() * 1_MeV), 0.;
  }
  */

  Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Identity();
  const auto D = state.get_measurement().meas_dim;
  //const matrix_type<D, D> covariance = state.measurement_covariance<D>;
  //cov.topLeftCorner(5, 5) = *state.measurement_covariance;
  // Convert the covariance matrix from MeV
  // FIXME: This needs to handle the annulus case as well - currently the cov
  // is wrong for annulus surfaces
  for (int i = 0; i < cov.rows(); i++) {
    cov(i, 4) = cov(i, 4) / 1_MeV;
  }
  for (int i = 0; i < cov.cols(); i++) {
    cov(4, i) = cov(4, i) / 1_MeV;
  }
  
  Trk::ParticleHypothesis hypothesis = Trk::pion;
  int charge = -1; //(atlasParam.qop() > 0) ? 1 : 0;
  // convert hypotheses
  ATH_MSG_VERBOSE("Converting hypothesis");
  float mass = Trk::ParticleMasses::mass[hypothesis] * Acts::UnitConstants::MeV;
  Acts::PdgParticle absPdg = Acts::makeAbsolutePdgParticle( static_cast<Acts::PdgParticle>(
        m_pdgToParticleHypothesis.convert(hypothesis, charge)));
  Acts::ParticleHypothesis actsHypothesis{
    absPdg, mass, Acts::AnyCharge{std::abs(static_cast<float>(charge))}};

  return Acts::BoundTrackParameters(actsSurface, params, cov, actsHypothesis);
}

const Acts::Surface& TrackConversionTool::atlasIdToActsSurface(const Identifier &atlasID) {

  auto it = m_actsSurfaceMap.find(atlasID);
  if (it != m_actsSurfaceMap.end()) {
    return *it->second;
  }
  ATH_MSG_ERROR("No Acts surface corresponding to this ATLAS surface:");
  ATH_MSG_ERROR(atlasID);
  throw std::domain_error("No Acts surface corresponding to the ATLAS one");
}


StatusCode TrackConversionTool::convertTracks(const EventContext& eventContext, traccc::track_state_container_types::host& resolved_tracks, std::map<Identifier, Identifier>& atlasHumanIDtoIdentifier, std::map<std::uint64_t, Identifier>& detrayToAtlasMap){

  ActsTrk::MutableTrackContainer trackContainer;
  auto& trackStateContainer = trackContainer.trackStateContainer();
  SG::WriteHandle<ActsTrk::TrackContainer> trackContainerHandle (m_trackContainerKey, eventContext);

  Acts::GeometryContext tgContext = m_trackingGeometryTool->getGeometryContext(eventContext).context();
  
  std::size_t n_resolved_tracks = resolved_tracks.size();
  for (std::size_t i = 0; i < n_resolved_tracks; i++) {

    auto const& [fit_res, states] = resolved_tracks.at(i);
    ATH_MSG_VERBOSE("Track " << i << " has " << states.size() << " measurements");
    if (states.size() < 1) {
          ATH_MSG_WARNING("No space points associated with this track");
          continue;
    }

    //const auto& fit_param = fit_res.fit_params;
    
    // Create the MutableTrack and add the parameters
    auto actsTrack = trackContainer.getTrack(trackContainer.addTrack());
    
    actsTrack.chi2() = fit_res.chi2;
    actsTrack.nDoF() = fit_res.ndf;

    bool first_state = true;
    //int measurementsCount = 0;
    ATH_MSG_VERBOSE("Loop over track states.");
    for (auto const& st : states) {

      Acts::TrackStatePropMask mask = Acts::TrackStatePropMask::None;
      mask |= Acts::TrackStatePropMask::Smoothed;
      auto index = Acts::MultiTrajectoryTraits::kInvalid;

      if (!first_state) {
        index = actsTrack.tipIndex();
      }
      auto actsTSOS = trackStateContainer.getTrackState(trackStateContainer.addTrackState(mask, index));
      ATH_MSG_VERBOSE("TipIndex: " << actsTrack.tipIndex()
                                   << " TSOS index within trajectory: "
                                   << actsTSOS.index());
      actsTrack.tipIndex() = actsTSOS.index();
      
      const Acts::BoundTrackParameters parameters =  convertToActsParameters(st, atlasHumanIDtoIdentifier, detrayToAtlasMap);
      ATH_MSG_VERBOSE("Track parameters: " << parameters.parameters());
      // Sanity check on positions - TODO
      //if (!actsTrackParameterPositionCheck(parameters, *(tsos->trackParameters()), gctx)){
      //  failedIds.push_back(tsos->trackParameters()->associatedSurface().associatedDetectorElementIdentifier());
      //}

      if (first_state) {
        // This is the first track state, so we need to set the track parameters
        ATH_MSG_VERBOSE("First state of track.");
        actsTrack.parameters() = parameters.parameters();
        actsTrack.covariance() = *parameters.covariance();
        actsTrack.setReferenceSurface(parameters.referenceSurface().getSharedPtr());
        first_state = false;
      } else {
        
        try{
          ATH_MSG_VERBOSE("Next state of track.");
          actsTSOS.setReferenceSurface(parameters.referenceSurface().getSharedPtr());
          ATH_MSG_VERBOSE("Got reference surface.");
          // Since we're converting final tracks, let's assume they're smoothed
          actsTSOS.smoothed() = parameters.parameters();
          ATH_MSG_VERBOSE("Got smoothed parameters.");
          actsTSOS.smoothedCovariance() = *parameters.covariance();
          // Not yet implemented in MultiTrajectory.icc
          // actsTSOS.typeFlags() |= Acts::TrackStateFlag::ParameterFlag;
          if (!(actsTSOS.hasSmoothed() && actsTSOS.hasReferenceSurface())) {
            ATH_MSG_VERBOSE("TrackState does not have smoothed state ["
                            << actsTSOS.hasSmoothed()
                            << "] or reference surface ["
                            << actsTSOS.hasReferenceSurface() << "].");
          } else {
            ATH_MSG_VERBOSE("TrackState has smoothed state and reference surface.");
          }
        }catch (const std::exception& e){
          ATH_MSG_ERROR("Unable to convert TrackParameter with exception ["<<e.what()<<"]. Will be missing from ACTS track."
                        <<(actsTSOS.smoothed()));
        }
      }

    }
    ATH_MSG_VERBOSE("Done with states of this track.");
   
  }
  ATH_MSG_INFO("Wrote out "<< trackContainer.size() << " tracks.");
  
  std::unique_ptr<ActsTrk::TrackContainer> constTracksContainer = m_tracksBackendHandlesHelper.moveToConst(std::move(trackContainer),m_trackingGeometryTool->getGeometryContext(eventContext).context(), eventContext);  
  ATH_CHECK(trackContainerHandle.record(std::move(constTracksContainer)));

  return StatusCode::SUCCESS;
}


StatusCode TrackConversionTool::finalize() {
  return StatusCode::SUCCESS;
}



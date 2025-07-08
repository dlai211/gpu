
/*
  Copyright (C) 2002-2025 CERN for the benefit of the ATLAS collaboration
*/

#include "TrackConversionTool.h"
#include <Acts/EventData/TrackParameters.hpp>
#include <limits>
#include <stdexcept>

#include "ActsEvent/TrackContainer.h"
#include "ActsGeometry/ActsDetectorElement.h"
#include "ActsGeometry/ATLASSourceLink.h"
#include "AthenaBaseComps/AthMsgStreamMacros.h"
#include "Identifier/Identifier.h"
#include "traccc/edm/track_state.hpp"
#include "detray/tracks/bound_track_parameters.hpp"
#include "traccc/edm/track_state.hpp"
#include "AthenaDetrayConversion.h"

TrackConversionTool::TrackConversionTool(const std::string& algname, const std::string& name, const IInterface* ifc)
  : base_class(algname, name, ifc)
  , m_AthenaToDetrayMap(nullptr)
  , m_DetrayToAthenaMap(nullptr)
{}

StatusCode TrackConversionTool::initialize() {

  ATH_MSG_DEBUG("Initializing conversion tool");
  ATH_CHECK(m_ActsTracccTrackContainerKey.initialize());
  ATH_CHECK(m_ActsTracccSeedTrackContainerKey.initialize());

  ATH_CHECK(detStore()->retrieve(m_pixelID, "PixelID"));
  ATH_CHECK(detStore()->retrieve(m_stripID, "SCT_ID"));

  ATH_CHECK(m_xAODPixelClusterFromTracccClusterKey.initialize());
  ATH_CHECK(m_xAODStripClusterFromTracccClusterKey.initialize());
  ATH_CHECK(m_xAODSpacepointFromTracccClusterKey.initialize());

  ATH_CHECK(m_xAODPixelClusterFromInDetClusterKey.initialize());
  ATH_CHECK(m_xAODStripClusterFromInDetClusterKey.initialize());
  ATH_CHECK(m_xAODSpacepointFromInDetClusterKey.initialize());

  ATH_CHECK(m_chrono.retrieve());
  ATH_MSG_DEBUG("Retrieved timing tool.");

  if (!m_trackingGeometryTool.empty()) {
    ATH_CHECK(m_trackingGeometryTool.retrieve());
    m_trackingGeometry = m_trackingGeometryTool->trackingGeometry();

    m_trackingGeometry->visitSurfaces([&](const Acts::Surface *surface) {
      // find acts surface with the same detector element ID
      if (!surface){
        ATH_MSG_ERROR("Null Acts::Surface encountered");
        return;
      }
      const auto *actsElement = dynamic_cast<const ActsDetectorElement *>(
          surface->associatedDetectorElement());
      if (!actsElement){
        ATH_MSG_ERROR("Null Acts::DetectorElement encountered");
        return;
      }

      Acts::GeometryIdentifier actsID = surface->geometryId();
      auto [it, ok] =  m_actsSurfaceMap.insert({actsID, surface});
      if (!ok) {
        ATH_MSG_WARNING("ACTS ID " << actsID
                                    << " has two ACTS surfaces: "
                                    << it->second->geometryId() << " and "
                                    << surface->geometryId());
      }

      const auto *geoElement = actsElement->upstreamDetectorElement();
      const auto *detElem = dynamic_cast<const InDetDD::SiDetectorElement*>(geoElement);

      if (!geoElement){

        ATH_MSG_DEBUG("We have this ACTS detector element, but it does not have a corresponding Athena detector element: " << actsID);
        return;
      }
      if(!detElem){

        ATH_MSG_DEBUG("We have this Athena detector element, but it is not a SiDetector element: " << actsID);
        return;
      }

      int nPix = 0;
      int nStrip = 0;
      Identifier athenaID;

      if(detElem->isPixel()){
        nPix++;
        athenaID = detElem->identify();

      } else if(detElem->isSCT()){
        nStrip++;
        athenaID = m_stripID->module_id(detElem->identify());
        const IdentifierHash Strip_ModuleHash = m_stripID->wafer_hash(athenaID);

        // Extract the correct Strip_ModuleID
        int side = m_stripID->side(detElem->identify());
        athenaID = m_stripID->wafer_id(Strip_ModuleHash + side);
      }

      m_AtlasToActs[athenaID] = actsID;
    });
  }

  ATH_CHECK(m_tracksBackendHandlesHelper.initialize(ActsTrk::prefixFromTrackContainerName(m_ActsTracccTrackContainerKey.key())));
  ATH_CHECK(m_seedTracksBackendHandlesHelper.initialize(ActsTrk::prefixFromTrackContainerName(m_ActsTracccSeedTrackContainerKey.key())));

  ATH_MSG_DEBUG("Initializing complete");
  return StatusCode::SUCCESS;
}

void TrackConversionTool::setAthenaDetrayConversionMaps(
  std::unordered_map<uint64_t, Identifier> const & detray_to_athena_map,
  std::unordered_map<Identifier, uint64_t> const & athena_to_detray_map
  )
{
  m_AthenaToDetrayMap = &athena_to_detray_map;
  m_DetrayToAthenaMap = &detray_to_athena_map;
}

StatusCode TrackConversionTool::finalize() {
  return StatusCode::SUCCESS;
}

const Acts::Surface& TrackConversionTool::actsIdToActsSurface(Acts::GeometryIdentifier& actsID)
{
  auto it = m_actsSurfaceMap.find(actsID);
  if (it != m_actsSurfaceMap.end()) {
    return *it->second;
  }
  ATH_MSG_ERROR("No Acts surface corresponding to this ACTS id:");
  ATH_MSG_ERROR(actsID);
  throw std::domain_error("No Acts surface corresponding to the ACTS ID");
}


std::optional<Acts::BoundTrackParameters> TrackConversionTool::convertToActsParameters(const traccc::track_state<traccc::default_algebra>& state)
{
  // get the associated surface
  const traccc::measurement& measurement = state.get_measurement();
  const std::uint_least64_t detray_id = measurement.surface_link.value();

  using namespace Acts::UnitLiterals;
  std::shared_ptr<const Acts::Surface> actsSurface;
  Acts::BoundVector params; // = Acts::BoundVector::Random();

  Identifier const atlas_ID = findAthenaID(detray_id, m_DetrayToAthenaMap);
  Acts::GeometryIdentifier acts_ID = m_AtlasToActs[atlas_ID];
  // get the associated surface
  try {
    actsSurface = actsIdToActsSurface(acts_ID).getSharedPtr();
  } catch (const std::exception &e) {
    ATH_MSG_ERROR("Could not find ACTS detector surface for this ID: " << acts_ID);

    throw;  // Nothing we can do, so just pass exception on...
  }

  // Construct track parameters
  ATH_MSG_VERBOSE("Constructing track parameters for this state");
  auto atlasParam = state.smoothed(); // for d0,z0 use bound_local
  //ATH_MSG_INFO("\td0= " << atlasParam.bound_local()[0] << " z0=" << atlasParam.bound_local()[1] << " phi=" << atlasParam.phi() << " theta=" << atlasParam.theta() << " qoverp=" << atlasParam.qop());

  if (std::isnan(atlasParam.bound_local()[0])
          || std::isnan(atlasParam.bound_local()[1])
          || atlasParam.bound_local()[0] == 0.
          || atlasParam.bound_local()[1] == 0.
          )
  {
    return std::optional<Acts::BoundTrackParameters>();
  }else{
    params << atlasParam.bound_local()[0], atlasParam.bound_local()[1], atlasParam.phi(), atlasParam.theta(), atlasParam.qop(), atlasParam.time();
  }

  //auto position = measurement.local(); shoudl there be a transform here??

  // if (actsSurface->bounds().type() == Acts::SurfaceBounds::BoundsType::eAnnulus) {
  //   // Annulus surfaces are constructed differently in Acts/Trk so we need to
  //   // convert local coordinates
  //   auto result =  actsSurface->globalToLocal(gctx, position, atlasParameter.momentum());
  //   if (result.ok()) {
  //     params << atlasParam.l0, atlasParam.l1, atlasParam.phi0, atlasParam.theta, atlasParam.qop, atlasParam.t;
  //   } else {
  //     ATH_MSG_WARNING("Unable to convert annulus surface - globalToLocal failed");
  //   }
  // } else {
  //   params << atlasParam[Trk::locX], atlasParam[Trk::locY],
  //   atlasParam[Trk::phi0], atlasParam[Trk::theta],
  //   atlasParameter.charge() / (atlasParameter.momentum().mag() * 1_MeV), 0.;
  // }

  constexpr double GeVToMeV = 1000;
  Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Identity();
  cov *= (GeVToMeV*GeVToMeV);
  // const auto D = measurement.meas_dim;
  //const matrix_type<D, D> covariance = state.measurement_covariance<D>;
  //cov.topLeftCorner(5, 5) = *state.measurement_covariance;
  // Convert the covariance matrix from MeV
  // FIXME: This needs to handle the annulus case as well - currently the cov
  // is wrong for annulus surfaces
  // for (int i = 0; i < cov.rows(); i++) {
  //   cov(i, 4) = cov(i, 4) / 1_MeV;
  // }
  // for (int i = 0; i < cov.cols(); i++) {
  //   cov(4, i) = cov(4, i) / 1_MeV;
  // }


  int charge = (atlasParam.qop() > 0) ? 1 : -1;
  // convert hypotheses
  ATH_MSG_VERBOSE("Converting hypothesis");
  Trk::ParticleHypothesis hypothesis = Trk::pion;
  float mass = Trk::ParticleMasses::mass[hypothesis] * Acts::UnitConstants::MeV;
  Acts::PdgParticle absPdg = Acts::makeAbsolutePdgParticle(Acts::ePionPlus);
  Acts::ParticleHypothesis actsHypothesis{
    absPdg, mass, Acts::AnyCharge{std::abs(static_cast<float>(charge))}};

  return Acts::BoundTrackParameters(actsSurface, params, cov, actsHypothesis);
}

StatusCode TrackConversionTool::convertTracks(
      EventContext const & eventContext
    , traccc::track_state_container_types::host const & resolved_tracks
    , std::map<int,int> const & cluster_map
      , unsigned & nb_output_tracks)
{
  nb_output_tracks = 0;

  ActsTrk::MutableTrackContainer trackContainer;
  SG::WriteHandle<ActsTrk::TrackContainer> trackContainerHandle(m_ActsTracccTrackContainerKey, eventContext);

  Acts::GeometryContext tgContext = m_trackingGeometryTool->getGeometryContext(eventContext).context();

  if(m_doTruth){
    ATH_MSG_INFO("Will map truth");
    if(m_recoFromHits){

      ATH_MSG_INFO("Mapping to truth, retrieving cluster container keys: " << m_xAODPixelClusterFromTracccClusterKey.key() << "," << m_xAODStripClusterFromTracccClusterKey.key());
    }else{

      ATH_MSG_INFO("Mapping to truth, retrieving cluster container keys: " << m_xAODPixelClusterFromInDetClusterKey.key() << "," << m_xAODStripClusterFromInDetClusterKey.key());
    }
  }

  int pout = 0;
  int sout = 0;

  // Debug stats
  chi2_sum = 0.;
  chi2_min = std::numeric_limits<float>::max();
  chi2_max = std::numeric_limits<float>::min();;
  ndf_sum = 0.;
  ndf_min = std::numeric_limits<float>::max();
  ndf_max = std::numeric_limits<float>::min();;
  meas_sum = 0.;
  meas_min = std::numeric_limits<unsigned>::max();
  meas_max = std::numeric_limits<unsigned>::min();;

  excluded_ndf = 0;
  excluded_weird_sp = 0;
  excluded_no_sp = 0;

  for (std::size_t i = 0; i < resolved_tracks.size(); i++) {
    auto const& [fit_res, states] = resolved_tracks.at(i);
    // ATH_MSG_INFO("Track " << i << " has " << states.size() << " measurements");
    if (states.size() < 1) {
      // ATH_MSG_WARNING("excluding track " << i << " for no spacepoints");
      excluded_no_sp += 1;
      continue;
    }

    // In Acts ndf, aka nDoF, is unsigned int. This makes sure the number is
    // safe to cast; exclude the track otherwise.
    if (fit_res.trk_quality.ndf > static_cast<float>(std::numeric_limits<unsigned int>::max())
        || fit_res.trk_quality.ndf < static_cast<float>(std::numeric_limits<unsigned int>::min())) {
      // ATH_MSG_WARNING("excluding track " << i << " for out-of-bound ndf: " << fit_res.trk_quality.ndf);
      excluded_ndf += 1;
      continue;
    }

    // Create the MutableTrack and add the parameters
    auto actsTrack = trackContainer.makeTrack();
    bool valid_track = true;

    actsTrack.chi2() = fit_res.trk_quality.chi2;
    actsTrack.nDoF() = fit_res.trk_quality.ndf;

    Acts::TrackStatePropMask const mask = Acts::TrackStatePropMask::Smoothed;

    bool first_state = true;
    ATH_MSG_VERBOSE("Loop over track states.");
    for (auto const& st : states) {
      auto actsTSOS = actsTrack.appendTrackState(mask);

      std::optional<Acts::BoundTrackParameters> params_opt = convertToActsParameters(st);
      if (!params_opt.has_value()) {
        // if the measurement was "buggy", the whole track is dismissed.
        valid_track = false;
        break;
      }
      Acts::BoundTrackParameters const & parameters = params_opt.value();
      ATH_MSG_VERBOSE("Track parameters: " << parameters.parameters());

      // Sanity check on positions - TODO
      //if (!actsTrackParameterPositionCheck(parameters, *(tsos->trackParameters()), gctx)){
      //  failedIds.push_back(tsos->trackParameters()->associatedSurface().associatedDetectorElementIdentifier());
      //}

      if(m_doTruth){

        SG::ReadHandle<xAOD::PixelClusterContainer> inputPixelClusterContainer = m_recoFromHits ? SG::makeHandle( m_xAODPixelClusterFromTracccClusterKey, eventContext ) : SG::makeHandle( m_xAODPixelClusterFromInDetClusterKey, eventContext );
        ATH_CHECK( inputPixelClusterContainer.isValid() );
        const xAOD::PixelClusterContainer *inputPixelClusters = inputPixelClusterContainer.cptr();

        SG::ReadHandle<xAOD::StripClusterContainer> inputStripClusterContainer = m_recoFromHits ? SG::makeHandle( m_xAODStripClusterFromTracccClusterKey, eventContext ) : SG::makeHandle( m_xAODStripClusterFromInDetClusterKey, eventContext );
        ATH_CHECK( inputStripClusterContainer.isValid() );
        const xAOD::StripClusterContainer *inputStripClusters = inputStripClusterContainer.cptr();

        m_chrono->chronoStart("traccc linking to truth");
        // match measurement to the cluster container
        if(cluster_map.empty()){
          ATH_MSG_FATAL("Retrieving Traccc measurement to xAOD cluster map failed.");
        }
        const traccc::measurement& meas = st.get_measurement();
        ActsTrk::ATLASUncalibSourceLink el;
        int cl_index = -1;
        if (auto it = cluster_map.find(meas.measurement_id); it != cluster_map.end()) {
          cl_index = it->second;
        } else {
          std::runtime_error("cannot find measurement in cluster map");
        }

        if(meas.meas_dim == 2u){

          //for debugging purposes, print positions of the measurement and the matched xAOD cluster
          if(pout < 10){
            const xAOD::PixelCluster* theCluster = inputPixelClusters->at(cl_index);
            auto localPos = theCluster->localPosition<2>();
            ATH_MSG_DEBUG("Traccc pixel cluster at " << meas.local[0] << "," << meas.local[1] << " matched to Athena cluster at " << localPos.x() << "," << localPos.y());
            pout++;
          }
          el = ActsTrk::makeATLASUncalibSourceLink(inputPixelClusters,cl_index,eventContext);

        }else{

          //for debugging purposes, print positions of the measurement and the matched xAOD cluster
          if(sout < 10){
            const xAOD::StripCluster* theCluster = inputStripClusters->at(cl_index);
            auto localPos = theCluster->localPosition<1>();
            ATH_MSG_DEBUG("Traccc strip cluster at " << meas.local[0] << "," << meas.local[1] << " matched to Athena cluster at " << localPos.x() );
            sout++;
          }
          el = ActsTrk::makeATLASUncalibSourceLink(inputStripClusters,cl_index,eventContext);

        }

        m_chrono->chronoStop("traccc linking to truth");
        actsTSOS.setUncalibratedSourceLink(Acts::SourceLink(el));
      }

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
          ATH_MSG_ERROR("Unable to convert TrackParameter with exception ["<<e.what()<<"]. Will be missing from ACTS track." <<(actsTSOS.smoothed()));
        }
      }

    }

    ATH_MSG_VERBOSE("Done with states of this track.");
    if (!valid_track) {
      // ATH_MSG_WARNING("excluding track " << i << " for weird spacepoints");
      excluded_weird_sp += 1;
      trackContainer.removeTrack(actsTrack.index());
    }

    // Debug stats
    chi2_sum += fit_res.trk_quality.chi2;
    if (fit_res.trk_quality.chi2 < chi2_min) chi2_min = fit_res.trk_quality.chi2;
    if (fit_res.trk_quality.chi2 > chi2_max) chi2_max = fit_res.trk_quality.chi2;
    ndf_sum += fit_res.trk_quality.ndf;
    if (fit_res.trk_quality.ndf < ndf_min) ndf_min = fit_res.trk_quality.ndf;
    if (fit_res.trk_quality.ndf > ndf_max) ndf_max = fit_res.trk_quality.ndf;
    meas_sum += states.size();
    if (states.size() < meas_min) meas_min = states.size();
    if (states.size() > meas_max) meas_max = states.size();

  }

  // ATH_MSG_INFO("Wrote out "<< trackContainer.size() << " tracks"
  //   << " out of " << resolved_tracks.size()
  //   << ", excluded:"
  //   << " no sp: " << excluded_no_sp
  //   << ", weird sp: " << excluded_weird_sp
  //   << ", ndf: " << excluded_ndf
  //   );
  nb_output_tracks = trackContainer.size();

  // Debug stats
  // ATH_MSG_INFO("chi2 min: " << chi2_min << ", max: " << chi2_max << ", avg: " << chi2_sum / static_cast<float>(trackContainer.size()));
  // ATH_MSG_INFO("ndf min: " << ndf_min << ", max: " << ndf_max << ", avg: " << ndf_sum / static_cast<float>(trackContainer.size()));
  // ATH_MSG_INFO("meas min: " << meas_min << ", max: " << meas_max << ", avg: " << static_cast<float>(meas_sum) / static_cast<float>(trackContainer.size()));

  std::unique_ptr<ActsTrk::TrackContainer> constTracksContainer =
      m_tracksBackendHandlesHelper.moveToConst(
        std::move(trackContainer),
        m_trackingGeometryTool->getGeometryContext(eventContext).context(),
        eventContext);
  ATH_CHECK(trackContainerHandle.record(std::move(constTracksContainer)));

  return StatusCode::SUCCESS;
}

StatusCode TrackConversionTool::convertSeeds(
    EventContext const & eventContext
  , traccc::edm::seed_collection::host const & seed_tracks
  , traccc::edm::spacepoint_collection::host const & traccc_spacepoints
  , traccc::measurement_collection_types::host const & traccc_measurements
  , std::map<int,int> const & cluster_map)
{

  if(cluster_map.empty()){
    ATH_MSG_FATAL("You have provided an empty cluster map, will not be able to do truth matching.");
  }
  if(m_recoFromHits){
    ATH_MSG_INFO("Mapping to truth, retrieving cluster container keys: " << m_xAODPixelClusterFromTracccClusterKey.key() << "," << m_xAODStripClusterFromTracccClusterKey.key());
  }else{
    ATH_MSG_INFO("Mapping to truth, retrieving cluster container keys: " << m_xAODPixelClusterFromInDetClusterKey.key() << "," << m_xAODStripClusterFromInDetClusterKey.key());
  }

  SG::ReadHandle<xAOD::SpacePointContainer> inputSpacepointContainer = m_recoFromHits ? SG::makeHandle( m_xAODSpacepointFromTracccClusterKey, eventContext ) : SG::makeHandle( m_xAODSpacepointFromInDetClusterKey, eventContext );
  ATH_CHECK( inputSpacepointContainer.isValid() );

  SG::WriteHandle<ActsTrk::SeedContainer> seedHandle = SG::makeHandle( m_ActsTracccSeedTrackContainerKey, eventContext );
  ATH_MSG_INFO( "Created Traccc ACTS Seed container with key " << m_ActsTracccSeedTrackContainerKey );
  ATH_CHECK( seedHandle.record( std::make_unique< ActsTrk::SeedContainer >() ) );
  ActsTrk::SeedContainer *seedContainer = seedHandle.ptr();

  // conversion of seeds to xAOD ACTS seeds
  seedContainer->reserve(seed_tracks.size());
  for (size_t st = 0; st < seed_tracks.size(); st++) {
      ATH_MSG_DEBUG("Seed number " << st);
      const auto& seed = seed_tracks.at(st);
      std::vector<unsigned int> sp_traccc_index{seed.bottom_index(),seed.middle_index(),seed.top_index()};
      std::vector<unsigned int> sp_container_index;
      for (int sp = 0; sp < 3; sp++) {

          const auto& sp_traccc = traccc_spacepoints.at(sp_traccc_index[sp]);
          ATH_MSG_DEBUG("Traccc sp: " << sp_traccc.x() << "," << sp_traccc.y() << "," << sp_traccc.z());
          unsigned int meas_traccc_index = sp_traccc.measurement_index_1();
          const traccc::measurement& meas = traccc_measurements.at(meas_traccc_index);
          int cl_index = -1;
          if (auto it = cluster_map.find(meas.measurement_id); it != cluster_map.end()) {
            cl_index = it->second;
          } else {
            std::runtime_error("cannot find measurement in cluster map");
          }

          const xAOD::SpacePoint* xaod_sp = inputSpacepointContainer->at(cl_index);
          ATH_MSG_DEBUG("Matched to: " << xaod_sp->x() << "," << xaod_sp->y() << "," << xaod_sp->z());

          sp_container_index.push_back(cl_index);

      }


      std::unique_ptr< ActsTrk::Seed > actsTracccSeed =
      std::make_unique< ActsTrk::Seed >(*(inputSpacepointContainer->at(sp_container_index[0])),
                *(inputSpacepointContainer->at(sp_container_index[1])),
                *(inputSpacepointContainer->at(sp_container_index[2])));

      seedContainer->push_back(std::move(actsTracccSeed));
  }

  ATH_MSG_INFO( "Wrote Traccc ACTS Seed container with key " << m_ActsTracccSeedTrackContainerKey << " and size: " << seedContainer->size());

  return StatusCode::SUCCESS;
}
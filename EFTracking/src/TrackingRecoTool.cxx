
/*
  Copyright (C) 2002-2023 CERN for the benefit of the ATLAS collaboration
*/

#include "StoreGate/DataHandle.h"
#include "TrackingRecoTool.h"
#include "traccc/definitions/math.hpp"
#include "traccc/fitting/kalman_filter/kalman_fitter.hpp"
#include <format>
#include <cmath>

#include <iostream>
#include <sstream>
#include <fstream>

template <typename scalar_t>
using unit = detray::unit<scalar_t>;

using transform3 = detray::dtransform3D<traccc::default_algebra>;
using point3 = detray::dpoint3D<traccc::default_algebra>;
using point2 = detray::dpoint2D<traccc::default_algebra>;

using host_detector_type = detray::detector<>;
using device_detector_type = detray::detector<detray::default_metadata,
                             detray::device_container_types>;

using b_field_t = covfie::field<detray::bfield::const_bknd_t>;
using rk_stepper_type = detray::rk_stepper<b_field_t::view_t,
                        typename host_detector_type::algebra_type,
                        detray::constrained_step<>>;
using host_navigator_type = detray::navigator<const host_detector_type>;
using host_fitter_type = traccc::kalman_fitter<rk_stepper_type, host_navigator_type>;
using device_navigator_type = detray::navigator<const device_detector_type>;
using device_fitter_type = traccc::kalman_fitter<rk_stepper_type, device_navigator_type>;

typename host_detector_type::geometry_context context{};
typename traccc::host::combinatorial_kalman_filter_algorithm::config_type m_finding_cfg;
typename traccc::host::kalman_fitting_algorithm::config_type m_fitting_cfg;
typename traccc::clustering_config m_cluster_cfg;

TrackingRecoTool::TrackingRecoTool(const std::string& algname, const std::string& name, const IInterface* ifc) :
  base_class(algname, name, ifc) {}

StatusCode TrackingRecoTool::initialize() {

  ATH_MSG_DEBUG("Initializing reco tool");
  
  // Read the detector.
  detray::io::detector_reader_config reader_cfg{};
//   reader_cfg.add_file("/eos/project/a/atlas-eftracking/GPU/ITk_data/ATLAS-P2-RUN4-03-00-00/ITk_DetectorBuilder_geometry.json");
//   reader_cfg.add_file("/eos/project/a/atlas-eftracking/GPU/ITk_data/ATLAS-P2-RUN4-03-00-00/ITk_DetectorBuilder_surface_grids.json");

  reader_cfg.add_file("/eos/user/j/jlai/itk_data/new/ITk_DetectorBuilder_geometry.json");
  reader_cfg.add_file("/eos/user/j/jlai/itk_data/new/ITk_DetectorBuilder_surface_grids.json");

  auto [host_det, _] = detray::io::read_detector<host_detector_type>(host_mr, reader_cfg);
  m_detector = std::move(host_det);

//   traccc::io::read_detector_description(m_dd,"/eos/project/a/atlas-eftracking/GPU/ITk_data/ATLAS-P2-RUN4-03-00-00/ITk_DetectorBuilder_geometry.json","/eos/project/a/atlas-eftracking/GPU/ITk_data/ATLAS-P2-RUN4-03-00-00/ITk_digitization_config.json", traccc::data_format::json);
  traccc::io::read_detector_description(m_dd,"/eos/user/j/jlai/itk_data/new/ITk_DetectorBuilder_geometry.json","/eos/user/j/jlai/itk_data/new/ITk_digitization_config.json", traccc::data_format::json);
  
  m_bFieldInZ = 1.99724;

  m_n_event = 0;

  ATH_MSG_DEBUG("Setting up configs");
  m_seedfinder.zMin = -3000;
  m_seedfinder.zMax = 3000;
  m_seedfinder.rMax = 320;
  m_seedfinder.collisionRegionMin = -200;
  m_seedfinder.collisionRegionMax = 200;
  m_seedfinder.minPt = 900.;
  m_seedfinder.cotThetaMax = 27.2899;
  m_seedfinder.deltaRMin = 20;
  m_seedfinder.deltaRMax = 280;
  m_seedfinder.impactMax = 2.;
  m_seedfinder.sigmaScattering = 2.0;
  m_seedfinder.maxPtScattering = 10e6;
  m_seedfinder.maxSeedsPerSpM = 5;

  m_finding_cfg.min_track_candidates_per_track = 3;
  m_finding_cfg.max_track_candidates_per_track = 10;
  m_finding_cfg.min_step_length_for_next_surface = 0.5;
  m_finding_cfg.max_step_counts_for_next_surface = 100;
  m_finding_cfg.chi2_max = 10;
  m_finding_cfg.max_num_branches_per_seed = 5;
  m_finding_cfg.max_num_skipping_per_cand = 2;
  m_finding_cfg.propagation.stepping.min_stepsize = 1e-4f;
  m_finding_cfg.propagation.stepping.rk_error_tol = 1e-4f;
  m_finding_cfg.propagation.stepping.step_constraint = 1.f;
  m_finding_cfg.propagation.stepping.path_limit = 5.f;
  m_finding_cfg.propagation.stepping.max_rk_updates = 10000u;
  m_finding_cfg.propagation.stepping.use_mean_loss = true;
  m_finding_cfg.propagation.stepping.use_eloss_gradient = false;
  m_finding_cfg.propagation.stepping.use_field_gradient = false;
  m_finding_cfg.propagation.stepping.do_covariance_transport = true;

  m_fitting_cfg.propagation.navigation.min_mask_tolerance = 1e-5f;
  m_fitting_cfg.propagation.navigation.max_mask_tolerance = 1.f;
  m_fitting_cfg.propagation.navigation.overstep_tolerance = -100.f;
  m_fitting_cfg.propagation.navigation.search_window[0] = 0u;
  m_fitting_cfg.propagation.navigation.search_window[1] = 0u;

 
  ATH_MSG_DEBUG("Initializing complete");
  return StatusCode::SUCCESS;

}

// Comparator used for sorting cells. This sorting is one of the assumptions
// made in the clusterization algorithm
struct cell_order{
    bool operator()(hitInfo& lhs, hitInfo& rhs){
        if (lhs.detray_id != rhs.detray_id){
            return(lhs.detray_id < rhs.detray_id);
        }
        else if (lhs.channel1 != rhs.channel1) {
            return (lhs.channel1 < rhs.channel1);
        } else {
            return (lhs.channel0 < rhs.channel0);
        }
    }
};

/// Comparison / ordering operator for measurements
struct measurement_sort_comp{
    bool operator()(clusterInfo& lhs, clusterInfo& rhs){

        if (lhs.detray_id != rhs.detray_id) {
            return lhs.detray_id < rhs.detray_id;
        } else if (lhs.localPosition[0] != rhs.localPosition[0]) {
            return lhs.localPosition[0] < rhs.localPosition[0];
        } else if (lhs.localPosition[1] != rhs.localPosition[1]) {
            return lhs.localPosition[1] < rhs.localPosition[1];
        } 
        return false;
    }
};

traccc::track_state_container_types::host TrackingRecoTool::recoClustersToTracks(traccc::spacepoint_collection_types::host& spacepoints_per_event, traccc::measurement_collection_types::host& measurements_per_event){

    uint64_t n_modules = 0;
    uint64_t n_measurements = 0;
    uint64_t n_measurements_cuda = 0;
    uint64_t n_spacepoints = 0;
    uint64_t n_spacepoints_cuda = 0;
    uint64_t n_seeds = 0;
    uint64_t n_seeds_cuda = 0;
    uint64_t n_params = 0;
    uint64_t n_params_cuda = 0;
    uint64_t n_found_tracks = 0;
    uint64_t n_found_tracks_cuda = 0;
    uint64_t n_fitted_tracks = 0;
    uint64_t n_fitted_tracks_cuda = 0;
    uint64_t n_ambiguity_free_tracks = 0;

    const traccc::vector3 B{0, 0, 2 * detray::unit<traccc::scalar>::T};
    auto field = detray::bfield::create_const_field(B);    

    // Read the surface transforms
    auto surface_transforms = traccc::io::alt_read_geometry(m_detector);
    // Copy detector desc to the device.
    device_detector = detray::get_buffer(detray::get_data(m_detector), device_mr, copy);
    stream.synchronize();
    auto det_view = detray::get_data(device_detector);
 
    // Seeding and track finding algorithms
    traccc::seeding_algorithm sa(m_seedfinder,
                                {m_seedfinder},
                                 m_seedfilter, host_mr);

    traccc::cuda::seeding_algorithm sa_cuda{m_seedfinder,
                                           {m_seedfinder},
                                            m_seedfilter,
                                            mr,
                                            async_copy,
                                            stream};


    traccc::track_params_estimation tp(host_mr);
    traccc::cuda::track_params_estimation tp_cuda{mr, async_copy, stream};
    
    // Finding algorithm object
    traccc::host::combinatorial_kalman_filter_algorithm host_finding(m_finding_cfg);
    traccc::cuda::finding_algorithm<rk_stepper_type, device_navigator_type> device_finding(m_finding_cfg, mr, async_copy, stream);

    traccc::host::kalman_fitting_algorithm host_fitting(m_fitting_cfg, host_mr);
    traccc::cuda::fitting_algorithm<device_fitter_type> device_fitting(m_fitting_cfg, mr, async_copy, stream);

    traccc::greedy_ambiguity_resolution_algorithm host_ambiguity_resolution{};

    // Instantiate host containers/collections
    traccc::seeding_algorithm::output_type seeds;
    traccc::track_params_estimation::output_type params;
    traccc::track_candidate_container_types::host track_candidates;
    traccc::track_state_container_types::host track_states;
    traccc::track_state_container_types::host resolved_track_states_ar;

    traccc::seed_collection_types::buffer seeds_cuda_buffer(0, *(mr.host));
    traccc::bound_track_parameters_collection_types::buffer params_cuda_buffer(0, *mr.host);

    traccc::track_candidate_container_types::buffer
            track_candidates_cuda_buffer{{{}, *(mr.host)},
                                         {{}, *(mr.host), mr.host}};

    traccc::track_state_container_types::buffer track_states_cuda_buffer{ {{}, *(mr.host)}, {{}, *(mr.host), mr.host}};
    traccc::device::container_d2h_copy_alg<traccc::track_candidate_container_types>
        copy_track_candidates(mr, copy);
    traccc::device::container_d2h_copy_alg<traccc::track_state_container_types>
        copy_track_states(mr, copy);
   

    ///////////////////////////////
    //  GPU track reconstrction  //
    ///////////////////////////////
    //----------------------------
    //     Seeding algorithm
    //----------------------------

    ATH_MSG_INFO("Start of track seeding");
    traccc::spacepoint_collection_types::buffer spacepoints_cuda_buffer(spacepoints_per_event.size(), mr.main);
    async_copy(vecmem::get_data(spacepoints_per_event),spacepoints_cuda_buffer);

    traccc::measurement_collection_types::buffer measurements_cuda_buffer(measurements_per_event.size(), mr.main);
    async_copy(vecmem::get_data(measurements_per_event),
                       measurements_cuda_buffer);

    seeds_cuda_buffer = sa_cuda(spacepoints_cuda_buffer);
    stream.synchronize();

    traccc::seed_collection_types::host seeds_cuda;
    async_copy(seeds_cuda_buffer, seeds_cuda)->wait();    

    //----------------------------
    //    Track params estimation
    //----------------------------

    ATH_MSG_INFO("Start of track param estimation");
    params_cuda_buffer = tp_cuda(spacepoints_cuda_buffer, seeds_cuda_buffer, {0.f, 0.f, m_bFieldInZ});
    stream.synchronize();
    traccc::bound_track_parameters_collection_types::host params_cuda;
    async_copy(params_cuda_buffer, params_cuda)->wait();
    
    //----------------------------
    // Track finding and fitting
    //----------------------------

    ATH_MSG_INFO("Start of track finding");
    track_candidates_cuda_buffer = device_finding(det_view, field, measurements_cuda_buffer,
                                   params_cuda_buffer);
    stream.synchronize();
    auto track_candidates_cuda = copy_track_candidates(track_candidates_cuda_buffer);
    
    ATH_MSG_INFO("Start of track fitting");
    track_states_cuda_buffer = device_fitting(det_view, field, track_candidates_cuda_buffer);
    stream.synchronize();
    auto track_states_cuda = copy_track_states(track_states_cuda_buffer);

    //n_modules += sp_out.modules.size();
    n_measurements += measurements_per_event.size();
    n_spacepoints += spacepoints_per_event.size();
    
    n_seeds_cuda += seeds_cuda.size();
    n_params_cuda += params_cuda.size();
    n_found_tracks_cuda += track_candidates_cuda.size();
    n_fitted_tracks_cuda += track_states_cuda.size();

    // CPU
    if(m_doRecoOnCPU){
        seeds = sa(spacepoints_per_event);
        params = tp(spacepoints_per_event, seeds, {0.f, 0.f, m_bFieldInZ});
        track_candidates = host_finding(m_detector, field, vecmem::get_data(measurements_per_event), vecmem::get_data(params));
        track_states = host_fitting(m_detector, field, traccc::get_data(track_candidates));
        ATH_MSG_INFO("Start of ambiguity resolving");
        resolved_track_states_ar = host_ambiguity_resolution(track_states);

        n_seeds += seeds.size();
        n_params += params.size();
        n_found_tracks += track_candidates.size();
        n_fitted_tracks += track_states.size();
        n_ambiguity_free_tracks += resolved_track_states_ar.size();
        
    }

    //----------------------------
    //      Statistics
    //----------------------------
    
    ATH_MSG_INFO("==> Statistics of track reco... " );
    ATH_MSG_INFO("- created  " << n_measurements << " measurements     " );
    ATH_MSG_INFO("- created  " << n_spacepoints << " spacepoints     " );
    ATH_MSG_INFO("- created  " << n_seeds << " seeds" );
    ATH_MSG_INFO("- created (cuda) " << n_seeds_cuda << " seeds" );
    ATH_MSG_INFO("- created  " << n_params << " initial parameters" );
    ATH_MSG_INFO("- created (cuda) " << n_params_cuda << " initial parameters" );
    ATH_MSG_INFO("- found    " << n_found_tracks << " tracks" );
    ATH_MSG_INFO("- found (cuda)   " << n_found_tracks_cuda << " tracks"  );
    ATH_MSG_INFO("- fitted   " << n_fitted_tracks << " tracks" );
    ATH_MSG_INFO("- fitted (cuda)  " << n_fitted_tracks_cuda << " tracks" );
    ATH_MSG_INFO("- resolved   " << n_ambiguity_free_tracks << " tracks" );

    m_modules += n_modules;
    m_measurements += n_measurements;
    m_measurements_cuda += n_measurements_cuda;
    m_spacepoints += n_spacepoints;
    m_spacepoints_cuda += n_spacepoints_cuda;
    m_seeds += n_seeds;
    m_seeds_cuda += n_seeds_cuda;
    m_params += n_params;
    m_params_cuda += n_params_cuda;
    m_found_tracks += n_found_tracks;
    m_found_tracks_cuda += n_found_tracks_cuda;
    m_fitted_tracks += n_fitted_tracks;
    m_fitted_tracks_cuda += n_fitted_tracks_cuda;
    m_ambiguity_free_tracks += n_ambiguity_free_tracks;

    return track_states_cuda;
}


traccc::track_state_container_types::host TrackingRecoTool::doRecoFromClusters(std::vector<hitInfo>& detray_hits, std::vector<clusterInfo>& detray_clusters) {

    // ATH_MSG_INFO("Start of hits -> track reconstruction");
    // ATH_MSG_INFO("Number of hits to reco: " << detray_hits.size());

    ATH_MSG_INFO("Start of hits saving");
    ATH_MSG_INFO("Number of hits: " << detray_hits.size());

    std::sort(detray_hits.begin(), detray_hits.end(), cell_order());
    
    // Read the surface transforms
    auto surface_transforms = traccc::io::alt_read_geometry(m_detector);
    // Copy detector desc to the device.
    device_detector = detray::get_buffer(detray::get_data(m_detector), device_mr, copy);
    stream.synchronize();
    auto det_view = detray::get_data(device_detector);
    
    // Instantiate host containers/collections
    traccc::edm::silicon_cell_collection::host cells_per_event{host_mr};
    traccc::silicon_detector_description::host host_det_descr{host_mr};

    read_cells(cells_per_event, detray_hits);
    m_cells += cells_per_event.size();

    // auto context = Gaudi::Hive::currentContext();
    // std::size_t event_n = context.eventID().event_number();
    std::cout << "writing cell data" << std::endl;
    traccc::io::write(m_n_event,"/eos/user/j/jlai/g200/gpu/G-200/traccc-athena/run",traccc::data_format::csv, vecmem::get_data(cells_per_event),vecmem::get_data(m_dd),false);
        
    ATH_MSG_INFO("Start of cluster saving");
    ATH_MSG_INFO("Number of clusters: " << detray_clusters.size());

    // Print the first 10 clusters for debugging
    int count = 0;
    for (const auto& cluster : detray_clusters) {  
        if (count <= 10) {
            ATH_MSG_INFO("Cluster " << count + 1
                        << " - atlas_id: " << cluster.atlas_id
                        << ", geometry_id: " << cluster.detray_id
                        << ", globalPosition: (" 
                        << cluster.globalPosition.x() << ", "
                        << cluster.globalPosition.y() << ", " 
                        << cluster.globalPosition.z() << ")"
                        << ", localPosition: (" 
                        << cluster.localPosition.x() << ", "
                        << cluster.localPosition.y() << ")"
                        << ", localCov(0,0): " << cluster.localCov[0]
                        << ", localCov(1,1): " << cluster.localCov[1]
                        << ", isPixel: " << cluster.pixel);
            count++;
        } else {
            break;
        }
    }


    // ATH_MSG_INFO("Start of cluster -> track reconstruction");
    // ATH_MSG_INFO("Number of clusters to reco: " << detray_clusters.size());
    std::sort(detray_clusters.begin(), detray_clusters.end(), measurement_sort_comp());

    traccc::spacepoint_collection_types::host spacepoints_per_event{&host_mr};
    traccc::measurement_collection_types::host measurements_per_event{&host_mr};

    read_spacepoints(spacepoints_per_event, detray_clusters, true);
    read_measurements(measurements_per_event, detray_clusters, true);

    make_test_data(m_n_event,detray_clusters);

    // auto context = Gaudi::Hive::currentContext();
    // std::size_t event_n = context.eventID().event_number();
    traccc::io::write(m_n_event,"/eos/user/j/jlai/g200/ITk_data/ITk_hit_data/",traccc::data_format::binary, vecmem::get_data(measurements_per_event));
    traccc::io::write(m_n_event,"/eos/user/j/jlai/g200/ITk_data/ITk_hit_data/",traccc::data_format::binary, vecmem::get_data(spacepoints_per_event));
    // traccc::io::write(m_n_event,"/eos/user/n/nribaric/tracccData/meas_test/athena_clusters/",traccc::data_format::binary, vecmem::get_data(measurements_per_event));
    // traccc::io::write(m_n_event,"/eos/user/n/nribaric/tracccData/meas_test/athena_clusters/",traccc::data_format::binary, vecmem::get_data(spacepoints_per_event));

    ATH_MSG_INFO("Number of clusters to reco: " << measurements_per_event.size());

    return {};
    // traccc::track_state_container_types::host track_states = recoClustersToTracks(spacepoints_per_event, measurements_per_event);
    // return track_states;

}



// traccc::track_state_container_types::host TrackingRecoTool::doRecoFromClusters(std::vector<clusterInfo>& detray_clusters){

    
//     ATH_MSG_INFO("Start of cluster saving");
//     ATH_MSG_INFO("Number of clusters: " << detray_clusters.size());

//     // ATH_MSG_INFO("Start of cluster -> track reconstruction");
//     // ATH_MSG_INFO("Number of clusters to reco: " << detray_clusters.size());
//     std::sort(detray_clusters.begin(), detray_clusters.end(), measurement_sort_comp());

//     traccc::spacepoint_collection_types::host spacepoints_per_event{&host_mr};
//     traccc::measurement_collection_types::host measurements_per_event{&host_mr};

//     read_spacepoints(spacepoints_per_event, detray_clusters, false);
//     read_measurements(measurements_per_event, detray_clusters, false);

//     make_test_data(m_n_event,detray_clusters);

//     // auto context = Gaudi::Hive::currentContext();
//     // std::size_t event_n = context.eventID().event_number();
//     traccc::io::write(m_n_event,"/eos/user/j/jlai/g200/ITk_data/ITk_hit_data/",traccc::data_format::binary, vecmem::get_data(measurements_per_event));
//     traccc::io::write(m_n_event,"/eos/user/j/jlai/g200/ITk_data/ITk_hit_data/",traccc::data_format::binary, vecmem::get_data(spacepoints_per_event));
//     // traccc::io::write(m_n_event,"/eos/user/n/nribaric/tracccData/meas_test/athena_clusters/",traccc::data_format::binary, vecmem::get_data(measurements_per_event));
//     // traccc::io::write(m_n_event,"/eos/user/n/nribaric/tracccData/meas_test/athena_clusters/",traccc::data_format::binary, vecmem::get_data(spacepoints_per_event));

//     ATH_MSG_INFO("Number of clusters to reco: " << measurements_per_event.size());

//     return {};
//     // traccc::track_state_container_types::host track_states = recoClustersToTracks(spacepoints_per_event, measurements_per_event);
//     // return track_states;

// }

void TrackingRecoTool::read_measurements(traccc::measurement_collection_types::host& measurements, std::vector<clusterInfo>& detray_clusters, bool do_strip){

  std::map<traccc::geometry_id, unsigned int> m;
  int n_strip = 0;

  std::pair<uint64_t,detray::geometry::barcode> last_sf;
  std::multimap<uint64_t,detray::geometry::barcode> sf_seen;

  for(std::vector<clusterInfo>::size_type i = 0; i < detray_clusters.size();i++){

    clusterInfo cluster = detray_clusters[i];
    if(do_strip && cluster.pixel){continue;}

    uint64_t geometry_id = cluster.detray_id;
    const point3 gl_pos{cluster.globalPosition[0],cluster.globalPosition[1],cluster.globalPosition[2]};
    const auto& sf = detray::geometry::barcode{geometry_id};
    const detray::tracking_surface<detray::detector<> > surface{m_detector, sf};
    const auto loc_pos = surface.transform(context).point_to_local(gl_pos);
    cluster.localPosition[0] = loc_pos[0];
    cluster.localPosition[1] = loc_pos[1];
    
    // Construct the measurement object.
    traccc::measurement meas;
    std::array<typename transform3::size_type, 2u> indices{0u, 0u};
    meas.meas_dim = 0u;
    for (unsigned int ipar = 0; ipar < 2u; ++ipar) {
        if (((cluster.local_key) & (1 << (ipar + 1))) != 0) {

            switch (ipar) {
                case 0: {
                    meas.local[0] = loc_pos[0]; //cluster.localPosition[0];
                    meas.variance[0] = 0.0025;
                    indices[meas.meas_dim++] = ipar;
                }; break;
                case 1: {
                    meas.local[1] = loc_pos[1]; //cluster.localPosition[1];
                    meas.variance[1] = 0.0025;
                    indices[meas.meas_dim++] = ipar;
                }; break;
            }
	    }
    }

    meas.subs.set_indices(indices);
    meas.surface_link = detray::geometry::barcode{geometry_id};

    // Keeps measurement_id for ambiguity resolution
    meas.measurement_id = i;
    measurements.push_back(meas);
    n_strip++;
  }

}

void TrackingRecoTool::read_spacepoints(traccc::spacepoint_collection_types::host& spacepoints, std::vector<clusterInfo>& detray_clusters, bool do_strip){

  traccc::measurement_collection_types::host measurements;
  read_measurements(measurements, detray_clusters, do_strip);

  std::map<traccc::geometry_id, unsigned int> m;
  int n_strip = 0;
  for(std::vector<clusterInfo>::size_type i = 0; i < detray_clusters.size();i++){
    clusterInfo cluster = detray_clusters[i];
    if(do_strip && cluster.pixel){continue;}
    
    // Construct the global 3D position of the spacepoint.
    const point3 pos{cluster.globalPosition[0], cluster.globalPosition[1], cluster.globalPosition[2]};

    // Construct the local 3D(2D) position of the measurement.
    traccc::measurement meas;
    meas = measurements[i];

    spacepoints.push_back({pos, meas});
    n_strip++;

  }
}

// traccc::track_state_container_types::host TrackingRecoTool::doRecoFromHits(std::vector<hitInfo>& detray_hits){

//     // ATH_MSG_INFO("Start of hits -> track reconstruction");
//     // ATH_MSG_INFO("Number of hits to reco: " << detray_hits.size());

//     ATH_MSG_INFO("Start of hits saving");
//     ATH_MSG_INFO("Number of hits: " << detray_hits.size());

//     std::sort(detray_hits.begin(), detray_hits.end(), cell_order());
    
//     // Read the surface transforms
//     auto surface_transforms = traccc::io::alt_read_geometry(m_detector);
//     // Copy detector desc to the device.
//     device_detector = detray::get_buffer(detray::get_data(m_detector), device_mr, copy);
//     stream.synchronize();
//     auto det_view = detray::get_data(device_detector);
    
//     // Instantiate host containers/collections
//     traccc::edm::silicon_cell_collection::host cells_per_event{host_mr};
//     traccc::silicon_detector_description::host host_det_descr{host_mr};

//     read_cells(cells_per_event, detray_hits);
//     m_cells += cells_per_event.size();

//     // auto context = Gaudi::Hive::currentContext();
//     // std::size_t event_n = context.eventID().event_number();
//     std::cout << "writing cell data" << std::endl;
//     traccc::io::write(m_n_event,"/eos/user/j/jlai/g200/ITk_data/ITk_hit_data/",traccc::data_format::csv, vecmem::get_data(cells_per_event),vecmem::get_data(m_dd),false);

//     return {};
// }


    // ATH_MSG_INFO("Read cells " << cells_per_event.size());
    
    // // Create device copy of input collections
    // traccc::edm::silicon_cell_collection::buffer cells_buffer(cells_per_event.size(), mr.main);
    // async_copy(vecmem::get_data(cells_per_event), cells_buffer);
    // ATH_MSG_DEBUG("created cuda cell buffer");

    // traccc::host::clusterization_algorithm::output_type measurements_per_event;
    // traccc::host::silicon_pixel_spacepoint_formation_algorithm::output_type spacepoints_per_event;

    // // Instantiate cuda containers/collections
    // traccc::measurement_collection_types::buffer measurements_cuda_buffer(0, *mr.host);
    // traccc::spacepoint_collection_types::buffer spacepoints_cuda_buffer(0, *mr.host);
    // ATH_MSG_DEBUG("created cuda meas and sp buffer");

    // traccc::cuda::clusterization_algorithm ca_cuda(mr, async_copy, stream, m_cluster_cfg);
    // traccc::cuda::measurement_sorting_algorithm ms_cuda(async_copy, stream);
    // traccc::cuda::spacepoint_formation_algorithm<traccc::default_detector::device> sf_cuda(mr, async_copy, stream);
    // ATH_MSG_DEBUG("created algos");

    // traccc::silicon_detector_description::data host_dd_data{vecmem::get_data(m_dd)};
    // traccc::silicon_detector_description::buffer device_dd{static_cast<traccc::silicon_detector_description::buffer::size_type>(m_dd.size()),device_mr};
    // copy.setup(device_dd)->wait();
    // copy(host_dd_data, device_dd)->wait();
    // ATH_MSG_DEBUG("copied det descr");


    // //--------------------------------------------
    // //    Clusterization and spacepoint formation
    // //--------------------------------------------
    // // CPU

    // ATH_MSG_INFO("Start of clusterization and spacepoint formation");

    // traccc::host::clusterization_algorithm ca(host_mr);
    // measurements_per_event = ca(vecmem::get_data(cells_per_event), host_dd_data);

    // traccc::host::silicon_pixel_spacepoint_formation_algorithm sf(host_mr);
    // spacepoints_per_event = sf(m_detector, vecmem::get_data(measurements_per_event));

    // // GPU
    // measurements_cuda_buffer = ca_cuda(cells_buffer, device_dd);
    // ms_cuda(measurements_cuda_buffer);
    // stream.synchronize();

    // spacepoints_cuda_buffer = sf_cuda(det_view, measurements_cuda_buffer);
    // stream.synchronize();

    // traccc::measurement_collection_types::host measurements_per_event_cuda;
    // traccc::spacepoint_collection_types::host spacepoints_per_event_cuda;
    // copy(measurements_cuda_buffer, measurements_per_event_cuda)->wait();
    // copy(spacepoints_cuda_buffer, spacepoints_per_event_cuda)->wait();


    // // if(m_doRecoOnCPU){
    // //     traccc::host::clusterization_algorithm ca(host_mr);
    // //     measurements_per_event = ca(vecmem::get_data(cells_per_event), host_dd_data);

    // //     traccc::host::silicon_pixel_spacepoint_formation_algorithm sf(host_mr);
    // //     spacepoints_per_event = sf(m_detector, vecmem::get_data(measurements_per_event));
    // // }

    // ATH_MSG_INFO("measurements: " << measurements_per_event.size());
    // ATH_MSG_INFO("spacepoints: " << spacepoints_per_event.size());

    // // write_traccc_data(m_n_event,spacepoints_per_event);

    // // this is because ATLAS data event numbering is not continuous and does not start with 0
    // // this then upsets traccc :)
    // // m_n_event++;

    // // traccc::io::write(event_n,"/eos/user/n/nribaric/tracccData/meas_test/traccc_clusters/",traccc::data_format::binary, vecmem::get_data(measurements_per_event));
    // // traccc::io::write(event_n,"/eos/user/n/nribaric/tracccData/meas_test/traccc_clusters/",traccc::data_format::binary, vecmem::get_data(spacepoints_per_event));

    // // ATH_MSG_INFO("Adding strip clusters");
    // // read_spacepoints(spacepoints_per_event, detray_clusters, true);
    // // read_measurements(measurements_per_event, detray_clusters, true);

    // // ATH_MSG_INFO("measurements: " << measurements_per_event.size());
    // // ATH_MSG_INFO("spacepoints: " << spacepoints_per_event.size());

    // traccc::track_state_container_types::host track_states = recoClustersToTracks(spacepoints_per_event_cuda, measurements_per_event_cuda);
    // return track_states;
  
// }


void TrackingRecoTool::read_cells(traccc::edm::silicon_cell_collection::host& cells, std::vector<hitInfo>& detray_hits) {

    cells.resize(detray_hits.size()+1);

    // If there is a detector description object, build a map of geometry IDs
    // to indices inside the detector description. We do not use ACTS IDs, so we map to detray IDs
    std::map<uint64_t, unsigned int> geomIdMap;
    
    for (unsigned int i = 0; i < m_dd.acts_geometry_id().size(); ++i) {
        geomIdMap[m_dd.geometry_id()[i].value()] = i;
    }
    
    // Fill the output containers with the ordered cells and modules.
    for(std::vector<hitInfo>::size_type i = 0; i < detray_hits.size(); i++){

        hitInfo hit = detray_hits[i];
        unsigned int ddIndex = 0;
        uint64_t geometry_id = hit.detray_id;
        if (hit.detray_id == 0 or geometry_id == 0) {
            ATH_MSG_INFO("you have not mapped this ID (detray/atlas): " << hit.detray_id << "," << hit.atlas_id);
            continue;
        }


        auto it = geomIdMap.find(geometry_id);
        if (it == geomIdMap.end()) {
            ATH_MSG_INFO("not sure what to do with this ID (detray/atlas): " << geometry_id << "," << hit.atlas_id);
            continue;
        }
        ddIndex = it->second;
        
        // Helper lambda for setting up SoA cells.
        auto assign = [ddIndex](auto soa_cell, hitInfo& hit) {
            soa_cell.channel0() = hit.channel0;
            soa_cell.channel1() = hit.channel1;
            soa_cell.activation() = hit.value;
            soa_cell.time() = hit.timestamp;
            soa_cell.module_index() = ddIndex;
        };

        //const std::size_t cellIndex = cells.size();
        assign(cells.at(i), hit);
        
    } 
    
}

std::vector<hitInfo> TrackingRecoTool::make_fake_data(traccc::silicon_detector_description::host& dd){

    std::ofstream t_cells_file("atlas_hits_diagonal.txt",std::ofstream::out);

    std::vector<uint64_t> v = {256986819764241407,256986822448597247,265994027072050047,275001274376760703};
    ATH_MSG_INFO("Making fake data " << dd.acts_geometry_id().size());
    std::vector<hitInfo> fake_hits;
    for (unsigned int i = 0; i < dd.acts_geometry_id().size(); ++i) {

        if(dd.dimensions()[i] != 2){continue;}
        if(std::find(v.begin(), v.end(), dd.acts_geometry_id()[i]) == v.end()) {continue;}
       
        float x_pitch = dd.pitch_x()[i];
        float y_pitch = dd.pitch_y()[i];
        float x_ref = dd.reference_x()[i];
        float y_ref = dd.reference_y()[i];
        int x_seg = (x_ref/(-0.5))/x_pitch;
        int y_seg = (y_ref/(-0.5))/y_pitch;
        ATH_MSG_INFO(dd.acts_geometry_id()[i] << ": " << x_ref << "," << x_pitch << "," << x_seg);
        ATH_MSG_INFO(y_ref << "," << y_pitch << "," << y_seg);
        

        // make fake data at the module corners
        // hitInfo thishit;
        // thishit.detray_id = dd.acts_geometry_id()[i];
        // thishit.channel0 = 0;
        // thishit.channel1 = 0;
        // thishit.value = 13;
        // thishit.timestamp = 8;
        // fake_hits.push_back(thishit);

        // make fake data along one  side
        // for(int d = 0; d < y_seg; d+=3){

        //     hitInfo thishit;
        //     thishit.detray_id = dd.acts_geometry_id()[i];
        //     thishit.channel0 = 0;
        //     thishit.channel1 = d;
        //     thishit.value = 13;
        //     thishit.timestamp = 8;
        //     fake_hits.push_back(thishit);
        //     t_cells_file << dd.acts_geometry_id()[i] << "," << 0 << "," << d << "\n";

        // }

        // make fake data along the diagonal
        //std::cout << "printig this mod: " << x_seg << "," << y_seg << std::endl;
        int diagonal_length = std::min(x_seg, y_seg);
        for (int d = 0; d < diagonal_length; d += 3) {

            hitInfo thishit;
            // Print the first pair
            thishit.detray_id = dd.acts_geometry_id()[i];
            thishit.channel0 = d;
            thishit.channel1 = d;
            thishit.value = 13;
            thishit.timestamp = 8;
            fake_hits.push_back(thishit);

            t_cells_file << dd.acts_geometry_id()[i] << "," << d << "," << d << "\n";
            //std::cout << "(" << i << ", " << i << ")" << std::endl;

            // Check if there's another element in the diagonal to print the second pair
            if (d + 1 < diagonal_length) {
                //std::cout << "(" << i + 1 << ", " << i + 1 << ")" << std::endl;
                thishit.detray_id = dd.acts_geometry_id()[i];
                thishit.channel0 = d+1;
                thishit.channel1= d+1;
                thishit.value = 13;
                thishit.timestamp = 8;
                fake_hits.push_back(thishit);
                t_cells_file << dd.acts_geometry_id()[i] << "," << d+1 << "," << d+1 << "\n";
            }
            
        }
    }

    t_cells_file.close();
    return fake_hits;
}

StatusCode TrackingRecoTool::finalize() {

    ATH_MSG_INFO("==> Statistics of traccc track reconstruction:" );
    ATH_MSG_INFO("- read from " << m_modules << " modules" );
    ATH_MSG_INFO("- read     " << m_cells << " cells");
    ATH_MSG_INFO("- created  " << m_measurements << " measurements     " );
    ATH_MSG_INFO("- created  " << m_spacepoints << " spacepoints     " );
    ATH_MSG_INFO("- created  " << m_seeds << " seeds" );
    ATH_MSG_INFO("- created (cuda) " << m_seeds_cuda << " seeds" );
    ATH_MSG_INFO("- created  " << m_params << " initial parameters" );
    ATH_MSG_INFO("- created (cuda) " << m_params_cuda << " initial parameters" );
    ATH_MSG_INFO("- found    " << m_found_tracks << " tracks" );
    ATH_MSG_INFO("- found (cuda)   " << m_found_tracks_cuda << " tracks"  );
    ATH_MSG_INFO("- fitted   " << m_fitted_tracks << " tracks" );
    ATH_MSG_INFO("- fitted (cuda)  " << m_fitted_tracks_cuda << " tracks" );
    ATH_MSG_INFO("- resolved   " << m_ambiguity_free_tracks << " tracks" );


    return StatusCode::SUCCESS;
}

void TrackingRecoTool::write_traccc_data(int m_n_event,traccc::spacepoint_collection_types::host& spacepoints){

   
    std::string padded = std::format("{:09}",m_n_event);

    std::ofstream t_ms_file("event"+padded+"-traccc-measurements.txt",std::ofstream::out);
    std::ofstream t_sp_file("event"+padded+"-traccc-spacepoints.txt",std::ofstream::out);

    for(long unsigned int i = 0; i < spacepoints.size(); i++){
        
        auto sp = spacepoints[i];
        auto meas = sp.meas;
        auto local = meas.local;
        auto global = sp.global;
        auto id = meas.surface_link.value();
        
        t_sp_file << id << "," << global[0] << "," << global[1] << "," << global[2] << "\n";
        t_ms_file << id << "," << local[0] << "," << local[1] << "\n";
    }
    
    t_ms_file.close();
    t_sp_file.close();

}    


void TrackingRecoTool::make_test_data(int m_n_event,std::vector<clusterInfo>& detray_clusters){

    
    std::string padded = std::format("{:09}",m_n_event);

    std::ofstream hitfile;
    hitfile.open("event"+padded+"-hits.csv");
    hitfile << "particle_id,geometry_id,tx,ty,tz,tt,tpx,tpy,tpz,te,deltapx,deltapy,deltapz,deltae,index" << "\n";

    std::ofstream measfile;
    measfile << std::fixed << std::setprecision(10);
    measfile.open("event"+padded+"-measurements.csv");
    measfile << "measurement_id,geometry_id,local_key,local0,local1,phi,theta,time,var_local0,var_local1,var_phi,var_theta,var_time" << "\n";

    std::ofstream mapfile;
    mapfile.open("event"+padded+"-measurement-simhit-map.csv");
    mapfile << "measurement_id,hit_id" << "\n";

    for(std::vector<clusterInfo>::size_type i = 0; i < detray_clusters.size(); i++){

        clusterInfo &thiscluster = detray_clusters[i];
        if(thiscluster.detray_id == 0){ATH_MSG_INFO("You have an unmatched cluster! " << i);}

        hitfile << i << "," << thiscluster.detray_id << "," << thiscluster.globalPosition[0] << "," << thiscluster.globalPosition[1] << "," << thiscluster.globalPosition[2];
        hitfile << "," << "0,0,0,0,0,0,0,0,0,0" << "\n";

        measfile << i << "," << thiscluster.detray_id << "," << thiscluster.local_key << "," << thiscluster.localPosition[0] << "," << thiscluster.localPosition[1] << ",";
        // measfile << "0,0,0,0.0025,0.0025,0,0,0" << "\n";
        measfile << "0,0,0," << thiscluster.localCov[0] << "," << thiscluster.localCov[1] << ",0,0,0" << "\n";
        mapfile << i << "," << i << "\n";

    }
    hitfile.close();
    measfile.close();
    mapfile.close();

}




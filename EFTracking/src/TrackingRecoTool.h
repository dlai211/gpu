/*
  Copyright (C) 2002-2025 CERN for the benefit of the ATLAS collaboration
*/

#ifndef EFTRACKING_TRACKINGRECOTOOL_H
#define EFTRACKING_TRACKINGRECOTOOL_H

// Athena includes
#include "AthenaBaseComps/AthAlgTool.h"
#include "AthenaBaseComps/AthMsgStreamMacros.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/IChronoStatSvc.h"
#include "InDetIdentifier/PixelID.h"
#include "InDetIdentifier/SCT_ID.h"
#include "InDetRawData/PixelRDO_Container.h"
#include "InDetRawData/SCT_RDO_Container.h"
#include "InDetReadoutGeometry/SiDetectorDesign.h"
#include "InDetReadoutGeometry/SiDetectorManager.h"
#include "PixelReadoutGeometry/PixelDetectorManager.h"
#include "SCT_ReadoutGeometry/SCT_DetectorManager.h"


// Project include(s).
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/silicon_detector_description.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/device/container_d2h_copy_alg.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/csv/cell.hpp"
#include "traccc/io/write.hpp"

// VecMem include(s).
#include <memory>
#include <traccc/definitions/common.hpp>
#include <traccc/utils/bfield.hpp>
#include <vecmem/memory/binary_page_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// Detray include(s).
#include "detray/core/detector.hpp"
#include "detray/io/frontend/detector_reader.hpp"

// Traccc alg includes
#include "traccc/io/read_bfield.hpp"
#include "traccc/clusterization/clustering_config.hpp"
#include "traccc/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "ActsEvent/MeasurementToTruthParticleAssociation.h"
#include "traccc/utils/seed_generator.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
#include "traccc/clusterization/clusterization_algorithm.hpp"

#include "EFTracking/ITrackingRecoTool.h"
#include "EFTracking/ITrackingHitInputTool.h"
#include "EFTracking/ITrackConversionTool.h"

#include <format>
#include "AthenaDetrayConversion.h"
#include "xAODInDetMeasurement/PixelClusterContainer.h"

template <typename TRAIT> class TrackingRecoTool
    : public extends<AthAlgTool, ITrackingRecoTool>
{
public:
    TrackingRecoTool(const std::string&, const std::string&, const IInterface*);
    virtual ~TrackingRecoTool() { ; }

    StatusCode initialize() override;
    StatusCode finalize()   override;
    StatusCode doRecoFromClusters(const EventContext& evtcontext,std::vector<clusterInfo>& detray_clusters, std::map<int,int>& cluster_map)  override;
    StatusCode doRecoFromHits(const EventContext& evtcontext)  override;
    StatusCode doRecoAndWriteFromHits(const EventContext& evtcontext) override;
    StatusCode makeTracccStandaloneData(int m_n_event,std::vector<clusterInfo>& detray_clusters, std::string const& outDir) override;
    StatusCode makeTracccStandaloneData(int m_n_event,const EventContext& evtcontext, std::string const& outDir) override;
    void setAthenaDetrayConversionMaps(
        std::unordered_map<uint64_t, Identifier> const & detray_to_athena_map,
        std::unordered_map<Identifier, uint64_t> const & athena_to_detray_map
        ) override;

private:

    Gaudi::Property<bool> m_doWriteData{this,"writeData",false,"If you also want to create input for traccc standalone"};

    Gaudi::Property<bool> m_checkSeeds{this,"checkSeeds",true,"If you want to check seeding performance"};
    Gaudi::Property<bool> m_doTruthSeeding{this, "doTruthSeeding", true, "Start with track parameters obtained from truth particles."};
    Gaudi::Property<std::string> m_filesDir{this, "filesDir", "/eos/project/a/atlas-eftracking/GPU/ITk_data/", "full path to location with detector files"};

    SG::ReadHandleKey<ActsTrk::MeasurementToTruthParticleAssociation>  m_pixelClustersToTruth{this, "PixelClustersToTruthParticlesMap", "ITkTracccPixelClustersToTruthParticles", "Association map from pixel measurements to generator particles."};
    SG::ReadHandleKey<xAOD::PixelClusterContainer> m_xAODPixelClusterFromInDetClusterKey{this, "xAODPixelClusterFromInDetClusterKey","xAODPixelClustersFromInDetCluster","InDet cluster->xAOD PixelClusters Container"};

    SG::ReadHandleKey<PixelRDO_Container> m_pixelRDOKey  { this, "PixelRDO", "ITkPixelRDOs" };
    SG::ReadHandleKey<SCT_RDO_Container> m_stripRDOKey  { this, "StripRDO", "ITkStripRDOs" };

    StatusCode clustersToTracks(
          EventContext const & evtcontext
        , traccc::measurement_collection_types::buffer const & measurements_gpu_buffer
        , std::map<int,int> cluster_map
    );

    StatusCode makeTruthSeeds(traccc::bound_track_parameters_collection_types::host& truth_seeds, traccc::measurement_collection_types::host& measurements_host_buffer);
    // Traccc alg configuration(s).
    IntegerProperty m_minTrackCandidates{this, "minTrackCandidates", 3., "min track candidates per track"};
    IntegerProperty m_maxTrackCandidates{this, "maxTrackCandidates", 100., "max track candidates per track"};
    FloatProperty m_minStepLengthToSurface{this, "minStepLengthToSurface", 0.5, "min step length for next surface"};
    IntegerProperty m_maxStepCountsToSurface{this, "maxStepCountsToSurface", 100., "max step counts for next surface"};
    FloatProperty m_maxChi2{this, "maxChi2", 30., "max chi-square allowed for branching"};
    IntegerProperty m_maxBranchesPerSeed{this, "maxBranchesPerSeed", 500, "max number of branches an inital seed can have at each step"};
    IntegerProperty m_maxNumSkip{this, "maxNumSkip", 3, "max number of allowed skipped steps per candidate"};
    IntegerProperty m_targetCellsPerPartition{this, "targetCellsPerPartition", 1024, "adapt to different GPU capabilities"};

    traccc::seedfinder_config m_seedfinder;
    traccc::seedfilter_config m_seedfilter;

    covfie::field<traccc::inhom_bfield_backend_t<traccc::scalar>> m_inhom_field;
    traccc::bfield m_inhom_host_field;
    traccc::bfield m_inhom_gpu_field;

    std::unique_ptr<const Acts::Logger> m_logger;
    std::unordered_map<std::uint64_t, Identifier> const * m_DetrayToAthenaMap;
    std::unordered_map<Identifier, std::uint64_t> const * m_AthenaToDetrayMap;

    const PixelID* m_pixelID;
    const SCT_ID*  m_stripID;
    const InDetDD::PixelDetectorManager* m_pixelManager{nullptr};
    const InDetDD::SCT_DetectorManager* m_stripManager{nullptr};

    typename traccc::host::combinatorial_kalman_filter_algorithm::config_type m_finding_cfg;
    typename traccc::host::kalman_fitting_algorithm::config_type m_fitting_cfg;
    typename traccc::clustering_config m_cluster_cfg;
    typename traccc::host::greedy_ambiguity_resolution_algorithm::config_type m_resolution_config;

    void read_measurements(traccc::measurement_collection_types::host& measurements, std::vector<clusterInfo>& detray_clusters);
    void read_hits(traccc::edm::silicon_cell_collection::host& cells,std::vector<hitInfo>& detray_hits);
    StatusCode read_cells(traccc::edm::silicon_cell_collection::host& cells, const EventContext& evtcontext);

    std::vector<hitInfo> make_fake_data(traccc::silicon_detector_description::host& dd);
    void write_traccc_data(int m_n_event,traccc::edm::spacepoint_collection::host& spacepoints, traccc::measurement_collection_types::host& measurements);
    void make_test_data(int m_n_event,std::vector<clusterInfo>& detray_clusters);

    // Memory resources used by the application.
    std::unique_ptr<vecmem::host_memory_resource> m_host_mr;
    std::unique_ptr<vecmem::copy> m_host_copy;
    std::unique_ptr<typename TRAIT::gpu_host_memory_resource_type> m_gpu_host_mr;
    std::unique_ptr<typename TRAIT::device_memory_resource_type> m_device_mr;
    std::unique_ptr<traccc::memory_resource> m_mr;
    std::unique_ptr<vecmem::binary_page_memory_resource> m_cached_host_mr;
    std::unique_ptr<vecmem::binary_page_memory_resource> m_cached_device_mr;

    // Device-specific types used.
    std::unique_ptr<typename TRAIT::stream_type> m_stream;
    std::unique_ptr<typename TRAIT::copy_type> m_copy;
    std::unique_ptr<typename TRAIT::async_copy_type> m_async_copy;
    std::unique_ptr<traccc::device::container_d2h_copy_alg<traccc::track_state_container_types>> m_copy_track_states;

    traccc::default_detector::host::buffer_type m_device_detector;
    traccc::default_detector::host::view_type m_det_view;
    std::unique_ptr<traccc::default_detector::host> m_detector;
    std::unique_ptr<traccc::silicon_detector_description::host> m_det_descr;
    std::unique_ptr<traccc::silicon_detector_description::buffer> m_device_det_descr;

    ServiceHandle<IChronoStatSvc>   m_chrono{this, "ChronoStatService", "ChronoStatSvc"};
    ToolHandle<ITrackConversionTool> m_cnvTool {this, "TrackConversionTool", "EFTracking/TrackConversionTool", "Track Convrsion Tool"};
    ToolHandle<ITrackingHitInputTool> m_hitInputTool {this,"HitInputTool","EFTracking/TrackingHitInputTool","Hit Input Tool"};

    uint64_t m_modules = 0;
    uint64_t m_cells = 0;
    uint64_t m_measurements_cpu = 0;
    uint64_t m_measurements_gpu = 0;
    uint64_t m_spacepoints_cpu = 0;
    uint64_t m_spacepoints_gpu = 0;
    uint64_t m_seeds_cpu = 0;
    uint64_t m_seeds_gpu = 0;
    uint64_t m_params_cpu = 0;
    uint64_t m_params_gpu = 0;
    uint64_t m_found_tracks_cpu = 0;
    uint64_t m_found_tracks_gpu = 0;
    uint64_t m_fitted_tracks_cpu = 0;
    uint64_t m_fitted_tracks_gpu = 0;
    uint64_t m_ambiguity_free_tracks_cpu = 0;
    uint64_t m_ambiguity_free_tracks_gpu = 0;
    uint64_t m_output_tracks = 0;

    std::unique_ptr<typename TRAIT::clusterization_algorithm_type> m_clusterization_alg_gpu;
    std::unique_ptr<typename TRAIT::measurement_sorting_algorithm_type> m_meas_sort_alg_gpu;
    std::unique_ptr<typename TRAIT::spacepoint_formation_algorithm_type> m_spacepoint_form_alg_gpu;
    std::unique_ptr<typename TRAIT::seeding_algorithm_type> m_seeding_alg_gpu;
    std::unique_ptr<typename TRAIT::track_params_estimation_type> m_track_params_estimation_gpu;
    std::unique_ptr<typename TRAIT::finding_algorithm_type> m_finding_alg_gpu;
    std::unique_ptr<typename TRAIT::fitting_algorithm_type> m_fitting_alg_gpu;

    // some functions that are specialized by template parameter
    void create_copy_object(std::unique_ptr<typename TRAIT::copy_type> &);
    void create_async_copy_object(std::unique_ptr<typename TRAIT::async_copy_type> &);
    void get_stream_for_algo(TRAIT::stream_for_algo_type* &);
    std::string get_init_message();
    traccc::bfield make_gpu_bfield(
          traccc::bfield const & host_bfield
        , TRAIT::stream_for_algo_type & stream);
};

template <typename scalar_t>
using unit = detray::unit<scalar_t>;

using transform3 = detray::dtransform3D<traccc::default_algebra>;
using point3 = detray::dpoint3D<traccc::default_algebra>;
using point2 = detray::dpoint2D<traccc::default_algebra>;

template<typename TRAIT> TrackingRecoTool<TRAIT>::TrackingRecoTool(
    const std::string& algname, const std::string& name, const IInterface* ifc)
    : base_class(algname, name, ifc)
    , m_AthenaToDetrayMap(nullptr)
    , m_DetrayToAthenaMap(nullptr)
{}

template<typename TRAIT> void TrackingRecoTool<TRAIT>::setAthenaDetrayConversionMaps(
    std::unordered_map<uint64_t, Identifier> const & detray_to_athena_map,
    std::unordered_map<Identifier, uint64_t> const & athena_to_detray_map
    )
{
  m_AthenaToDetrayMap = &athena_to_detray_map;
  m_DetrayToAthenaMap = &detray_to_athena_map;
  m_cnvTool->setAthenaDetrayConversionMaps(detray_to_athena_map, athena_to_detray_map);
  m_hitInputTool->setAthenaDetrayConversionMaps(detray_to_athena_map, athena_to_detray_map);
}

template<typename TRAIT> StatusCode TrackingRecoTool<TRAIT>::initialize()
{
    ATH_MSG_INFO(get_init_message());

    m_host_copy = std::make_unique<vecmem::copy>();
    m_host_mr = std::make_unique<vecmem::host_memory_resource>();
    m_gpu_host_mr = std::make_unique<typename TRAIT::gpu_host_memory_resource_type>();
    m_device_mr = std::make_unique<typename TRAIT::device_memory_resource_type>();
    m_cached_host_mr = std::make_unique<vecmem::binary_page_memory_resource>(*m_gpu_host_mr);
    m_cached_device_mr = std::make_unique<vecmem::binary_page_memory_resource>(*m_device_mr);
    m_mr = std::make_unique<traccc::memory_resource>(*m_cached_device_mr, m_cached_host_mr.get());

    m_stream = std::make_unique<typename TRAIT::stream_type>();
    typename TRAIT::stream_for_algo_type* p_stream_obj;
    get_stream_for_algo(p_stream_obj);

    create_copy_object(m_copy);
    create_async_copy_object(m_async_copy);

    ATH_MSG_DEBUG("Retrieving timing tool");
    ATH_CHECK(m_chrono.retrieve());

    ATH_MSG_DEBUG("Initializing Conversion Tool");
    ATH_CHECK(m_cnvTool.retrieve());

    ATH_MSG_DEBUG("Initializing Input Tool for cluster conversion");
    ATH_CHECK(m_hitInputTool.retrieve());

    ATH_CHECK(m_xAODPixelClusterFromInDetClusterKey.initialize());
    ATH_MSG_INFO("Got key for truth seeding");
    ATH_CHECK(m_pixelClustersToTruth.initialize());
    ATH_MSG_INFO("Got truth maps for truth seeding");

    ATH_CHECK(detStore()->retrieve(m_pixelID, "PixelID"));
    ATH_CHECK(m_pixelRDOKey.initialize());

    ATH_CHECK(detStore()->retrieve(m_stripID, "SCT_ID"));
    ATH_CHECK(m_stripRDOKey.initialize());

    ATH_CHECK(detStore()->retrieve(m_pixelManager));
    ATH_CHECK(detStore()->retrieve(m_stripManager));

    ATH_MSG_INFO("Creating magnetic field");
    traccc::io::read_bfield(m_inhom_field, m_filesDir+"ITk_bfield.cvf", traccc::data_format::binary);
    m_inhom_host_field = traccc::bfield{std::move(m_inhom_field)};
    m_inhom_gpu_field = make_gpu_bfield(m_inhom_host_field, *p_stream_obj);

    // Read the detector.
    detray::io::detector_reader_config reader_cfg{};
    reader_cfg.add_file(m_filesDir+"ITk_DetectorBuilder_geometry.json");
    reader_cfg.add_file(m_filesDir+"ITk_detector_material.json");
    reader_cfg.add_file(m_filesDir+"ITk_DetectorBuilder_surface_grids.json");
    #if NDEBUG
    reader_cfg.do_check(false);
    ATH_MSG_INFO("Loading detray detector: checks disabled");
    #else
    ATH_MSG_INFO("Loading detray detector: checks enabled");
    #endif
    auto [host_det, _] = detray::io::read_detector<traccc::default_detector::host>(*m_host_mr, reader_cfg);
    m_detector = std::make_unique<traccc::default_detector::host>(std::move(host_det));

    ATH_MSG_INFO("Loading traccc detector");
    m_det_descr = std::make_unique<traccc::silicon_detector_description::host>(*m_host_mr);
    traccc::io::read_detector_description(*m_det_descr,
        m_filesDir+"ITk_DetectorBuilder_geometry.json",
        m_filesDir+"ITk_digitization_config_with_strips.json",
        traccc::data_format::json);

    ATH_MSG_DEBUG("Setting up configs");
    // m_seedfinder.zMin = -3000.f * unit<float>::mm;
    // m_seedfinder.zMax = 3000.f * unit<float>::mm;
    // m_seedfinder.rMax = 320.f * unit<float>::mm;
    // m_seedfinder.rMin = 33.f * unit<float>::mm;
    // m_seedfinder.collisionRegionMin = -200* unit<float>::mm;
    // m_seedfinder.collisionRegionMax = 200* unit<float>::mm;
    // m_seedfinder.minPt = 500.f * unit<float>::MeV;
    // m_seedfinder.cotThetaMax = 27.2899f;
    // m_seedfinder.deltaRMin = 20* unit<float>::mm;
    // m_seedfinder.deltaRMax = 280* unit<float>::mm;
    // m_seedfinder.impactMax = 2.f * unit<float>::mm;
    // m_seedfinder.sigmaScattering = 100.0f;
    // m_seedfinder.maxPtScattering = 10.f * unit<float>::GeV;
    // m_seedfinder.maxSeedsPerSpM = 1;
    // m_seedfinder.radLengthPerSeed = 0.05f;

    m_finding_cfg.max_num_branches_per_seed = 3;
    m_finding_cfg.max_num_branches_per_surface = 5;
    m_finding_cfg.min_track_candidates_per_track = 7;
    m_finding_cfg.max_track_candidates_per_track = 20;
    m_finding_cfg.min_step_length_for_next_surface = 0.5f * detray::unit<float>::mm;
    m_finding_cfg.max_step_counts_for_next_surface = 100;
    m_finding_cfg.chi2_max = 10.f;
    m_finding_cfg.max_num_skipping_per_cand = 3;

    m_finding_cfg.propagation.stepping.min_stepsize = 1e-4f* unit<float>::mm;
    m_finding_cfg.propagation.stepping.rk_error_tol = 1e-4f* unit<float>::mm;
    m_finding_cfg.propagation.stepping.step_constraint = std::numeric_limits<float>::max();
    m_finding_cfg.propagation.stepping.path_limit = 5.f * unit<float>::m;
    m_finding_cfg.propagation.stepping.max_rk_updates = 10000u;
    m_finding_cfg.propagation.stepping.use_mean_loss = true;
    m_finding_cfg.propagation.stepping.use_eloss_gradient = false;
    m_finding_cfg.propagation.stepping.use_field_gradient = false;
    m_finding_cfg.propagation.stepping.do_covariance_transport = true;
    m_finding_cfg.propagation.navigation.overstep_tolerance = -300.f * unit<float>::um;

    m_fitting_cfg.propagation.navigation.min_mask_tolerance = 1e-5f * unit<float>::mm;
    m_fitting_cfg.propagation.navigation.max_mask_tolerance = 3.f * unit<float>::mm;
    m_fitting_cfg.propagation.navigation.overstep_tolerance = -300.f * unit<float>::um;
    m_fitting_cfg.propagation.navigation.search_window[0] = 0u;
    m_fitting_cfg.propagation.navigation.search_window[1] = 0u;

    m_logger = Acts::getDefaultLogger("EFTrackingGPU", Acts::Logging::INFO);

    m_clusterization_alg_gpu = std::make_unique<typename TRAIT::clusterization_algorithm_type>(
        *m_mr, *m_async_copy, *p_stream_obj, m_cluster_cfg,
        m_logger->clone("DeviceClusteringAlg"));
    m_meas_sort_alg_gpu = std::make_unique<typename TRAIT::measurement_sorting_algorithm_type>(
        *m_mr, *m_async_copy, *p_stream_obj,
        m_logger->clone("DeviceMeasurementSortingAlg"));
    m_spacepoint_form_alg_gpu = std::make_unique<typename TRAIT::spacepoint_formation_algorithm_type>(
        *m_mr, *m_async_copy, *p_stream_obj);
    m_seeding_alg_gpu = std::make_unique<typename TRAIT::seeding_algorithm_type>(m_seedfinder,
        m_seedfinder, m_seedfilter, *m_mr, *m_async_copy, *p_stream_obj,
        m_logger->clone("DeviceSeedingAlg"));
    m_track_params_estimation_gpu = std::make_unique<typename TRAIT::track_params_estimation_type>(
        *m_mr, *m_async_copy, *p_stream_obj,
        m_logger->clone("DeviceTrackParamEstAlg"));
    m_finding_alg_gpu = std::make_unique<typename TRAIT::finding_algorithm_type>(
        m_finding_cfg, *m_mr, *m_async_copy, *p_stream_obj,
        m_logger->clone("DeviceFindingAlg"));
    m_fitting_alg_gpu = std::make_unique<typename TRAIT::fitting_algorithm_type>(
        m_fitting_cfg, *m_mr, *m_async_copy, *p_stream_obj,
        m_logger->clone("DeviceFittingAlg"));

    m_copy_track_states = std::make_unique<
        traccc::device::container_d2h_copy_alg<
            traccc::track_state_container_types>>(*m_mr, *m_async_copy,
                m_logger->clone("TrackStateD2HCopyAlg"));

    m_device_det_descr = std::make_unique<traccc::silicon_detector_description::buffer>(
        static_cast<traccc::silicon_detector_description::buffer::size_type>(
            m_det_descr->size()),
        *m_device_mr
        );
    typename TRAIT::copy_type::event_type const copy_rv = (*m_copy)(vecmem::get_data(*m_det_descr), *m_device_det_descr);
    if (!copy_rv) {
        throw std::runtime_error("copy detector description to gpu failed");
    }
    m_device_detector = detray::get_buffer(detray::get_data(*m_detector), *m_device_mr, *m_copy);
    m_det_view = detray::get_data(m_device_detector);

    ATH_MSG_INFO("Initializing complete");
    return StatusCode::SUCCESS;
}

// Comparator used for sorting cells. This sorting is one of the assumptions
// made in the clusterization algorithm
struct cell_order {
    bool operator()(const traccc::io::csv::cell& lhs,
                    const traccc::io::csv::cell& rhs) const {
        if (lhs.channel1 != rhs.channel1) {
            return (lhs.channel1 < rhs.channel1);
        } else {
            return (lhs.channel0 < rhs.channel0);
        }
    }
};

// Comparison / ordering operator for measurements
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

template<typename TRAIT> StatusCode TrackingRecoTool<TRAIT>::doRecoFromClusters(
    const EventContext& evtcontext,
    std::vector<clusterInfo>& detray_clusters,
    std::map<int,int>& cluster_map)
{
    ATH_MSG_DEBUG("Reconstructing tracks from clusters: " << detray_clusters.size());
    m_chrono->chronoStart("track reconstruction: clusters to Acts tracks");

    m_chrono->chronoStart("1 cpu cluster sorting");
    std::sort(detray_clusters.begin(), detray_clusters.end(), measurement_sort_comp());
    m_chrono->chronoStop("1 cpu cluster sorting");

    m_chrono->chronoStart("2 read measurements");
    traccc::measurement_collection_types::host measurements_host_buffer{m_host_mr.get()};
    read_measurements(measurements_host_buffer, detray_clusters);
    m_chrono->chronoStop("2 read measurements");

    m_chrono->chronoStart("3 host to gpu data transfer");
    traccc::measurement_collection_types::buffer measurements_gpu_buffer(
        measurements_host_buffer.size(), m_mr->main);
    (*m_async_copy)(vecmem::get_data(measurements_host_buffer),
        measurements_gpu_buffer)->wait();
    m_chrono->chronoStop("3 host to gpu data transfer");

    ATH_CHECK(clustersToTracks(evtcontext, measurements_gpu_buffer, cluster_map));

    m_chrono->chronoStop("track reconstruction: clusters to Acts tracks");

    ATH_MSG_DEBUG("Reconstruction done");
    return StatusCode::SUCCESS;
}

template <typename TRAIT> StatusCode TrackingRecoTool<TRAIT>::doRecoFromHits(
    const EventContext& evtcontext)
{
    // ATH_MSG_INFO("Reconstructing tracks from hits");
    m_chrono->chronoStart("track reconstruction: RDO to Acts tracks");

    // Instantiate host containers/collections
    traccc::edm::silicon_cell_collection::host cells_host_buffer{*m_host_mr};

    // ATH_MSG_INFO("reading cells");
    m_chrono->chronoStart("1 read cells");
    ATH_CHECK(read_cells(cells_host_buffer, evtcontext));
    m_chrono->chronoStop("1 read cells");

    m_cells += cells_host_buffer.size();

    // ATH_MSG_INFO("host to gpu data transfer: " << cells_host_buffer.size() << " bytes");

    // stops in clustersToTracks
    m_chrono->chronoStart("track reconstruction: host memory to host memory");

    m_chrono->chronoStart("2 host to gpu data transfer");
    traccc::edm::silicon_cell_collection::buffer cells_gpu_buffer(
        cells_host_buffer.size(), m_mr->main);
    vecmem::copy::event_type rv = (*m_async_copy)(vecmem::get_data(
        cells_host_buffer), cells_gpu_buffer);
    if (!rv) {
      throw std::runtime_error("cells copy failed");
    }
    m_chrono->chronoStop("2 host to gpu data transfer");

    // ATH_MSG_DEBUG("clusterization #cells: " << m_cells);

    // stops in clustersToTracks
    m_chrono->chronoStart("track reconstruction: gpu work");

    m_chrono->chronoStart("3 gpu clusterization");
    traccc::measurement_collection_types::buffer measurements_gpu_buffer =
        (*m_clusterization_alg_gpu)(cells_gpu_buffer, *m_device_det_descr);
    (*m_meas_sort_alg_gpu)(measurements_gpu_buffer);
    m_chrono->chronoStop("3 gpu clusterization");
    m_measurements_gpu += (*m_copy).get_size(measurements_gpu_buffer);

    std::map<int,int> cluster_map;
    ATH_CHECK(clustersToTracks(evtcontext, measurements_gpu_buffer, cluster_map));

    m_chrono->chronoStop("track reconstruction: RDO to Acts tracks");

    // ATH_MSG_INFO("Reconstruction done");
    return StatusCode::SUCCESS;
}

template <typename TRAIT> StatusCode TrackingRecoTool<TRAIT>::doRecoAndWriteFromHits(
    EventContext const& eventContext)
{
    // Instantiate host containers/collections
    traccc::edm::silicon_cell_collection::host cells_host_buffer{*m_host_mr};

    ATH_MSG_DEBUG("reading cells when truth info is required");
    m_chrono->chronoStart("1 read cells");
    if(read_cells(cells_host_buffer, eventContext).isFailure()){
        throw std::runtime_error("reading of RDO cells failed");
    }
    m_chrono->chronoStop("1 read cells");

    traccc::host::clusterization_algorithm::output_type measurements_host_buffer;
    traccc::host::silicon_pixel_spacepoint_formation_algorithm::output_type spacepoints_host_buffer{*m_cached_host_mr};

    traccc::host::clusterization_algorithm ca(*m_host_mr,m_logger->clone("HostClusteringAlg"));
    traccc::host::sparse_ccl_algorithm cc(*m_host_mr);
    traccc::edm::silicon_cluster_collection::host traccc_clusters = cc(
        vecmem::get_data(cells_host_buffer));

    m_chrono->chronoStart("3 clusterization");
    measurements_host_buffer = ca(vecmem::get_data(cells_host_buffer), vecmem::get_data(*m_det_descr));
    m_chrono->chronoStop("3 clusterization");

    ATH_MSG_DEBUG("spacepoint formation");
    traccc::host::silicon_pixel_spacepoint_formation_algorithm sf(*m_host_mr);

    m_chrono->chronoStart("3.5 cluster container creation");
    std::map<int,int> cluster_map = m_hitInputTool->convertClusters(eventContext, traccc_clusters, measurements_host_buffer, cells_host_buffer);
    m_chrono->chronoStop("3.5 cluster container creation");

    ATH_MSG_DEBUG("Copying measurements over to gpu");
    traccc::measurement_collection_types::buffer measurements_gpu_buffer(
        static_cast<unsigned int>(measurements_host_buffer.size()),
        m_mr->main);
    (*m_async_copy).setup(measurements_gpu_buffer)->wait();
    (*m_async_copy)(vecmem::get_data(measurements_host_buffer),
            measurements_gpu_buffer)->wait();

    ATH_CHECK(clustersToTracks(eventContext, measurements_gpu_buffer, cluster_map));

    ATH_MSG_DEBUG("Reconstruction done");
    return StatusCode::SUCCESS;

}

template<typename TRAIT> StatusCode TrackingRecoTool<TRAIT>::clustersToTracks(
      EventContext const & evtcontext
    , traccc::measurement_collection_types::buffer const & measurements_gpu_buffer
    , std::map<int,int> cluster_map)
{

    // if there are no hits, then no tracking can be done
    // we do not pass empty containers to traccc, but still have to create empty Athena containers
    if((m_async_copy->get_size(measurements_gpu_buffer)) == 0){
        // starts in calling function
        m_chrono->chronoStop("track reconstruction: gpu work");
        // starts in calling function
        m_chrono->chronoStop("track reconstruction: host memory to host memory");

        traccc::track_state_container_types::host track_candidates;
        unsigned nb_output_tracks = 0;
        ATH_CHECK(m_cnvTool->convertTracks(evtcontext, track_candidates,
            cluster_map, nb_output_tracks));

        if(m_checkSeeds){

            traccc::edm::seed_collection::host seeds_host_buffer{*m_cached_host_mr};
            traccc::measurement_collection_types::host measurements_host_buffer{
                m_cached_host_mr.get()};
            traccc::edm::spacepoint_collection::host spacepoints_host_buffer{
                *m_cached_host_mr};
            ATH_CHECK(m_cnvTool->convertSeeds(evtcontext, seeds_host_buffer,
            spacepoints_host_buffer, measurements_host_buffer, cluster_map));
        }

        return StatusCode::SUCCESS;

    }

    // ATH_MSG_INFO("spacepoint formation");
    m_chrono->chronoStart("4 gpu spacepoint formation");
    typename TRAIT::spacepoint_formation_algorithm_type::output_type spacepoints_gpu_buffer =
        (*m_spacepoint_form_alg_gpu)(m_det_view, measurements_gpu_buffer);
    m_chrono->chronoStop("4 gpu spacepoint formation");
    m_spacepoints_gpu += (*m_copy).get_size(spacepoints_gpu_buffer);

    // ATH_MSG_INFO("seeding");
    m_chrono->chronoStart("5 traccc seeding");
    typename TRAIT::seeding_algorithm_type::output_type seeds_gpu_buffer =
        (*m_seeding_alg_gpu)(spacepoints_gpu_buffer);
    m_chrono->chronoStop("5 traccc seeding");
    m_seeds_gpu += (*m_copy).get_size(seeds_gpu_buffer);

    if(m_checkSeeds){

        traccc::edm::seed_collection::host seeds_host_buffer{*m_cached_host_mr};
        traccc::measurement_collection_types::host measurements_host_buffer{
            m_cached_host_mr.get()};
        traccc::edm::spacepoint_collection::host spacepoints_host_buffer{
            *m_cached_host_mr};

        m_chrono->chronoStart("5.5 convert traccc seeds");
        ATH_MSG_DEBUG("Copying over the measurements and seeds from device");
        (*m_async_copy)(seeds_gpu_buffer, seeds_host_buffer)->wait();
        (*m_async_copy)(spacepoints_gpu_buffer, spacepoints_host_buffer)->wait();
        (*m_async_copy)(measurements_gpu_buffer, measurements_host_buffer)->wait();
        ATH_CHECK(m_cnvTool->convertSeeds(evtcontext, seeds_host_buffer,
            spacepoints_host_buffer, measurements_host_buffer, cluster_map));
        m_chrono->chronoStop("5.5 convert traccc seeds");

    }

    // ATH_MSG_INFO("track params estimation");
    m_chrono->chronoStart("6 gpu track params estimation");
    typename TRAIT::track_params_estimation_type::output_type track_params_gpu_buffer =
        (*m_track_params_estimation_gpu)(measurements_gpu_buffer,
            spacepoints_gpu_buffer, seeds_gpu_buffer,
            {0.f, 0.f, m_seedfinder.bFieldInZ});
    m_chrono->chronoStop("6 gpu track params estimation");
    m_params_gpu += (*m_copy).get_size(track_params_gpu_buffer);

    // ATH_MSG_INFO("track finding");
    m_chrono->chronoStart("7 gpu track finding");
    typename TRAIT::finding_algorithm_type::output_type track_candidates_gpu_buffer =
        (*m_finding_alg_gpu)(m_det_view, m_inhom_gpu_field, measurements_gpu_buffer,
            track_params_gpu_buffer);
    m_chrono->chronoStop("7 gpu track finding");
    m_found_tracks_gpu += (*m_copy).get_size(track_candidates_gpu_buffer);

    // ATH_MSG_INFO("track fitting");
    m_chrono->chronoStart("8 gpu track fitting");
    typename TRAIT::fitting_algorithm_type::output_type track_states_gpu_buffer =
        (*m_fitting_alg_gpu)(m_det_view, m_inhom_gpu_field,
        {track_candidates_gpu_buffer, measurements_gpu_buffer});
    m_chrono->chronoStop("8 gpu track fitting");
    // m_fitted_tracks_gpu += (*m_copy).get_size(track_states_gpu_buffer);
    m_fitted_tracks_gpu += track_states_gpu_buffer.items.size();

    // starts in calling function
    m_chrono->chronoStop("track reconstruction: gpu work");

    // ATH_MSG_INFO("copying data from gpu to host");
    m_chrono->chronoStart("9 gpu to host data transfer");
    auto track_states_gpu = (*m_copy_track_states)(track_states_gpu_buffer);
    m_chrono->chronoStop("9 gpu to host data transfer");

    // starts in calling function
    m_chrono->chronoStop("track reconstruction: host memory to host memory");

    // ATH_MSG_INFO("converting traccc tracks to AOD");
    m_chrono->chronoStart("10 traccc tracks conversion to Acts tracks");
    unsigned nb_output_tracks = 0;
    ATH_CHECK(m_cnvTool->convertTracks(evtcontext, track_states_gpu,
        cluster_map, nb_output_tracks));
    m_output_tracks += nb_output_tracks;
    m_chrono->chronoStop("10 traccc tracks conversion to Acts tracks");

    return StatusCode::SUCCESS;
}

template<typename TRAIT> void TrackingRecoTool<TRAIT>::read_measurements(
    traccc::measurement_collection_types::host& measurements,
    std::vector<clusterInfo>& detray_clusters)
{
  for(std::vector<clusterInfo>::size_type i = 0; i < detray_clusters.size();i++){

    clusterInfo cluster = detray_clusters[i];

    uint64_t geometry_id = cluster.detray_id;
    const auto& sf = detray::geometry::barcode{geometry_id};
    const detray::tracking_surface surface{*m_detector, sf};
    cluster.localPosition[0] = cluster.localPosition.x();
    cluster.localPosition[1] = cluster.localPosition.y();

    ATH_MSG_DEBUG("Traccc measurement at index " << i << ": " << cluster.localPosition[0] << "," << cluster.localPosition[1]);

    // Construct the measurement object.
    traccc::measurement meas;
    std::array<detray::dsize_type<traccc::default_algebra>, 2u> indices{0u, 0u};
    meas.meas_dim = 0u;
    for (unsigned int ipar = 0; ipar < 2u; ++ipar) {
        if (((cluster.local_key) & (1 << (ipar + 1))) != 0) {

            switch (ipar) {
                case 0: {
                    meas.local[0] = cluster.localPosition.x();
                    meas.variance[0] = 0.0025;
                    indices[meas.meas_dim++] = ipar;
                }; break;
                case 1: {
                    meas.local[1] = cluster.localPosition.y();
                    meas.variance[1] = 0.0025;
                    indices[meas.meas_dim++] = ipar;
                }; break;
            }
	    }
    }

    meas.subs.set_indices(indices);
    meas.surface_link = detray::geometry::barcode{geometry_id};

    // Keeps measurement_id for ambiguity resolution
    meas.measurement_id = cluster.measurement_id;
    measurements.push_back(meas);
  }
}


template <typename TRAIT> StatusCode TrackingRecoTool<TRAIT>::read_cells(
    traccc::edm::silicon_cell_collection::host& cells,
    const EventContext& evtcontext)
{
    std::map<std::uint64_t, std::vector<traccc::io::csv::cell> > result;

    ATH_MSG_DEBUG("Reading pixel hits");
    auto pixelRDOHandle = SG::makeHandle(m_pixelRDOKey, evtcontext);
    int nPix = 0;
    for (const InDetRawDataCollection<PixelRDORawData>* pixel_rdoCollection : *pixelRDOHandle) {
        if (pixel_rdoCollection == nullptr) { continue; }

        for (const PixelRDORawData* pixelRawData : *pixel_rdoCollection) {

            Identifier rdoId = pixelRawData->identify();

            // get the det element from the det element collection
            const InDetDD::SiDetectorElement* sielement = m_pixelManager->getDetectorElement(rdoId); assert(sielement);
            const Identifier Pixel_ModuleID = sielement->identify();
            InDetDD::SiCellId id = sielement->cellIdFromIdentifier(rdoId);

            uint64_t const geometry_id = findDetrayID(Pixel_ModuleID, m_AthenaToDetrayMap);

            ATH_MSG_DEBUG("Doing this module: " << geometry_id);

            result[geometry_id].push_back({
                geometry_id,
                0,
                static_cast<uint32_t>(id.phiIndex()),
                static_cast<uint32_t>(id.etaIndex()),
                static_cast<float>(pixelRawData->getToT()),
                8 // timestamp is not used
            });

            nPix++;
        }

    }
    ATH_MSG_DEBUG("Read " << nPix << " pixel hits");

    ATH_MSG_DEBUG("Reading strip hits");
    int nStrip = 0;
    auto stripRDOHandle = SG::makeHandle(m_stripRDOKey, evtcontext);
    for (const InDetRawDataCollection<SCT_RDORawData>* strip_Collection : *stripRDOHandle) {
        if (strip_Collection == nullptr) { continue; }
        for (const SCT_RDORawData* stripRawData : *strip_Collection) {

            const Identifier rdoId = stripRawData->identify();
            const InDetDD::SiDetectorElement* sielement = m_stripManager->getDetectorElement(rdoId);

            const Identifier strip_moduleID = m_stripID->module_id(sielement->identify()); //from wafer id to module id
            const IdentifierHash Strip_ModuleHash = m_stripID->wafer_hash(strip_moduleID);

            // Extract the correct Strip_ModuleID
            int side = m_stripID->side(sielement->identify());
            const Identifier Strip_ModuleID = m_stripID->wafer_id(Strip_ModuleHash + side);

            InDetDD::SiCellId id = sielement->cellIdFromIdentifier(rdoId);

            uint64_t const geometry_id = findDetrayID(Strip_ModuleID, m_AthenaToDetrayMap);

            for (int i = 0; i < stripRawData->getGroupSize(); i++) {

                result[geometry_id].push_back({
                    geometry_id,
                    0,
                    static_cast<uint32_t>(id.phiIndex()+i),
                    0,
                    1,
                    8 // timestamp is not used
                });
                nStrip++;
            }

        }

    }
    ATH_MSG_DEBUG("Read " << nStrip << " strip hits");

    // Sort the cells. Deduplication or not, they do need to be sorted.
    for (auto& [_, hits] : result) {
        std::sort(hits.begin(), hits.end(), ::cell_order());
    }

    ATH_MSG_DEBUG("Sorted the cells container");

    // If there is a detector description object, build a map of geometry IDs
    // to indices inside the detector description. We do not use ACTS IDs, so we map to detray IDs
    std::map<uint64_t, unsigned int> geomIdMap;

    for (unsigned int i = 0; i < m_det_descr->acts_geometry_id().size(); ++i) {
        geomIdMap[m_det_descr->geometry_id()[i].value()] = i;
    }

    // Fill the output containers with the ordered cells and modules.
    for (const auto& [geometry_id, cellz] : result) {

        // Figure out the index of the detector description object, for this
        // group of cells.
        unsigned int ddIndex = 0;

        auto it = geomIdMap.find(geometry_id);
        if (it == geomIdMap.end()) {
            throw std::runtime_error("Could not find geometry ID (" +
                                        std::to_string(geometry_id) +
                                        ") in the detector description");
        }

        ddIndex = it->second;

        // Add the cells to the output.
        for (auto& cell : cellz) {
            cells.push_back({cell.channel0, cell.channel1, cell.value,
                             cell.timestamp, ddIndex});
        }
    }

    return StatusCode::SUCCESS;
}

template <typename TRAIT> std::vector<hitInfo>
TrackingRecoTool<TRAIT>::make_fake_data(
    traccc::silicon_detector_description::host& dd)
{
    std::ofstream t_cells_file("atlas_hits_diagonal.txt",std::ofstream::out);

    std::vector<uint64_t> v = {256986819764241407, 256986822448597247,
        265994027072050047, 275001274376760703};
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

        ////////////////////////////////////////
        // make fake data at the module corners
        ///////////////////////////////////////

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

        /////////////////////////////////////
        // make fake data along the diagonal
        /////////////////////////////////////

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

template <typename TRAIT> StatusCode TrackingRecoTool<TRAIT>::finalize()
{
    ATH_MSG_INFO("==> Statistics of traccc track reconstruction:");
    // ATH_MSG_INFO("- read from      " << m_modules << " modules");
    ATH_MSG_INFO("- read           " << m_cells << " cells");
    // ATH_MSG_INFO("- created (cpu)  " << m_measurements_cpu << " measurements");
    ATH_MSG_INFO("- created (gpu)  " << m_measurements_gpu << " measurements");
    // ATH_MSG_INFO("- created (cpu)  " << m_spacepoints_cpu << " spacepoints");
    ATH_MSG_INFO("- created (gpu)  " << m_spacepoints_gpu << " spacepoints");
    // ATH_MSG_INFO("- created (cpu)  " << m_seeds_cpu << " seeds");
    ATH_MSG_INFO("- created (gpu)  " << m_seeds_gpu << " seeds");
    // ATH_MSG_INFO("- created (cpu)  " << m_params_cpu << " initial parameters");
    ATH_MSG_INFO("- created (gpu)  " << m_params_gpu << " initial parameters");
    // ATH_MSG_INFO("- found (cpu)    " << m_found_tracks_cpu << " tracks");
    ATH_MSG_INFO("- found (gpu)    " << m_found_tracks_gpu << " tracks");
    // ATH_MSG_INFO("- fitted (cpu)   " << m_fitted_tracks_cpu << " tracks");
    ATH_MSG_INFO("- fitted (gpu)   " << m_fitted_tracks_gpu << " tracks");
    ATH_MSG_INFO("- output (gpu)   " << m_output_tracks << " tracks");
    // ATH_MSG_INFO("- resolved (cpu) " << m_ambiguity_free_tracks_cpu << " tracks");
    // ATH_MSG_INFO("- resolved (gpu) " << m_ambiguity_free_tracks_gpu << " tracks");

    return StatusCode::SUCCESS;
}

template <typename TRAIT> void TrackingRecoTool<TRAIT>::write_traccc_data(
    int m_n_event, traccc::edm::spacepoint_collection::host& spacepoints,
    traccc::measurement_collection_types::host& measurements)
{
    std::string padded = std::format("{:09}",m_n_event);

    std::ofstream t_ms_file("event"+padded+"-traccc-measurements.txt",std::ofstream::out);
    std::ofstream t_sp_file("event"+padded+"-traccc-spacepoints.txt",std::ofstream::out);

    for(long unsigned int i = 0; i < spacepoints.size(); i++){

        auto sp = spacepoints[i];
        unsigned int meas_id = sp.measurement_index_1();
        const auto& meas = measurements.at(meas_id);
        auto local = meas.local;
        auto global = sp.global();
        auto id = meas.surface_link.value();

        t_sp_file << id << "," << global[0] << "," << global[1] << "," << global[2] << "\n";
        t_ms_file << id << "," << local[0] << "," << local[1] << "\n";
    }

    t_ms_file.close();
    t_sp_file.close();
}

template <typename TRAIT> StatusCode
TrackingRecoTool<TRAIT>::makeTracccStandaloneData(
    int m_n_event,const EventContext& evtcontext,
    std::string const& outDir)
{

    if(!std::filesystem::exists(outDir)){
        ATH_MSG_FATAL("You did not provide a valid file path for traccc standalone file writing!");
    }

    std::string padded = std::format("{:09}",m_n_event);

    traccc::edm::silicon_cell_collection::host cells_per_event{*m_host_mr};
    ATH_CHECK(read_cells(cells_per_event, evtcontext));

    traccc::io::write(m_n_event,outDir,traccc::data_format::csv, vecmem::get_data(cells_per_event),vecmem::get_data(*m_det_descr),false);

    return StatusCode::SUCCESS;
}

template <typename TRAIT> StatusCode
TrackingRecoTool<TRAIT>::makeTracccStandaloneData(
    int m_n_event, std::vector<clusterInfo>& detray_clusters, std::string const& outDir)
{

    if(!std::filesystem::exists(outDir)){
        ATH_MSG_FATAL("You did not provide a valid file path for traccc standalone file writing!");
    }

    ATH_MSG_INFO("Converting event number " << m_n_event << " to Traccc standalone data.");
    std::string padded = std::format("{:09}",m_n_event);

    std::ofstream hitfile;
    hitfile.open(outDir+"/event"+padded+"-hits.csv");
    hitfile << "particle_id,geometry_id,tx,ty,tz,tt,tpx,tpy,tpz,te,deltapx,deltapy,deltapz,deltae,index" << "\n";

    std::ofstream measfile;
    measfile.open(outDir+"/event"+padded+"-measurements.csv");
    measfile << "measurement_id,geometry_id,local_key,local0,local1,phi,theta,time,var_local0,var_local1,var_phi,var_theta,var_time" << "\n";

    std::ofstream mapfile;
    mapfile.open(outDir+"/event"+padded+"-measurement-simhit-map.csv");
    mapfile << "measurement_id,hit_id" << "\n";

    for(std::vector<clusterInfo>::size_type i = 0; i < detray_clusters.size(); i++){

        clusterInfo &thiscluster = detray_clusters[i];
        if(thiscluster.detray_id == 0){ATH_MSG_INFO("You have an unmatched cluster! " << i);}

        hitfile << i << "," << thiscluster.detray_id << "," << thiscluster.globalPosition[0] << "," << thiscluster.globalPosition[1] << "," << thiscluster.globalPosition[2];
        hitfile << "," << "0,0,0,0,0,0,0,0,0,0" << "\n";

        measfile << i << "," << thiscluster.detray_id << ",6," << thiscluster.localPosition[0] << "," << thiscluster.localPosition[1] << ",";
        measfile << "0,0,0,0.0025,0.0025,0,0,0" << "\n";
        mapfile << i << "," << i << "\n";

    }

    hitfile.close();
    measfile.close();
    mapfile.close();

    return StatusCode::SUCCESS;
}

template <typename TRAIT> StatusCode
TrackingRecoTool<TRAIT>::makeTruthSeeds(traccc::bound_track_parameters_collection_types::host& truth_seeds, traccc::measurement_collection_types::host& measurements_host_buffer)
{
    SG::ReadHandle<ActsTrk::MeasurementToTruthParticleAssociation> pixelClustersToTruthAssociation =
      SG::makeHandle(m_pixelClustersToTruth, Gaudi::Hive::currentContext());
    if (!pixelClustersToTruthAssociation.isValid()) {
        ATH_MSG_ERROR("No pixel clusters for key " << m_pixelClustersToTruth.key());
        //return StatusCode::FAILURE;
    }

    ATH_MSG_INFO("Getting truth particles associated to pixel clusters");

    SG::ReadHandle<xAOD::PixelClusterContainer> inputPixelClusterContainer = SG::makeHandle(
        m_xAODPixelClusterFromInDetClusterKey, Gaudi::Hive::currentContext());
    if (!inputPixelClusterContainer.isValid()) {
        ATH_MSG_FATAL("xAOD::PixelClusterContainer with key " << m_xAODPixelClusterFromInDetClusterKey.key()
                        << " is not available...");
    }
    const xAOD::PixelClusterContainer *inputPixelClusters = inputPixelClusterContainer.cptr();

    // make a map of stable truth particles (barcode) and first cluster associated to it
    // we loop over all the clusters and only collect the lowest cluster per truth particle
    std::map<int, const xAOD::PixelCluster*> barcode_cluster_map;
    std::map<int, const xAOD::TruthParticle*> barcode_tp_map;


    for (unsigned int idx = 0; idx < inputPixelClusters->size(); idx++) {
        const xAOD::PixelCluster *inputPixelCluster = inputPixelClusters->at(idx);
        for (const xAOD::TruthParticle *truth_particle : pixelClustersToTruthAssociation->at(idx)) {
            if (!truth_particle) continue;
            //int bc = truth_particle->barcode();
            // select stable, charged particles
            if(!(truth_particle->isCharged() && truth_particle->isStable())){continue;}
            if (barcode_cluster_map.find((*truth_particle).barcode()) != barcode_cluster_map.end()) {
                ATH_MSG_DEBUG("Found this truth particle before: " << (*truth_particle).barcode());
                const xAOD::PixelCluster* prev_cl =  barcode_cluster_map[(*truth_particle).barcode()];
                ATH_MSG_DEBUG("Getting cluster position");
                auto prev_globalPos = prev_cl->globalPosition();
                float prev_dist = std::sqrt(std::pow(static_cast<float>(prev_globalPos.x()),2) +
                                    std::pow(static_cast<float>(prev_globalPos.y()),2) +
                                    std::pow(static_cast<float>(prev_globalPos.z()),2) );

                auto this_globalPos = inputPixelCluster->globalPosition();
                float this_dist = std::sqrt(std::pow(static_cast<float>(this_globalPos.x()),2) +
                                    std::pow(static_cast<float>(this_globalPos.y()),2) +
                                    std::pow(static_cast<float>(this_globalPos.z()),2) );
                ATH_MSG_DEBUG("Old cluster " << prev_globalPos.x() << "," << prev_globalPos.y() << "," << prev_globalPos.z() << ", dist " << prev_dist << " and new one " << this_globalPos.x() << "," << this_globalPos.y() << "," << this_globalPos.z() << " dist: " << prev_dist << "/" << this_dist);
                if(this_dist < prev_dist){
                    ATH_MSG_DEBUG("Will replace the old cluster");
                    barcode_cluster_map[(*truth_particle).barcode()] = inputPixelCluster;
                }
            }else{
                ATH_MSG_DEBUG("Adding this truth particle: " << (*truth_particle).barcode());
                barcode_cluster_map[(*truth_particle).barcode()] = inputPixelCluster;
                barcode_tp_map[(*truth_particle).barcode()] = truth_particle;

            }
        }
    }

    // Standard deviations for seed track parameters
    static constexpr std::array<traccc::scalar, traccc::e_bound_size> stddevs =
        {1e-4f * traccc::unit<traccc::scalar>::mm,
            1e-4f * traccc::unit<traccc::scalar>::mm,
            1e-3f,
            1e-3f,
            1e-4f / traccc::unit<traccc::scalar>::GeV,
            1e-4f * traccc::unit<traccc::scalar>::ns};

    traccc::seed_generator<traccc::default_detector::host> sg(*m_detector,
                                                                stddevs);

    // loop over the map and create initial point and starting parameters
    ATH_MSG_INFO("Constructing " << barcode_cluster_map.size() << " truth seeds: initial parameters + first measurement");
    for(const auto& [barcode, cluster] : barcode_cluster_map){

        auto idHash = cluster->identifierHash();
        traccc::measurement meas;

        // find the measurement that corresponds to this cluster
        for(long unsigned int m = 0; m < measurements_host_buffer.size(); m++){

            const auto& this_meas = measurements_host_buffer.at(m);
            const std::uint_least64_t detray_id = this_meas.surface_link.value();
            Identifier const atlas_ID = findAthenaID(detray_id, m_DetrayToAthenaMap);
            const IdentifierHash Pixel_ModuleHash = m_pixelID->wafer_hash(atlas_ID);
            if(idHash != Pixel_ModuleHash){continue;}
            auto localPos = cluster->localPosition<2>();
            float dist = std::sqrt( std::pow(static_cast<float>(localPos.x()) - static_cast<float>(this_meas.local[0]),2)+ std::pow(static_cast<float>(localPos.y()) - static_cast<float>(this_meas.local[1]),2));
            if(dist==0){
                meas = this_meas;
            }
        }

        const xAOD::TruthParticle* tp = barcode_tp_map[barcode];
        ATH_MSG_DEBUG("This truth particle info: pdgID " << tp->pdgId() << " charge " << tp->charge());
        const detray::free_track_parameters<traccc::default_algebra> free_param(
            {cluster->globalPosition().x(),cluster->globalPosition().y(),cluster->globalPosition().z()},
            0.f, {tp->px(),tp->py(),tp->pz()}, tp->charge());
        auto seed_params = sg(meas.surface_link, free_param,
                            traccc::detail::particle_from_pdg_number<traccc::scalar>(tp->pdgId()));
        truth_seeds.push_back(seed_params);
    }

    ATH_MSG_DEBUG("Done truth seeding: " << truth_seeds.size());

  return StatusCode::SUCCESS;
}

#endif

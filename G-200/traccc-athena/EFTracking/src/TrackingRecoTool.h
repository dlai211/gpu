/*
  Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration
*/

#ifndef EFTRACKING_TRACKINGRECOTOOL_H
#define EFTRACKING_TRACKINGRECOTOOL_H

// Athena includes
#include "AthenaBaseComps/AthAlgTool.h"
#include "GaudiKernel/ServiceHandle.h"

// Project include(s).
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/cuda/clusterization/clusterization_algorithm.hpp"
#include "traccc/clusterization/clustering_config.hpp"
#include "traccc/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
#include "traccc/cuda/clusterization/measurement_sorting_algorithm.hpp"
#include "traccc/cuda/seeding/spacepoint_formation_algorithm.hpp"
#include "traccc/cuda/finding/finding_algorithm.hpp"
#include "traccc/cuda/fitting/fitting_algorithm.hpp"
#include "traccc/cuda/seeding/seeding_algorithm.hpp"
#include "traccc/cuda/seeding/track_params_estimation.hpp"
#include "traccc/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/device/container_d2h_copy_alg.hpp"
#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/io/read_geometry.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"
#include "traccc/io/read_geometry.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/read_spacepoints.hpp"
#include "traccc/io/read_measurements.hpp"
#include "traccc/io/read_cells.hpp"
#include "traccc/io/write.hpp"
#include "traccc/geometry/detector.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// Detray include(s).
#include "detray/core/detector.hpp"
#include "detray/core/detector_metadata.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"

#include "EFTracking/ITrackingRecoTool.h"

class TrackingRecoTool : public extends<AthAlgTool, ITrackingRecoTool> {

  public:
        TrackingRecoTool(const std::string&, const std::string&, const IInterface*);
        virtual ~TrackingRecoTool() { ; }

        virtual StatusCode initialize() override;
        virtual StatusCode finalize()   override;
        virtual traccc::track_state_container_types::host doRecoFromClusters(std::vector<clusterInfo>& detray_clusters)  override;
        virtual traccc::track_state_container_types::host doRecoFromHits(std::vector<hitInfo>& detray_hits,std::vector<clusterInfo>& detray_clusters)  override;
        

  private:
         Gaudi::Property<bool> m_doRecoOnCPU{this,"recoOnCPU",true,"If you also want to perform track reconstruction on CPU"};
         Gaudi::Property<bool> m_doWriteData{this,"writeData",true,"If you also want to create input for traccc standalone"};

         int m_n_event;
         std::map<std::uint64_t, int> m_DetrayToAtlasMap;
         std::map<int, std::uint64_t> m_AtlasToDetrayMap;

         virtual traccc::track_state_container_types::host recoClustersToTracks(traccc::spacepoint_collection_types::host& spacepoints_per_event, traccc::measurement_collection_types::host& measurements_per_event);
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
         float m_bFieldInZ;

         void read_spacepoints(traccc::spacepoint_collection_types::host& spacepoints, std::vector<clusterInfo>& detray_clusters, bool do_strip);
         void read_measurements(traccc::measurement_collection_types::host& measurements, std::vector<clusterInfo>& detray_clusters, bool do_strip);
         void read_hits(traccc::edm::silicon_cell_collection::host& cells,std::vector<hitInfo>& detray_hits);
         void read_cells(traccc::edm::silicon_cell_collection::host& cells, std::vector<hitInfo>& detray_hits);

         std::vector<hitInfo> make_fake_data(traccc::silicon_detector_description::host& dd);
         void write_traccc_data(int m_n_event,traccc::spacepoint_collection_types::host& spacepoints);
         void make_test_data(int m_n_event,std::vector<clusterInfo>& detray_clusters);
         
         
         // Memory resources used by the application.
         vecmem::host_memory_resource host_mr;
         vecmem::cuda::host_memory_resource cuda_host_mr;
         vecmem::cuda::managed_memory_resource mng_mr;
         vecmem::cuda::device_memory_resource device_mr;
         traccc::memory_resource mr{device_mr, &cuda_host_mr};

         // CUDA types used.
         traccc::cuda::stream stream;
         vecmem::cuda::copy copy;
         vecmem::cuda::async_copy async_copy{stream.cudaStream()};

         using host_detector_type = detray::detector<>;
         host_detector_type::buffer_type device_detector;
         host_detector_type::view_type det_view;
         host_detector_type m_detector{host_mr};
         traccc::silicon_detector_description::host m_dd{host_mr};


         uint64_t m_modules = 0;
         uint64_t m_cells = 0;
         uint64_t m_measurements = 0;
         uint64_t m_measurements_cuda = 0;
         uint64_t m_spacepoints = 0;
         uint64_t m_spacepoints_cuda = 0;
         uint64_t m_seeds = 0;
         uint64_t m_seeds_cuda = 0;
         uint64_t m_params = 0;
         uint64_t m_params_cuda = 0;
         uint64_t m_found_tracks = 0;
         uint64_t m_found_tracks_cuda = 0;
         uint64_t m_fitted_tracks = 0;
         uint64_t m_fitted_tracks_cuda = 0;
         uint64_t m_ambiguity_free_tracks = 0;

};

#endif

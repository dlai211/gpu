#ifndef EFTRACKING_TRACKINGRECOTOOLCUDA_H
#define EFTRACKING_TRACKINGRECOTOOLCUDA_H

#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/copy.hpp>

#include "traccc/geometry/detector.hpp"

#include "traccc/cuda/clusterization/clusterization_algorithm.hpp"
#include "traccc/cuda/clusterization/measurement_sorting_algorithm.hpp"
#include "traccc/cuda/seeding/spacepoint_formation_algorithm.hpp"
#include "traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/cuda/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/cuda/seeding/seeding_algorithm.hpp"
#include "traccc/cuda/seeding/track_params_estimation.hpp"

#include "TrackingRecoTool.h"

struct TrackingRecoToolTraitCUDA
{
    typedef vecmem::cuda::host_memory_resource gpu_host_memory_resource_type;
    typedef vecmem::cuda::device_memory_resource device_memory_resource_type;
    typedef traccc::cuda::stream stream_type;
    typedef traccc::cuda::stream stream_for_algo_type;
    typedef vecmem::cuda::copy copy_type;
    typedef vecmem::cuda::async_copy async_copy_type;
    typedef traccc::cuda::clusterization_algorithm clusterization_algorithm_type;
    typedef traccc::cuda::measurement_sorting_algorithm measurement_sorting_algorithm_type;
    typedef traccc::cuda::spacepoint_formation_algorithm<traccc::default_detector::device>
      spacepoint_formation_algorithm_type;
    typedef traccc::cuda::seeding_algorithm seeding_algorithm_type;
    typedef traccc::cuda::track_params_estimation track_params_estimation_type;
    typedef traccc::cuda::combinatorial_kalman_filter_algorithm finding_algorithm_type;
    typedef traccc::cuda::kalman_fitting_algorithm fitting_algorithm_type;
};

typedef TrackingRecoTool<TrackingRecoToolTraitCUDA> TrackingRecoToolCUDA;

#endif

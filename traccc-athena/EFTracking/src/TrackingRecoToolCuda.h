#ifndef EFTRACKING_TRACKINGRECOTOOLCUDA_H
#define EFTRACKING_TRACKINGRECOTOOLCUDA_H

#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/copy.hpp>

#include "detray/propagator/rk_stepper.hpp"
// #include "detray/detectors/bfield.hpp"

#include "traccc/geometry/detector.hpp"
#include "traccc/fitting/kalman_filter/kalman_fitter.hpp"

#include "traccc/cuda/clusterization/clusterization_algorithm.hpp"
#include "traccc/cuda/clusterization/measurement_sorting_algorithm.hpp"
#include "traccc/cuda/seeding/spacepoint_formation_algorithm.hpp"
#include "traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/cuda/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/cuda/seeding/seeding_algorithm.hpp"
#include "traccc/cuda/seeding/track_params_estimation.hpp"

#include "TrackingRecoTool.h"

using scalar_type = traccc::default_detector::device::scalar_type;
using b_field_t = covfie::field<traccc::const_bfield_backend_t<scalar_type>>;
using device_detector_type = traccc::default_detector::device;
using rk_stepper_type = detray::rk_stepper<b_field_t::view_t,
  traccc::default_algebra, detray::constrained_step<scalar_type>>;
using device_navigator_type = detray::navigator<const device_detector_type>;
using device_fitter_type = traccc::kalman_fitter<rk_stepper_type, device_navigator_type>;

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

#ifndef EFTRACKING_TRACKINGRECOTOOLSYCL_H
#define EFTRACKING_TRACKINGRECOTOOLSYCL_H

#include <vecmem/memory/sycl/device_memory_resource.hpp>
#include <vecmem/memory/sycl/host_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/sycl/async_copy.hpp>
#include <vecmem/utils/sycl/copy.hpp>
#include <vecmem/utils/sycl/queue_wrapper.hpp>

#include <traccc/sycl/utils/queue_wrapper.hpp>
#include <traccc/sycl/clusterization/clusterization_algorithm.hpp>
#include <traccc/sycl/clusterization/measurement_sorting_algorithm.hpp>
#include <traccc/sycl/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp>
#include <traccc/sycl/seeding/seeding_algorithm.hpp>
#include <traccc/sycl/seeding/track_params_estimation.hpp>
#include <traccc/sycl/finding/combinatorial_kalman_filter_algorithm.hpp>
#include <traccc/sycl/fitting/kalman_fitting_algorithm.hpp>

#include "TrackingRecoTool.h"

struct TrackingRecoToolTraitSYCL
{
    typedef vecmem::sycl::host_memory_resource gpu_host_memory_resource_type;
    typedef vecmem::sycl::device_memory_resource device_memory_resource_type;
    typedef vecmem::sycl::queue_wrapper stream_type;
    typedef traccc::sycl::queue_wrapper stream_for_algo_type;
    typedef vecmem::sycl::copy copy_type;
    typedef vecmem::sycl::async_copy async_copy_type;
    typedef traccc::sycl::clusterization_algorithm clusterization_algorithm_type;
    typedef traccc::sycl::measurement_sorting_algorithm measurement_sorting_algorithm_type;
    typedef traccc::sycl::silicon_pixel_spacepoint_formation_algorithm
        spacepoint_formation_algorithm_type;
    typedef traccc::sycl::seeding_algorithm seeding_algorithm_type;
    typedef traccc::sycl::track_params_estimation track_params_estimation_type;
    typedef traccc::sycl::combinatorial_kalman_filter_algorithm finding_algorithm_type;
    typedef traccc::sycl::kalman_fitting_algorithm fitting_algorithm_type;
};

typedef TrackingRecoTool<TrackingRecoToolTraitSYCL> TrackingRecoToolSYCL;

#endif

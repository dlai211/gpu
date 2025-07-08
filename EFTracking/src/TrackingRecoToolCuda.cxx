#include "TrackingRecoToolCuda.h"
#include "TrackingRecoTool.h"
#include <memory>
#include "traccc/cuda/utils/make_bfield.hpp"

template<> void TrackingRecoTool<TrackingRecoToolTraitCUDA>::create_copy_object(
    std::unique_ptr<TrackingRecoToolTraitCUDA::copy_type> & copy_obj)
{
    copy_obj = std::make_unique<TrackingRecoToolTraitCUDA::copy_type>();
}

template<> void TrackingRecoTool<TrackingRecoToolTraitCUDA>::create_async_copy_object(
    std::unique_ptr<TrackingRecoToolTraitCUDA::async_copy_type> & async_copy_obj)
{
    async_copy_obj = std::make_unique<TrackingRecoToolTraitCUDA::async_copy_type>(
        m_stream->cudaStream());
}

template<> void TrackingRecoTool<TrackingRecoToolTraitCUDA>::get_stream_for_algo(
    TrackingRecoToolTraitCUDA::stream_for_algo_type* & pobj)
{
    pobj = m_stream.get();
}

template<> std::string
    TrackingRecoTool<TrackingRecoToolTraitCUDA>::get_init_message()
{
    return "Initializing reco tool: CUDA implementation";
}

template<> traccc::bfield
    TrackingRecoTool<TrackingRecoToolTraitCUDA>::make_gpu_bfield(
          traccc::bfield const & host_bfield
        , TrackingRecoToolTraitCUDA::stream_for_algo_type &
        )
{
    return traccc::cuda::make_bfield(host_bfield);
}

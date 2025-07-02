#include "TrackingRecoToolSycl.h"
#include "TrackingRecoTool.h"
#include <memory>

template<> void TrackingRecoTool<TrackingRecoToolTraitSYCL>::create_copy_object(
    std::unique_ptr<TrackingRecoToolTraitSYCL::copy_type> & copy_obj)
{
    copy_obj = std::make_unique<TrackingRecoToolTraitSYCL::copy_type>(*m_stream);
}

template<> void TrackingRecoTool<TrackingRecoToolTraitSYCL>::create_async_copy_object(
    std::unique_ptr<TrackingRecoToolTraitSYCL::async_copy_type> & async_copy_obj)
{
    async_copy_obj = std::make_unique<TrackingRecoToolTraitSYCL::async_copy_type>(*m_stream);
}

template<> void TrackingRecoTool<TrackingRecoToolTraitSYCL>::get_stream_for_algo(
    TrackingRecoToolTraitSYCL::stream_for_algo_type* & pobj)
{
    static traccc::sycl::queue_wrapper tqw{m_stream->queue()};
    pobj = &tqw;
}

template<> std::string
    TrackingRecoTool<TrackingRecoToolTraitSYCL>::get_init_message()
{
    return "Initializing reco tool: SYCL implementation";
}

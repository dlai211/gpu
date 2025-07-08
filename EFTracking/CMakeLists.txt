# Copyright (C) 2023-2025 CERN for the benefit of the ATLAS collaboration

# Set the package name.
atlas_subdir( EFTracking )

# ----
# Build a component library with an example algorithm that uses traccc.
# ----

if(BUILD_CUDA)
  atlas_add_component( EFTrackingCUDA
    src/TrackingRecoTool.h
    src/TrackingRecoToolCuda.h
    src/TrackingRecoToolCuda.cxx
    src/TrackConversionTool.h
    src/TrackConversionTool.cxx
    src/TrackingAlg.h
    src/TrackingAlg.cxx
    src/TrackingHitInputTool.h
    src/TrackingHitInputTool.cxx
    src/components/EFTrackingEntriesCuda.cxx
    src/AthenaDetrayConversion.h
    src/AthenaDetrayConversion.cxx
    PUBLIC_HEADERS EFTracking
    DEFINITIONS ALGEBRA_PLUGINS_INCLUDE_ARRAY
    LINK_LIBRARIES vecmem::core vecmem::cuda traccc::core traccc::cuda traccc::io
      detray::io detray::core traccc::device_common GaudiKernel AthenaBaseComps
      xAODEventInfo xAODTracking InDetPrepRawData InDetRawData
      InDetReadoutGeometry ReadoutGeometryBase GeoPrimitives
      PixelReadoutGeometryLib ActsToolInterfacesLib ActsToolInterfacesLib ActsCore
      ActsEventLib)
endif()

if(BUILD_SYCL)
  atlas_add_component( EFTrackingSYCL
    src/TrackingRecoTool.h
    src/TrackingRecoToolSycl.h
    src/TrackingRecoToolSycl.cxx
    src/TrackConversionTool.h
    src/TrackConversionTool.cxx
    src/TrackingAlg.h
    src/TrackingAlg.cxx
    src/TrackingHitInputTool.h
    src/TrackingHitInputTool.cxx
    src/components/EFTrackingEntriesSycl.cxx
    src/AthenaDetrayConversion.h
    src/AthenaDetrayConversion.cxx
    PUBLIC_HEADERS EFTracking
    DEFINITIONS ALGEBRA_PLUGINS_INCLUDE_ARRAY
    LINK_LIBRARIES vecmem::core vecmem::sycl traccc::core traccc::sycl traccc::io
      detray::io detray::core traccc::device_common GaudiKernel AthenaBaseComps
      xAODEventInfo xAODTracking InDetPrepRawData InDetRawData
      InDetReadoutGeometry ReadoutGeometryBase GeoPrimitives
      PixelReadoutGeometryLib ActsToolInterfacesLib ActsToolInterfacesLib ActsCore
      ActsEventLib)
endif()

# Install files from the package:
atlas_install_python_modules( python/*.py POST_BUILD_CMD ${ATLAS_FLAKE8} )

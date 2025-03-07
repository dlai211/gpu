# Copyright (C) 2023 CERN for the benefit of the ATLAS collaboration

# Set the package name.
atlas_subdir( EFTracking )

# Build a component library with an example algorithm that uses traccc.
atlas_add_component( EFTracking
   src/*.h src/*.cxx src/components/*.cxx
   PUBLIC_HEADERS EFTracking
   DEFINITIONS ALGEBRA_PLUGINS_INCLUDE_ARRAY
   LINK_LIBRARIES vecmem::core vecmem::cuda traccc::core traccc::cuda traccc::io detray::io detray::core traccc::device_common traccc::cuda
                  GaudiKernel AthenaBaseComps xAODEventInfo xAODTracking 
                  InDetPrepRawData InDetRawData InDetReadoutGeometry ReadoutGeometryBase GeoPrimitives PixelReadoutGeometryLib
                  ActsEventCnvLib ActsToolInterfacesLib ActsCore ActsEventLib)

# Install files from the package:
atlas_install_python_modules( python/*.py POST_BUILD_CMD ${ATLAS_FLAKE8} )

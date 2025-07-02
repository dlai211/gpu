# Copyright (C) 2002-2024 CERN for the benefit of the ATLAS collaboration

# Declare the package name:
atlas_subdir( GeoModelConversion )

find_package(detray)

# Component(s) in the package:
atlas_add_component( GeoModelConversion
                     src/*.cxx
                     src/components/*.cxx
                     LINK_LIBRARIES AthenaBaseComps xAODTracking GaudiKernel AthContainers InDetPrepRawData InDetRawData
      InDetReadoutGeometry ReadoutGeometryBase GeoPrimitives PixelReadoutGeometryLib
                      ActsToolInterfacesLib ActsCore ActsEventLib ActsPluginGeoModel ActsPluginDetray
                     vecmem::core detray::io detray::core)


atlas_install_python_modules( python/*.py POST_BUILD_CMD ${ATLAS_FLAKE8} )

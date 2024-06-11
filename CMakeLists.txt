cmake_minimum_required(VERSION 3.12)

project(Anisotropy VERSION 1.0)

include(CMakeDependentOption)
include(cmake/environment.cmake)

find_package(Sofa.Core REQUIRED)
sofa_find_package(Sofa.Simulation.Core REQUIRED)
sofa_find_package(Sofa.Component.SolidMechanics.FEM.Elastic REQUIRED)
sofa_find_package(Sofa.Topology REQUIRED)
sofa_find_package(Sofa.GL REQUIRED)

set(Ani_SRC "src/Anisotropy")

set(HEADER_FILES
    ${Ani_SRC}/config.h.in

    ${Ani_SRC}/TetrahedronAnisotropicForceField.h
    ${Ani_SRC}/TetrahedronAnisotropicForceField.inl
    )

set(SOURCE_FILES
    ${Ani_SRC}/initAnisotropy.cpp

    ${Ani_SRC}/TetrahedronAnisotropicForceField.cpp
    )

set(README_FILES
    README.md
    )

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_ANISOTROPY")
target_link_libraries(${PROJECT_NAME}
    Sofa.Core
    Sofa.Component.SolidMechanics.FEM.Elastic
    Sofa.GL
)

sofa_create_package_with_targets(
    PACKAGE_NAME Anisotropy
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR "Anisotropy"
    RELOCATABLE "plugins"
    )

#include(cmake/packaging.cmake)

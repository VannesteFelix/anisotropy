cmake_minimum_required(VERSION 3.12)

project(Anisotropy VERSION 1.0)

include(cmake/environment.cmake)

find_package(Sofa.Config REQUIRED)
sofa_find_package(Sofa.Simulation.Core REQUIRED)
sofa_find_package(Sofa.Component.SolidMechanics.FEM.Elastic REQUIRED)
sofa_find_package(Sofa.Topology REQUIRED)
sofa_find_package(Sofa.GL REQUIRED)
sofa_find_package(SoftRobots.Inverse REQUIRED)

set(Ani_SRC "src/Anisotropy/component/")

set(HEADER_FILES
    ${Ani_SRC}/initAnisotropy.h.in

    ${Ani_SRC}/forcefield/TetrahedronAnisotropicForceField.h
    ${Ani_SRC}/forcefield/TetrahedronAnisotropicForceField.inl

    ${Ani_SRC}/constraint/AnisotropyControlPointActuator.h
    ${Ani_SRC}/constraint/AnisotropyControlPointActuator.inl    
    )

set(SOURCE_FILES
    ${Ani_SRC}/initAnisotropy.cpp

    ${Ani_SRC}/forcefield/TetrahedronAnisotropicForceField.cpp

    ${Ani_SRC}/constraint/AnisotropyControlPointActuator.cpp
    )

set(README_FILES
    README.md
    )

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${README_FILES})

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_ANISOTROPY")
target_link_libraries(${PROJECT_NAME}
    SoftRobots
    SoftRobots.Inverse
    Sofa.Core
    Sofa.Component.SolidMechanics.FEM.Elastic
    Sofa.GL
)

sofa_create_package_with_targets(
    PACKAGE_NAME Anisotropy
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    RELOCATABLE "plugins"
    )

#include(cmake/packaging.cmake)

cmake_minimum_required (VERSION 3.16)

set(PROJ_NAME signals-cpp)
project(${PROJ_NAME})

set(CMAKE_CXX_STANDARD 17)

option(BUILD_TESTS "Build Tests" ON)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(manif-geom-cpp REQUIRED)

add_library(${PROJ_NAME} INTERFACE)
target_include_directories(${PROJ_NAME} INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(${PROJ_NAME} INTERFACE Eigen3::Eigen manif-geom-cpp)

if (BUILD_TESTS)
    set(UNIT_TEST unit-tests)
    add_executable(${UNIT_TEST}
        tests/MainTest.cpp
        tests/SignalTest.cpp
        tests/IntegrationTest.cpp
        tests/DynamicsTest.cpp
    )
    target_link_libraries(${UNIT_TEST}
        ${PROJ_NAME}
        boost_unit_test_framework
    )
    add_test(NAME ${UNIT_TEST} COMMAND ${UNIT_TEST})

    add_custom_command(
        TARGET ${UNIT_TEST}
        COMMENT "Run unit tests"
        POST_BUILD
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMAND ${UNIT_TEST}
    )
endif()

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/${PROJ_NAME}ConfigVersion.cmake"
    VERSION 0.0
    COMPATIBILITY AnyNewerVersion
)

install(TARGETS ${PROJ_NAME}
    EXPORT ${PROJ_NAME}Targets
    LIBRARY DESTINATION lib COMPONENT Runtime
    ARCHIVE DESTINATION lib COMPONENT Development
    RUNTIME DESTINATION bin COMPONENT Runtime
    PUBLIC_HEADER DESTINATION include COMPONENT Development
    BUNDLE DESTINATION bin COMPONENT Runtime
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/${PROJ_NAME}Config.cmake.in"
    "${PROJECT_BINARY_DIR}/${PROJ_NAME}Config.cmake"
    INSTALL_DESTINATION lib/cmake/${PROJ_NAME}
)

install(EXPORT ${PROJ_NAME}Targets DESTINATION lib/cmake/${PROJ_NAME})
install(FILES "${PROJECT_BINARY_DIR}/${PROJ_NAME}ConfigVersion.cmake"
              "${PROJECT_BINARY_DIR}/${PROJ_NAME}Config.cmake"
        DESTINATION lib/cmake/${PROJ_NAME})
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

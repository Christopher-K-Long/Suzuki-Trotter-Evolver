cmake_minimum_required(VERSION 3.15...4.0)
project(Suzuki-Trotter-Evolver VERSION 1.1.0
                               DESCRIPTION
                                   "A C++ header-only library for evolving states under the Schrödinger equation using first-order Suzuki-Trotter and computing switching functions."
                               HOMEPAGE_URL "https://github.com/Christopher-K-Long/Suzuki-Trotter-Evolver"
                               LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
include(GNUInstallDirs)

option(Suzuki-Trotter-Evolver_INSTALL_LIBRARY
       "Enable installing of Suzuki-Trotter-Evolver library into default locations."
       ${PROJECT_IS_TOP_LEVEL})
option(Suzuki-Trotter-Evolver_BUILD_TESTING
       "Build and run Suzuki-Trotter-Evolver unit tests."
       ${PROJECT_IS_TOP_LEVEL})
option(Suzuki-Trotter-Evolver_BUILD_DOCS
       "Build the documentation for Suzuki-Trotter-Evolver."
       ${PROJECT_IS_TOP_LEVEL})

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIRS})

if(Suzuki-Trotter-Evolver_BUILD_TESTING)
  find_package(Catch2)
  if(NOT Catch2_FOUND)
    message(WARNING "Catch2 not found, not building tests")
    set(Suzuki-Trotter-Evolver_BUILD_TESTING OFF)
  endif()
endif()

add_library(Suzuki-Trotter-Evolver INTERFACE)
target_include_directories(
    Suzuki-Trotter-Evolver
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
              $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(Suzuki-Trotter-Evolver INTERFACE Eigen3::Eigen)
target_compile_features(Suzuki-Trotter-Evolver INTERFACE cxx_std_20)

if(Suzuki-Trotter-Evolver_BUILD_TESTING)
  enable_testing()

  include(Catch)
  add_subdirectory(tests)
  add_subdirectory(examples)
endif()

if(Suzuki-Trotter-Evolver_BUILD_DOCS)
  add_subdirectory(docs)
endif()

if(Suzuki-Trotter-Evolver_INSTALL_LIBRARY)
    install(TARGETS Suzuki-Trotter-Evolver
            EXPORT ${PROJECT_NAME}_Targets
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    include(CMakePackageConfigHelpers)
    write_basic_package_version_file("Suzuki-Trotter-Evolver_ConfigVersion.cmake"
                                    VERSION ${PROJECT_VERSION}
                                    COMPATIBILITY SameMajorVersion)
    configure_package_config_file(
        "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}_Config.cmake.in"
        "${PROJECT_BINARY_DIR}/${PROJECT_NAME}_Config.cmake"
        INSTALL_DESTINATION
        ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)
    install(EXPORT ${PROJECT_NAME}_Targets
            FILE ${PROJECT_NAME}Targets.cmake
            NAMESPACE ${PROJECT_NAME}::
            DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)
    install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}_Config.cmake"
                "${PROJECT_BINARY_DIR}/${PROJECT_NAME}_ConfigVersion.cmake"
            DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/Suzuki-Trotter-Evolver DESTINATION include)
endif()

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
include(CPack)
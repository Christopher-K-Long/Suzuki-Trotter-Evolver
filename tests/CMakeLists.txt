cmake_minimum_required(VERSION 3.23)

execute_process(COMMAND bash version_consistency.sh
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                RESULT_VARIABLE _EXIT_CODE)

if (NOT _EXIT_CODE EQUAL 0)
    message(FATAL_ERROR "The versions numbers are inconsistent.")
endif()

add_executable(Suzuki-Trotter-Evolver_unit_tests)
target_sources(Suzuki-Trotter-Evolver_unit_tests
               PRIVATE
               src/unit_tests/test_UnitaryEvolver.cpp)

target_include_directories(
    Suzuki-Trotter-Evolver_unit_tests PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                                              $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>)

target_link_libraries(Suzuki-Trotter-Evolver_unit_tests Catch2::Catch2WithMain)

catch_discover_tests(Suzuki-Trotter-Evolver_unit_tests)

list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")
add_custom_command(
     TARGET Suzuki-Trotter-Evolver_unit_tests
     COMMENT "Running tests..."
     POST_BUILD 
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
     COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIG> --rerun-failed --output-on-failure
)
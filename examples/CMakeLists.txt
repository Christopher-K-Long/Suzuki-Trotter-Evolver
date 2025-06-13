cmake_minimum_required(VERSION 3.23)

file(GLOB files "*.cpp")
foreach(file ${files})
  get_filename_component(example ${file} NAME)
  add_executable("${example}.exe")
  target_sources("${example}.exe"
                 PRIVATE
                 ${file})

  target_include_directories(
    "${example}.exe" PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                             $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
  )

  add_custom_command(
    TARGET "${example}.exe"
    COMMENT "Checking example ${example} runs..."
    POST_BUILD 
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/examples
    COMMAND "${example}.exe"
  )
endforeach()
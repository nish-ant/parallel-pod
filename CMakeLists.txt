cmake_minimum_required(VERSION 3.0)
project(POD)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

set(CMAKE_CXX_FLAGS_RELEASE "-std=c++1z -O3")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++1z -g")

#file(GLOB_RECURSE POD_SRC "src/*.cpp")

set(POD_DEPENDENCIES_DIR ${CMAKE_SOURCE_DIR}/dependencies)
set(POD_INCLUDE_DIRS ${POD_DEPENDENCIES_DIR}/eigen ${POD_DEPENDENCIES_DIR}/ezOptionParser src)
include_directories(${POD_INCLUDE_DIRS})

if (CMAKE_BUILD_TYPE MATCHES Debug)
    message(" ")
    message("CMAKE IN DEBUG MODE")
    message(" ")
elseif(CMAKE_BUILD_TYPE MATCHES Release)
    message(" ")
    message("CMAKE IN RELEASE MODE")
    message(" ")
endif ()

set(UTILS_SRC "src/utils.cpp")
add_library(UTILS STATIC ${UTILS_SRC})

set(POD_SRC "src/pod.cpp")
add_executable(POD ${POD_SRC})
target_link_libraries(POD UTILS)
set_target_properties(POD PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set(REC_SRC "src/rec.cpp")
add_executable(REC ${REC_SRC})
target_link_libraries(REC UTILS)
set_target_properties(REC PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

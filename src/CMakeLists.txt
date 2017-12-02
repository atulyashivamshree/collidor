cmake_minimum_required(VERSION 3.0)
project(collidor)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfatal-errors -Wformat=5")

include_directories("/usr/include/eigen3")
set(FCL_PATH "/home/atulya/Documents/eecs587/fcl")
find_package( octomap REQUIRED)
#link_directories("/usr/local/lib")

#include_directories("${FCL_PATH}/include" "${FCL_PATH}/build/include" "${FCL_PATH}/build/test")
include_directories("${FCL_PATH}/include" "${FCL_PATH}/build/include" "${OCTOMAP_INCLUDE_DIRS}")
link_directories("${FCL_PATH}/build/lib")

add_library(test_fcl_utility_local test_fcl_utility.cpp)
target_link_libraries(test_fcl_utility_local fcl)

add_executable(load_files parser.cpp)
target_link_libraries(load_files fcl test_fcl_utility_local ccd "${OCTOMAP_LIBS}")


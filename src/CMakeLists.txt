cmake_minimum_required(VERSION 3.0)
project(collidor)
set(CMAKE_CXX_STANDARD 11)

include_directories("/usr/include/eigen3")
set(FCL_PATH "/home/atulya/Documents/eecs587/fcl")

include_directories("${FCL_PATH}/include" "${FCL_PATH}/build/include")
link_directories("${FCL_PATH}/build/lib")

add_executable(load_files parser.cpp test_fcl_utility.cpp)
target_link_libraries(load_files fcl octomap ccd)


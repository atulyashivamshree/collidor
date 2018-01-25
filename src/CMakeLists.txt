cmake_minimum_required(VERSION 3.0)
project(collidor)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wfatal-errors -Wformat=5")
enable_testing()

include_directories("/usr/include/eigen3")
set(FCL_PATH "/home/atulya/eecs587/fcl")

include_directories("${FCL_PATH}/include" "${FCL_PATH}/build/include" )
link_directories("${FCL_PATH}/build/lib")

add_executable(save_to_bvh utils/fcl_cdp_conversion.cpp)
target_link_libraries(save_to_bvh fcl ccd)

add_executable(compute_dist_bvh utils/compute_FCL_distance.cpp)
target_link_libraries(compute_dist_bvh fcl ccd)

add_executable(verify_dist utils/distance_verification.cpp)
target_link_libraries(verify_dist fcl ccd)

add_executable(test_vector test/test_vector_math.cpp)
add_executable(test_rss test/test_rss.cpp)

add_test(vector test_vector)
add_test(RSS test_rss)

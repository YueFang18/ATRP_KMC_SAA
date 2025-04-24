cmake_minimum_required(VERSION 3.19)
project(effPropaRecorder)

set(CMAKE_CXX_STANDARD 11)

find_package(Boost 1.76.0 COMPONENTS system filesystem REQUIRED)
# 如果需要boost的thread同时添加下面一行
find_package(Boost COMPONENTS thread REQUIRED)
if(NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find boost!")
endif()

include_directories(
        include
        ${Boost_INCLUDE_DIRS})

add_executable(effATRP_MMA src/eff_propa_reorder.cpp include/add_schedule.h include/my_rand.h)

target_link_libraries(effATRP_MMA ${Boost_LIBRARIES})
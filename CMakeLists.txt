cmake_minimum_required (VERSION 3.14)
project (GESMO)

set (CMAKE_CXX_STANDARD 17)

add_library(problems STATIC problems_gecco23.cpp)
target_include_directories(problems PUBLIC 
    include 
    external/fmt/include
    external/cxxopts/include 
    external/json/include 
)

add_executable(gsemo gsemo.cpp)
target_link_libraries(gsemo problems)

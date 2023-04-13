cmake_minimum_required (VERSION 3.14)
project (GESMO)

set (CMAKE_CXX_STANDARD 17)

add_library(problems STATIC src/problems.cpp)
target_include_directories(problems PUBLIC 
    include
    src 
    external/fmt/include
    external/cxxopts/include 
    external/json/include
    /usr/local/include
)

add_executable(gsemo src/main.cpp)
add_executable(agsemo src/agsemo.cpp)
target_link_libraries(gsemo problems)
target_link_libraries(agsemo problems)
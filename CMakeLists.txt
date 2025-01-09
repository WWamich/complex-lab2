cmake_minimum_required(VERSION 3.26)
project(calculusgit)

find_package(Python3 3.12 COMPONENTS Interpreter Development NumPy REQUIRED)
include_directories(${Python3_INCLUDE_DIRS} include)

include_directories("C:/Users/admir/CLionProjects/calculusgit/libs/matplotlib-cpp-master/")
add_executable(calculusgit main.cpp
)
target_link_libraries(calculusgit Python3::Python)



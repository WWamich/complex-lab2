cmake_minimum_required(VERSION 3.26)
project(complex-lab2)

find_package(Python3 3.12 COMPONENTS Interpreter Development NumPy REQUIRED)
include_directories(${Python3_INCLUDE_DIRS} include)

include_directories("C:/Users/admir/CLionProjects/complex-lab2/libs/matplotlib-cpp-master/")
add_executable(complex-lab2 main.cpp
)
target_link_libraries(complex-lab2 Python3::Python)



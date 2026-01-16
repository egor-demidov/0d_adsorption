cmake_minimum_required(VERSION 4.1)
project(0d_adsorption)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SUNDIALS CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(Ceres CONFIG REQUIRED)

add_executable(0d_adsorption main.cpp
        model.cpp
        model.h)

#set(CMAKE_CXX_FLAGS "-O3 -flto=auto ${CMAKE_CXX_FLAGS}")

target_link_libraries(0d_adsorption PRIVATE
        fmt::fmt
        Ceres::ceres
        SUNDIALS::cvode
        SUNDIALS::nvecserial
        SUNDIALS::sunmatrixdense
        SUNDIALS::sunlinsoldense
)

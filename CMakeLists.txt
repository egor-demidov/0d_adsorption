cmake_minimum_required(VERSION 4.1)
project(0d_adsorption)

set(CMAKE_CXX_STANDARD 20)

find_package(SUNDIALS CONFIG REQUIRED)

add_executable(0d_adsorption main.cpp)

target_link_libraries(0d_adsorption PRIVATE
        sundials_cvode
        sundials_nvecserial
        sundials_sunlinsoldense
        sundials_sunmatrixdense
)

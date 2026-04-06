cmake_minimum_required(VERSION 3.28)
project(0d_adsorption)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SUNDIALS CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(Ceres CONFIG REQUIRED)
find_package(nlohmann_json CONFIG REQUIRED)
find_package(RapidJSON CONFIG REQUIRED)

add_compile_definitions(_USE_MATH_DEFINES)

add_executable(0d_adsorption_fit_chained fit_chained.cpp
        model_chained.cpp)

add_executable(0d_adsorption_chained run_chained.cpp
        model_chained.cpp)

#set(CMAKE_CXX_FLAGS "-O3 -flto=auto ${CMAKE_CXX_FLAGS}")

if (${EMSCRIPTEN})
    set(CMAKE_EXE_LINKER_FLAGS "-s EXPORT_ES6=1 -s MODULARIZE=1 -sDISABLE_EXCEPTION_CATCHING=0 -sEXPORTED_RUNTIME_METHODS=getExceptionMessage -fexceptions ${CMAKE_EXE_LINKER_FLAGS}")
    set(EMSCRIPTEN_LINK_FLAGS "--bind")
    add_compile_definitions(EMSCRIPTEN)
endif ()

target_link_libraries(0d_adsorption_fit_chained PRIVATE
        fmt::fmt
        Ceres::ceres
        SUNDIALS::cvode
        SUNDIALS::nvecserial
        SUNDIALS::sunmatrixdense
        SUNDIALS::sunlinsoldense
        nlohmann_json::nlohmann_json
        rapidjson
        ${EMSCRIPTEN_LINK_FLAGS}
)

target_link_libraries(0d_adsorption_chained PRIVATE
        fmt::fmt
        SUNDIALS::cvode
        SUNDIALS::nvecserial
        SUNDIALS::sunmatrixdense
        SUNDIALS::sunlinsoldense
        nlohmann_json::nlohmann_json
)

cmake_minimum_required(VERSION 3.2)
project(DosageConvertor VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 11)
execute_process(COMMAND date OUTPUT_VARIABLE DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND whoami OUTPUT_VARIABLE USER OUTPUT_STRIP_TRAILING_WHITESPACE)

add_definitions(-DVERSION="${PROJECT_VERSION}" -DUSER="${USER}" -DDATE="${DATE}")

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()



find_package(ZLIB REQUIRED)
find_library(STATGEN_LIBRARY StatGen)
add_executable(DosageConvertor
        src/HaplotypeSet.cpp
        src/HaplotypeSet.h
        src/Main.cpp)

target_link_libraries(DosageConvertor ${STATGEN_LIBRARY} ${ZLIB_LIBRARIES})

install(TARGETS DosageConvertor RUNTIME DESTINATION bin)

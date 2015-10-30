cmake_minimum_required(VERSION 2.8)

####
# project settings
####

set(sample_VERSION_MAJOR 1)
set(sample_VERSION_MINOR 0)


project(lensingCalculator)

include_directories(include)

####
# sources
####

#set(SOURCES lensingCalculator.cpp myUtilities.cpp inputFunctions.cpp lensFitter.cpp pixelmapFunctions.cpp)
file(GLOB SOURCES "src/*.cpp")

####
# dependencies
####

find_package(GLAMER NO_MODULE REQUIRED)
include_directories(${GLAMER_INCLUDE_DIRS})

####
# omp
####
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#add_compile_options( "-mcmodel=medium" )
####
# build
####

add_executable(lensing_calculator ${SOURCES})

target_link_libraries(lensing_calculator ${GLAMER_LIBRARIES})
cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0048 NEW) # avoid annoing warnings
project(rof LANGUAGES CXX)

message(STATUS "Configuring project: " ${PROJECT_NAME})

###############################################################################

# NOTE: -fPIC is added so that the static library can be included into
#       shared libraries.
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -w -fPIC -O3 ")
set(CMAKE_BUILD_TYPE Release)

###############################################################################

## # Enable C99
## if (CMAKE_VERSION VERSION_LESS "3.1")
## 	set (CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c99")
## else ()
## 	set (CMAKE_CXX_STANDARD 99)
## endif ()

###############################################################################

find_package( OpenMP REQUIRED)
if(OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

###############################################################################

include_directories(${IIO_INCLUDE_DIR})

# rof executable
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
add_executable(rof main.cpp)
target_link_libraries(rof iio m)

###############################################################################

# Expose public includes to other subprojects through cache variable
set(RDPOF_INCLUDE_DIR ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)


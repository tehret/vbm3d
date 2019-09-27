cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0048 NEW) # avoid annoing warnings
project(iio LANGUAGES C)

message(STATUS "Configuring project: " ${PROJECT_NAME})

###############################################################################

# NOTE: -fPIC is added so that the static library can be included into
#       shared libraries.
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -w -fPIC")
set(CMAKE_BUILD_TYPE Release)

###############################################################################

# Enable C99
if (CMAKE_VERSION VERSION_LESS "3.1")
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
else ()
	set (CMAKE_C_STANDARD 99)
endif ()

###############################################################################

find_package(TIFF REQUIRED)
find_package(JPEG REQUIRED)
find_package(PNG  REQUIRED)

#include_directories(PUBLIC ${TIFF_INCLUDE_DIR}
#                    PUBLIC ${JPEG_INCLUDE_DIR}
#                    PUBLIC ${PNG_INCLUDE_DIRS}
#                    PUBLIC ${FFTW_INCLUDES})

add_library(iio STATIC iio.c)

target_link_libraries(iio png tiff jpeg)

###############################################################################

# Expose public includes to other subprojects through cache variable
set(IIO_INCLUDE_DIR ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "${PROJECT_NAME}: Include Directories" FORCE)


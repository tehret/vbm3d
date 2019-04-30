set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

add_library(iio_tvl1 SHARED iio.c)

add_executable(tvl1flow main.c)

link_directories(${CMAKE_BINARY_DIR}/lib)
target_link_libraries(tvl1flow iio_tvl1 tiff jpeg png m)

configure_file(tvl1flow.sh ${CMAKE_BINARY_DIR}/bin/tvl1flow.sh COPYONLY)

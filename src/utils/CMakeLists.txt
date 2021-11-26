set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/../Utilities)

add_library(utilities SHARED
	../Utilities/iio.c
	../Utilities/LibImages.cpp
	../Utilities/LibVideoT.cpp
	../Utilities/mt19937ar.c
	../Utilities/Utilities.cpp
	)


set_target_properties(utilities PROPERTIES
	SOVERSION ${PROJECT_VERSION})

install(TARGETS utilities
    LIBRARY DESTINATION lib
    COMPONENT libs)

add_executable(psnr compute_psnr.cpp)
add_executable(addnoise addnoise.cpp)

configure_file(VBM3Ddenoising_OF.sh ${CMAKE_BINARY_DIR}/bin/VBM3Ddenoising_OF.sh COPYONLY)
configure_file(VBM3Ddenoising_OF_multiscale.sh ${CMAKE_BINARY_DIR}/bin/VBM3Ddenoising_OF_multiscale.sh COPYONLY)
configure_file(VBM3Ddenoising_multiscale.sh ${CMAKE_BINARY_DIR}/bin/VBM3Ddenoising_multiscale.sh COPYONLY)

link_directories(${CMAKE_BINARY_DIR}/lib)
target_link_libraries(psnr      utilities png tiff jpeg)
target_link_libraries(addnoise  utilities png tiff jpeg)

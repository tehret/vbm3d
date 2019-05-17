set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/Utilities)

add_executable(VBM3Ddenoising main.cpp)

include_directories(${FFTW_INCLUDEs})


add_library(VBM3DLibrary
	vbm3d.cpp
    lib_transforms.cpp
	Utilities/iio.c
	Utilities/LibImages.cpp
	Utilities/mt19937ar.c
    Utilities/LibVideoT.cpp
    Utilities/Utilities.cpp
	)


set_target_properties(VBM3DLibrary PROPERTIES
	SOVERSION ${PROJECT_VERSION})

install(TARGETS VBM3DLibrary
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    COMPONENT libs)


link_directories(${CMAKE_BINARY_DIR}/lib)
target_link_libraries(VBM3Ddenoising VBM3DLibrary png ${FFTW_LIBRARIES})

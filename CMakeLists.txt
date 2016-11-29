cmake_minimum_required(VERSION 2.6)
project(mnf C CXX)

#BLAS library
if (USE_MKL_LIBRARIES)
	set(CMAKE_BLAS_LIBS_INIT mkl_intel_lp64 mkl_intel_thread mkl_core iomp5)
else()
	find_package(BLAS REQUIRED)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#mnf library
add_library(mnf SHARED src/mnf_c.cpp src/mnf_linebyline.cpp src/readimage.cpp)

#mnf executable
add_executable(mnf-bin src/main.cpp src/readimage.cpp)
set_target_properties(mnf-bin PROPERTIES OUTPUT_NAME mnf) #set mnf-bin binary name to mnf
target_link_libraries(mnf-bin mnf ${CMAKE_BLAS_LIBS_INIT} lapacke)

#install libraries, binaries and header files
install(TARGETS mnf DESTINATION lib)
install(TARGETS mnf-bin DESTINATION bin)
install(FILES include/mnf_linebyline.h include/mnf.h DESTINATION include)

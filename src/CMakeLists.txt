###########################################################
# PBD algorithm implementation                            #
# Thesis: Muscle Fibres Deformation using Particle System #
# Author: Martin Cervenka                                 #
# Version: 5/2019                                         #
###########################################################

cmake_minimum_required(VERSION 2.8)
project(simulation)

find_package(VTK 8.2.0 EXACT REQUIRED)
include(${VTK_USE_FILE})
# using GCC (simple and straigthforward ;))
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "-std=c++11 -pedantic-errors -Wall -Wextra -Werror -Wno-parentheses")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-exceptions -fopenmp -Ofast -lm")
# using visual c++ (notice how hard is it to setup and optimize this
# piece of crap (MSVC) compared with GCC ;/)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	SET(CMAKE_CXX_FLAGS_RELEASE
		"${CMAKE_CXX_FLAGS_RELEASE} /arch:AVX /Gy /Oi /Oy /Ox /GS-" CACHE STRING
		"Flags used by the compiler during the release builds." FORCE)
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /ZI /MP" CACHE STRING
		"Flags used by the compiler during the debug builds." FORCE)
	set(CMAKE_EXE_LINKER_FLAGS_RELEASE
		"${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO}" CACHE STRING
		"Flags used by the linker during the release builds." FORCE)
	set(CMAKE_SHARED_LINKER_FLAGS_RELEASE
		"${CMAKE_SHARED_LINKER_FLAGS_RELWITHDEBINFO}" CACHE STRING
		"Flags used by the linker during the release builds." FORCE)
endif()

add_executable(pbd
	main.cxx
	pbd.cxx
	graph.cxx
	point_3D.cxx
	data.cxx
	matrix_4x4.cxx
	collision/collision.cxx
	collision/voxel.cxx
)
add_executable(generator
	generator.cxx
)
target_link_libraries(pbd ${VTK_LIBRARIES})
target_link_libraries(generator ${VTK_LIBRARIES})
#target_compile_definitions(pbd PRIVATE NO_HUD)
#vtk_module_autoinit(TARGETS pbd MODULES ${VTK_LIBRARIES})

#-----------------------------------------------------------------------
# - Enforce an out-of-source builds before anything else
#
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "CRADLE++ requires an out-of-source build.")
    message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
    message(STATUS "CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there")
    message(FATAL_ERROR "in-source build detected")
endif()

cmake_minimum_required(VERSION 3.1)

project(CRADLE++ VERSION 1.0 LANGUAGES CXX)

set(CMAKE_MACOSX_RPATH OFF)

# # use, i.e. don't skip the full RPATH for the build tree
# set(CMAKE_SKIP_BUILD_RPATH FALSE)
#
# # when building, don't use the install RPATH already
# # (but later on when installing)
# set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
#
# set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#
# # add the automatically determined parts of the RPATH
# # which point to directories outside the build tree to the install RPATH
# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

# - Prepend our own CMake Modules to the search path
# NB: if our custom modules include others that we don't supply, those in
# the base path will be used, so watch for incompatibilities!!
#
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})
include(CMakeDependentOption)
include(GNUInstallDirs)
include(CMakeUninstallTarget)

#---------------------------------------------------------------------------------------
# compiler config
#---------------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    #add_compile_options("-Wconversion")
    add_compile_options("-pedantic")
    add_compile_options("-Wfatal-errors")
    add_compile_options("-pthread")
    add_compile_options("-Wno-reorder")
endif()

# Offer the user the choice of overriding the installation directories
# set(INSTALL_LIB_DIR ${PROJECT_BINARY_DIR}/lib CACHE PATH "Installation directory for libraries")
# set(INSTALL_BIN_DIR ${PROJECT_BINARY_DIR}/bin CACHE PATH "Installation directory for executables")
# set(INSTALL_INCLUDE_DIR ${PROJECT_BINARY_DIR}/include CACHE PATH
#   "Installation directory for header files")

add_library(Cradle SHARED src/ConfigParser.cc src/DecayChannel.cc src/DecayManager.cc src/DecayMode.cc src/Particle.cc src/SpectrumGenerator.cc src/ThreadPool.cc)
add_executable(CRADLE++ src/CRADLE++.cc)

find_package(Boost REQUIRED)
find_package(GSL REQUIRED)
find_package(Threads)
find_package(ROOT REQUIRED)

target_include_directories(Cradle PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(Cradle PUBLIC ${GSL_INCLUDE_DIRS})
target_include_directories(Cradle PUBLIC ${ROOT_INCLUDE_DIRS})

target_include_directories(
    Cradle
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
set(HEADER_BASE "${CMAKE_CURRENT_SOURCE_DIR}/include")


target_link_libraries(Cradle PUBLIC "-L/usr/local/lib -lgsl -lgslcblas")
target_link_libraries(Cradle PUBLIC ${Boost_LIBRARIES})
target_link_libraries(Cradle PUBLIC ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(Cradle PUBLIC ${ROOT_LIBRARIES})
target_link_libraries(CRADLE++ PRIVATE Cradle)


# add_custom_command(TARGET CRADLE++
#                    POST_BUILD
# 		   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:CRADLE++> ${PROJECT_BINARY_DIR}/bin/$<TARGET_FILE_NAME:CRADLE++>)
#
# install(TARGETS CRADLE++ DESTINATION bin)

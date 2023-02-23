cmake_minimum_required(VERSION 3.10)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build")

project(mtqr LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17 REQUIRED)

set(INCLUDES "include")
set(SOURCES "src")
set(TESTS "tests")
set(DATA "data")

include_directories(${INCLUDES})

find_package(GSL REQUIRED)
if (${GSL_FOUND})
	include_directories(${GSL_INCLUDE_DIR})
else ()
	message("Error: failed to locate GSL!")
endif()

find_package(Boost REQUIRED)
if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
else()
    message("Error: failed to locate Boost!")
endif()

add_library(mtqr STATIC ${SOURCES}/mtqr.cpp ${SOURCES}/data_management.cpp ${SOURCES}/monomial_transformation.cpp ${SOURCES}/vector_operations.cpp ${DATA}/BETAS.h ${DATA}/GL_NODES.h ${DATA}/GL_WEIGHTS.h)
set_target_properties(mtqr
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build"
)
target_include_directories(mtqr PUBLIC ${INCLUDES} ${DATA})

if (WIN32)
    set(DEPENDENCIES GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
else(UNIX)
	set(DEPENDENCIES m quadmath GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
endif ()

add_executable(drivers ${TESTS}/drivers/main.cpp)
include_directories(drivers ${INCLUDES})
target_link_libraries(drivers mtqr ${DEPENDENCIES})

add_executable(bessels ${TESTS}/bessels/main.cpp)
include_directories(bessels ${INCLUDES})
target_link_libraries(bessels mtqr ${DEPENDENCIES})
cmake_minimum_required(VERSION 3.10)

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
target_include_directories(mtqr PUBLIC ${INCLUDES} ${DATA})

if (WIN32)
    set(DEPENDENCIES GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
else(UNIX)
	set(DEPENDENCIES m quadmath GSL::gsl GSL::gslcblas ${Boost_LIBRARIES})
endif ()

add_executable(drivers ${TESTS}/drivers/main.cpp)
include_directories(drivers ${INCLUDES})
target_link_libraries(drivers mtqr ${DEPENDENCIES})

add_executable(monomial ${TESTS}/monomial/main.cpp)
include_directories(monomial ${INCLUDES})
target_link_libraries(monomial mtqr ${DEPENDENCIES})

add_executable(integration ${TESTS}/integration/main.cpp)
include_directories(integration ${INCLUDES})
target_link_libraries(integration mtqr ${DEPENDENCIES})

add_executable(bessel_1 ${TESTS}/bessel/bessel_1/main.cpp)
include_directories(bessel_1 ${INCLUDES})
target_link_libraries(bessel_1 mtqr ${DEPENDENCIES})

add_executable(bessel_2 ${TESTS}/bessel/bessel_2/main.cpp)
include_directories(bessel_2 ${INCLUDES})
target_link_libraries(bessel_2 mtqr ${DEPENDENCIES})

add_executable(bessel_3 ${TESTS}/bessel/bessel_3/main.cpp)
include_directories(bessel_3 ${INCLUDES})
target_link_libraries(bessel_3 mtqr ${DEPENDENCIES})
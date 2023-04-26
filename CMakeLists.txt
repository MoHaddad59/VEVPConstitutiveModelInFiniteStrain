cmake_minimum_required(VERSION 3.16)
project(IncSecantMFHFiniteStrain)

add_definitions(-DNPY_NO_DEPRECATED_API=NPY_1_22_API_VERSION)

include_directories(/usr/local/lib/python3.8/dist-packages/numpy/core/include/numpy)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)

set(PYTHON_MODULE_DIR ${../CMAKE_CURRENT_SOURCE_DIR})
set(PYTHONPATH "${PYTHON_MODULE_DIR}:${Python3_SITELIB}")
include_directories(${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS})

add_executable(ExecutableIncrSecantFiniteStrainMethod Testing.cpp VEVPFiniteStrainNS.cpp)

target_link_libraries(ExecutableIncrSecantFiniteStrainMethod  Python3::Python )

set_property(TARGET ExecutableIncrSecantFiniteStrainMethod PROPERTY ENVIRONMENT "PYTHONPATH=${PYTHONPATH}")

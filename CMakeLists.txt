cmake_minimum_required(VERSION 3.12)
project(ShallowWaterSimulation)

set(CMAKE_CXX_STANDARD 11)

# Set compiler flags
add_compile_options(-Wall -O3 -flto -g)

# Include directories
include_directories(
    include
    /opt/homebrew/Cellar/libomp/16.0.6/include
    /opt/homebrew/Cellar/openblas/0.3.23/include
)

# List of source files
file(GLOB_RECURSE SRCS "src/*.cpp")


# Find Boost
find_package(Boost REQUIRED COMPONENTS program_options)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
endif()


# Add executable
add_executable(main ${SRCS})

# Link libraries
target_link_libraries(main 
    blas 
    Boost::program_options
    # Uncomment the following line if OpenMP is needed
    # gomp
)

# Set the output directory for the executable
set_target_properties(main PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../output)


# Custom target for tests
function(add_custom_test TEST_NAME TEST_PARAMS)
    add_custom_target(${TEST_NAME}
        COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=1 $<TARGET_FILE:main> ${TEST_PARAMS}
        COMMAND ${CMAKE_COMMAND} -E echo "Running Gnuplot..."
        COMMAND gnuplot ${CMAKE_SOURCE_DIR}/plot/3d.gnu
        #COMMAND open ${CMAKE_SOURCE_DIR}/output/3d.png
        DEPENDS main
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/..
        COMMENT "Running ${TEST_NAME} Simulation and Gnuplot"
    )
    # Add the test cases for ctest
    add_test(NAME ${TEST_NAME} COMMAND $<TARGET_FILE:main> ${TEST_PARAMS})
endfunction()

# Define test parameters
set(TEST1_PARAMS --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 1 --choice 1)
set(TEST2_PARAMS --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 2 --choice 1)
set(TEST3_PARAMS --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 3 --choice 1)
set(TEST4_PARAMS --dt 0.1 --T 80 --Nx 100 --Ny 100 --ic 4 --choice 1)

# Add custom tests
add_custom_test(test1 "${TEST1_PARAMS}")
add_custom_test(test2 "${TEST2_PARAMS}")
add_custom_test(test3 "${TEST3_PARAMS}")
add_custom_test(test4 "${TEST4_PARAMS}")

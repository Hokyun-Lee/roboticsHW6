cmake_minimum_required(VERSION 3.0)

project(panda_control_real CXX)

add_definitions(-DDOF=7)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Franka REQUIRED)
find_package(Eigen3 REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

include_directories(include
${EIGEN3_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} 
    src/main.cpp
    src/controller.cpp
    src/motion_generator.cpp)

target_link_libraries(${PROJECT_NAME}
Franka::Franka pthread
)


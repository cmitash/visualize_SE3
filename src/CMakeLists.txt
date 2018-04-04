cmake_minimum_required(VERSION 2.8.3)
project(visualize_se3)
set(SRC_DIR ${PROJECT_SOURCE_DIR}/)

set(PCL_DIR /home/chaitanya/software/pcl-pcl-1.7.2/build)
find_package (PCL REQUIRED COMPONENTS io common visualization features) 

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  sensor_msgs 
  cv_bridge
  pcl_ros
  tf
  image_transport 
  image_geometry
  geometry_msgs
)
find_package(OpenCV REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_BUILD_TYPE Debug)

###################################
## catkin specific configuration ##
###################################
catkin_package(
  DEPENDS cv_bridge
  CATKIN_DEPENDS 
  image_transport
  image_geometry 
  roscpp
  pcl_ros
  tf
  geometry_msgs
)

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(${SRC_DIR})
include_directories(
  ${catkin_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
  ${PCL_INCLUDE_DIRS}
  ${CATKIN_DEVEL_PREFIX}/include
)


## Declare a C++ executable
add_executable(visualize_se3 ${SRC_DIR}/visualize_se3.cpp)
add_dependencies(visualize_se3 ${catkin_EXPORTED_TARGETS})

## Specify libraries to link a library or executable target against
target_link_libraries(visualize_se3
  ${catkin_LIBRARIES}
  ${OpenCV_LIBRARIES}
  ${PCL_LIBRARIES}
)
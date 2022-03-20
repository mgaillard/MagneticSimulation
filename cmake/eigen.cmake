if(TARGET dlib::dlib)
    return()
endif()

message(STATUS "Third-party (external): creating target 'eigen::eigen'")

# Source: https://stackoverflow.com/questions/65860094/how-to-add-eigen-library-to-a-cmake-c-project-via-fetchcontent

# Retrieve and activate Eigen3
include(FetchContent)
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW TRUE
)

set(EIGEN_BUILD_DOC OFF)
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_PKGCONFIG OFF)
set( OFF)

FetchContent_MakeAvailable(Eigen)

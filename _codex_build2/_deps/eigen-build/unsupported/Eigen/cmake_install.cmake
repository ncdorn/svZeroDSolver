# Install script for directory: /Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/AdolcForward"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/AlignedVector3"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/ArpackSupport"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/AutoDiff"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/BVH"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/EulerAngles"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/FFT"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/IterativeSolvers"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/KroneckerProduct"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/LevenbergMarquardt"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/MatrixFunctions"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/MPRealSupport"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/NNLS"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/NonLinearOptimization"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/NumericalDiff"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/OpenGLSupport"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/Polynomials"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/SparseExtra"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/SpecialFunctions"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/Splines"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/Tensor"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/TensorSymmetry"
    "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/ThreadPool"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/build/_deps/eigen-src/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/svz/repos/svZeroDSolver/_codex_build2/_deps/eigen-build/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()


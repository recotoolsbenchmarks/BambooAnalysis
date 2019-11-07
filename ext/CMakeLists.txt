cmake_minimum_required( VERSION 2.8 )

project(BambooExt LANGUAGES CXX)

include(GNUInstallDirs)
set(CMAKE_INSTALL_LIBDIR "lib")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

if(EXISTS "$ENV{CMAKE_PREFIX_PATH}/include/boost")
    message(STATUS "Will use $ENV{CMAKE_PREFIX_PATH} as base path for boost")
    set(BOOST_ROOT $ENV{CMAKE_PREFIX_PATH})
    set(Boost_NO_BOOST_CMAKE ON)
endif()

find_package(ROOT REQUIRED COMPONENTS MathCore Hist ROOTVecOps)
find_package(Boost REQUIRED)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(BambooRandom SHARED src/bamboorandom.cc)
target_include_directories(BambooRandom PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(BambooRandom ROOT::MathCore)
install(TARGETS BambooRandom LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_library(BinnedValues SHARED src/BinnedValues.cc src/BinnedValuesJSONParser.cc src/BTagCalibrationStandalone.cpp)
target_compile_definitions(BinnedValues PRIVATE NO_CHECK_ALL_FORMULAS)
target_include_directories(BinnedValues PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(BinnedValues ROOT::Hist Boost::boost)
install(TARGETS BinnedValues LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_library(BambooLumiMask SHARED src/LumiMask.cc)
target_include_directories(BambooLumiMask PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(BambooLumiMask Boost::boost)
install(TARGETS BambooLumiMask LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/jetclasses_filenames.txt FILENAMES_CMSJET_INC REGEX ".*/interface/.*")
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/jetclasses_filenames.txt FILENAMES_CMSJET_SRC REGEX ".*/src/.*")
string(REGEX REPLACE "[^;]+/interface/([^;]+)" "${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include/\\1" INC_CMSJET "${FILENAMES_CMSJET_INC}")
string(REGEX REPLACE "[^;]+/src/([^;]+)" "${CMAKE_CURRENT_BINARY_DIR}/CMSJet/src/\\1" SRC_CMSJET "${FILENAMES_CMSJET_SRC}")
add_custom_target(CMSJetSources COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/getjetclasses.sh COMMENT "Copying jet classes from CMSSW" BYPRODUCTS ${INC_CMSJET} ${SRC_CMSJET})
install(FILES ${INC_CMSJET} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(JMEObjects SHARED ${SRC_CMSJET} src/JMESystematicsCalculators.cc)
add_dependencies(JMEObjects CMSJetSources)
target_include_directories(JMEObjects PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/CMSJet/src> $<INSTALL_INTERFACE:include>)
target_link_libraries(JMEObjects ROOT::Hist ROOT::ROOTVecOps Boost::boost BambooRandom)
install(TARGETS JMEObjects LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_library(RoccoR SHARED src/RoccoR.cc src/RochesterCorrectionCalculator.cc)
target_include_directories(RoccoR PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(RoccoR ROOT::ROOTVecOps Boost::boost BambooRandom)
install(TARGETS RoccoR LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

## ML extensions

find_package(TensorflowC)
if(TensorflowC_FOUND)
  add_library(BambooTensorflowC SHARED src/bambootensorflowc.cc)
  target_compile_features(BambooTensorflowC PRIVATE cxx_std_17)
  target_include_directories(BambooTensorflowC
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
  target_link_libraries(BambooTensorflowC TensorflowC)
  install(TARGETS BambooTensorflowC
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  message(STATUS "Building Tensorflow C API integration (from ${TENSORFLOWC_INCLUDE_DIR})")
else()
  message(STATUS "Tensorflow C API integration will not be built. Specify CMAKE_PREFIX_PATH with the correct prefix to include it")
endif()

find_package(Torch)
if(Torch_FOUND)
  add_library(BambooTorch SHARED src/bambootorch.cc)
  target_compile_features(BambooTorch PRIVATE cxx_std_17)
  target_include_directories(BambooTorch
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
  target_link_libraries(BambooTorch torch)
  install(TARGETS BambooTorch
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  message(STATUS "Building libTorch integration with libTorch from ${TORCH_INSTALL_PREFIX}")
else()
  message(STATUS "libTorch integration will not be built. Specify CMAKE_PREFIX_PATH with the libTorch prefix to include it")
endif()

find_package(lwtnn)
if(lwtnn_FOUND)
  add_library(BambooLwtnn SHARED src/bamboolwtnn.cc)
  target_compile_features(BambooLwtnn PRIVATE cxx_std_17)
  target_include_directories(BambooLwtnn
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
  target_link_libraries(BambooLwtnn lwtnn::lwtnn)
  install(TARGETS BambooLwtnn
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
  message(STATUS "Building lwtnn integration with lwtnn from ${lwtnn_INSTALL_DIR}")
else()
  message(STATUS "lwtnn integration will not be built")
endif()

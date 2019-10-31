cmake_minimum_required( VERSION 2.8 )

project(BambooExt LANGUAGES CXX)

include(GNUInstallDirs)
set(CMAKE_INSTALL_LIBDIR "lib")

if(EXISTS "$ENV{CMAKE_PREFIX_PATH}/include/boost")
    message(STATUS "Will use $ENV{CMAKE_PREFIX_PATH} as base path for boost")
    set(BOOST_ROOT $ENV{CMAKE_PREFIX_PATH})
    set(Boost_NO_BOOST_CMAKE ON)
endif()

find_package(ROOT REQUIRED COMPONENTS Hist ROOTVecOps)
find_package(Boost REQUIRED)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(BinnedValues SHARED src/BinnedValues.cc src/BinnedValuesJSONParser.cc)
target_compile_features(BinnedValues PRIVATE cxx_std_17)
target_include_directories(BinnedValues PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(BinnedValues ROOT::Hist Boost::boost)
install(TARGETS BinnedValues LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_library(BambooLumiMask SHARED src/LumiMask.cc)
target_compile_features(BambooLumiMask PRIVATE cxx_std_17)
target_include_directories(BambooLumiMask PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(BambooLumiMask Boost::boost)
install(TARGETS BambooLumiMask LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/jetclasses_filenames.txt FILENAMES_CMSJET_INC REGEX ".*/interface/.*")
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/jetclasses_filenames.txt FILENAMES_CMSJET_SRC REGEX ".*/src/.*")
string(REGEX REPLACE "[^;]+/interface/([^;]+)" "${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include/\\1" INC_CMSJET "${FILENAMES_CMSJET_INC}")
string(REGEX REPLACE "[^;]+/src/([^;]+)" "${CMAKE_CURRENT_BINARY_DIR}/CMSJet/src/\\1" SRC_CMSJET "${FILENAMES_CMSJET_SRC}")
add_custom_target(CMSJetSources COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/getjetclasses.sh COMMENT "Copying jet classes from CMSSW" BYPRODUCTS ${INC_CMSJET} ${SRC_CMSJET})
install(FILES ${INC_CMSJET} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(JMEObjects SHARED ${SRC_CMSJET} src/JMESystematicsCalculator.cc)
add_dependencies(JMEObjects CMSJetSources)
target_compile_features(JMEObjects PRIVATE cxx_std_17)
target_include_directories(JMEObjects PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/CMSJet/src> $<INSTALL_INTERFACE:include>)
target_link_libraries(JMEObjects ROOT::Hist ROOT::ROOTVecOps Boost::boost)
install(TARGETS JMEObjects LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_library(RoccoR SHARED src/RoccoR.cc src/RochesterCorrectionCalculator.cc)
target_compile_features(RoccoR PRIVATE cxx_std_17)
target_include_directories(RoccoR PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(RoccoR ROOT::ROOTVecOps Boost::boost)
install(TARGETS RoccoR LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

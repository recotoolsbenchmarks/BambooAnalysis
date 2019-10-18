cmake_minimum_required( VERSION 2.8 )

project(BambooExt LANGUAGES CXX)

include(GNUInstallDirs)
set(CMAKE_INSTALL_LIBDIR "lib")

if(EXISTS "$ENV{CMAKE_PREFIX_PATH}/include/boost")
    message(STATUS "Will use $ENV{CMAKE_PREFIX_PATH} as base path for boost")
    set(BOOST_ROOT $ENV{CMAKE_PREFIX_PATH})
    set(Boost_NO_BOOST_CMAKE ON)
endif()

find_package(ROOT REQUIRED COMPONENTS Hist)
find_package(Boost REQUIRED)

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

set(SRC_CMSJET
  CMSJet/src/JetCorrectionUncertainty.cc
  CMSJet/src/JetCorrectorParameters.cc
  CMSJet/src/JetCorrectorParametersHelper.cc
  CMSJet/src/JetResolution.cc
  CMSJet/src/JetResolutionObject.cc
  CMSJet/src/SimpleJetCorrectionUncertainty.cc
  CMSJet/src/SimpleJetCorrector.cc
  CMSJet/src/FactorizedJetCorrector.cc
  CMSJet/src/FactorizedJetCorrectorCalculator.cc
  CMSJet/src/formulaEvaluatorBase.cc
  CMSJet/src/formulaConstantEvaluator.cc
  CMSJet/src/formulaParameterEvaluator.cc
  CMSJet/src/formulaVariableEvaluator.cc
  CMSJet/src/FormulaEvaluator.cc
  )
set(CMSJET_HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include/FactorizedJetCorrectorCalculator.h
  ${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include/FactorizedJetCorrector.h
  ${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include/FormulaEvaluator.h
  ${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include/JetCorrectionUncertainty.h
  ${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include/JetCorrectorParameters.h
  ${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include/JetCorrectorParametersHelper.h
  ${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include/JetResolution.h
  ${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include/JetResolutionObject.h
  ${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include/SimpleJetCorrectionUncertainty.h
  ${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include/SimpleJetCorrector.h
  ${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include/Utilities.h
  )
add_custom_target(CMSJetSources COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/getjetclasses.sh COMMENT "Copying jet classes from CMSSW" BYPRODUCTS ${SRC_CMSJET} ${CMSJET_HEADERS})
install(FILES ${CMSJET_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(JMEObjects SHARED ${SRC_CMSJET} src/JMESystematicsCalculator.cc ${CMSJET_HEADERS})
add_dependencies(JMEObjects CMSJetSources)
target_compile_features(JMEObjects PRIVATE cxx_std_17)
target_include_directories(JMEObjects PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(JMEObjects ROOT::Hist Boost::boost)
install(TARGETS JMEObjects LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_library(RoccoR SHARED src/RoccoR.cc src/RochesterCorrectionCalculator.cc)
target_compile_features(RoccoR PRIVATE cxx_std_17)
target_include_directories(RoccoR PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(RoccoR Boost::boost)
install(TARGETS RoccoR LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

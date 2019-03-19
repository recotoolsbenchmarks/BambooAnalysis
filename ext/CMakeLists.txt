cmake_minimum_required( VERSION 2.8 )

project(BambooExt LANGUAGES CXX)

include(GNUInstallDirs)
set(CMAKE_INSTALL_LIBDIR "lib")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

find_package( ROOT REQUIRED )
find_package( Boost REQUIRED ) # property_tree is header-only

message(STATUS "Found Boost with include dir ${Boost_INCLUDE_DIR} (change by passing e.g. '-DBoost_NO_BOOST_CMAKE=ON -DBOOST_ROOT=$(scram tool tag boost BOOST_BASE)')")

include_directories(${ROOT_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/interface)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(BinnedValues SHARED src/BinnedValues.cc src/BinnedValuesJSONParser.cc)
target_compile_features(BinnedValues PRIVATE cxx_std_17)
target_include_directories(BinnedValues PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(BinnedValues ${Boost_LIBRARIES} ${ROOT_LIBRARIES})
install(TARGETS BinnedValues
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

add_library(BambooLumiMask SHARED src/LumiMask.cc)
target_compile_features(BambooLumiMask PRIVATE cxx_std_17)
target_include_directories(BambooLumiMask PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(BambooLumiMask ${Boost_LIBRARIES})
install(TARGETS BambooLumiMask
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

set(SRC_CMSJET
  src/JetCorrectionUncertainty.cc
  src/JetCorrectorParameters.cc
  src/JetCorrectorParametersHelper.cc
  src/JetResolution.cc
  src/JetResolutionObject.cc
  src/SimpleJetCorrectionUncertainty.cc
  src/SimpleJetCorrector.cc
  src/FactorizedJetCorrector.cc
  src/FactorizedJetCorrectorCalculator.cc
  src/formulaEvaluatorBase.cc
  src/formulaConstantEvaluator.cc
  src/formulaParameterEvaluator.cc
  src/formulaVariableEvaluator.cc
  src/FormulaEvaluator.cc
  )
add_library(JMEObjects SHARED ${SRC_CMSJET} src/JMESystematicsCalculator.cc)
target_compile_features(JMEObjects PRIVATE cxx_std_17)
target_include_directories(JMEObjects PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
target_link_libraries(JMEObjects ${ROOT_LIBRARIES})
install(TARGETS JMEObjects
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

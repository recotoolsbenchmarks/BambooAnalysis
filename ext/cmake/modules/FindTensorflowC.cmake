# Try to find the Tensorflow C bindings
# This will define
#  TENSORFLOWC_FOUND
#  TENSORFLOWC_INCLUDE_DIRS
#  TENSORFLOWC_LIBRARIES

find_path(TENSORFLOWC_INCLUDE_DIR tensorflow/c/c_api.h)
find_library(TENSORFLOWC_LIBRARY tensorflow)
find_library(TENSORFLOWC_FRAMEWORK_LIBRARY tensorflow_framework)
set(TENSORFLOWC_LIBRARIES ${TENSORFLOWC_LIBRARY} ${TENSORFLOWC_FRAMEWORK_LIBRARY})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set TENSORFLOWC_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(TensorflowC DEFAULT_MSG TENSORFLOWC_LIBRARIES TENSORFLOWC_INCLUDE_DIR)
mark_as_advanced(TENSORFLOWC_INCLUDE_DIR TENSORFLOWC_LIBRARY )

add_library(TensorflowC SHARED IMPORTED)
set_property(TARGET TensorflowC PROPERTY IMPORTED_LOCATION "${TENSORFLOWC_LIBRARY}")
set_property(TARGET TensorflowC PROPERTY INTERFACE_LINK_LIBRARIES "${TENSORFLOWC_FRAMEWORK_LIBRARY}")
set_property(TARGET TensorflowC PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${TENSORFLOWC_INCLUDE_DIR}")

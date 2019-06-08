INCLUDE(test_macros)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/fem/src)

CONFIGURE_FILE(case2layers_baseline.sif case2layers_baseline.sif COPYONLY)

file(COPY ELMERSOLVER_STARTINFO square.grd DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

ADD_ELMER_TEST(ExtrudedMeshLayers)
ADD_ELMER_LABEL(ExtrudedMeshLayers quick)

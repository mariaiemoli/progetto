SET(SRC_DIR ${PACS_SOURCE_DIR}/src)
SET(HEADERS_DIR ${PACS_SOURCE_DIR}/include)

INCLUDE_DIRECTORIES( ${PACS_SOURCE_DIR}/include)
LINK_DIRECTORIES(${PACS_SOURCE_DIR}/src)


EXECUTE_PROCESS( COMMAND ln -s ${CMAKE_CURRENT_SOURCE_DIR}/Bifurcation ${CMAKE_CURRENT_BINARY_DIR}/Bifurcation )



ADD_EXECUTABLE( Test2 Test2.cc ${SRC_DIR}/BC.cc ${SRC_DIR}/BCHandler.cc ${SRC_DIR}/DarcyFractured.cc ${SRC_DIR}/Exporter.cc ${SRC_DIR}/FractureData.cc ${SRC_DIR}/FractureHandler.cc ${SRC_DIR}/FractureIntersect.cc ${SRC_DIR}/FracturesSet.cc ${SRC_DIR}/IntersectData.cc ${SRC_DIR}/LevelSetData.cc ${SRC_DIR}/LevelSetHandler.cc ${SRC_DIR}/MatrixBifurcationHandler.cc ${SRC_DIR}/MediumData.cc ${SRC_DIR}/MeshHandler.cc ${SRC_DIR}/Parser.cc ${SRC_DIR}/PointData.cc ${SRC_DIR}/StringUtility.cc ${SRC_DIR}/TriangleData.cc ${SRC_DIR}/QuadrilaterData.cc ${SRC_DIR}/Geometry.cc ${SRC_DIR}/UsefulFunctions.cc ${SRC_DIR}/XFEMOperators.cc  ${HEADERS_DIR}/BC.h ${HEADERS_DIR}/BCHandler.h ${HEADERS_DIR}/DarcyFractured.h ${HEADERS_DIR}/Exporter.h ${HEADERS_DIR}/FractureData.h ${HEADERS_DIR}/FractureHandler.h ${HEADERS_DIR}/FractureIntersect.h ${HEADERS_DIR}/FracturesSet.h ${HEADERS_DIR}/IntersectData.h ${HEADERS_DIR}/LevelSetData.h ${HEADERS_DIR}/LevelSetHandler.h ${HEADERS_DIR}/MatrixBifurcationHandler.h ${HEADERS_DIR}/MediumData.h ${HEADERS_DIR}/MeshHandler.h ${HEADERS_DIR}/Parser.h ${HEADERS_DIR}/PointData.h ${HEADERS_DIR}/QuadrilaterData.h ${HEADERS_DIR}/StringUtility.h ${HEADERS_DIR}/TriangleData.h ${HEADERS_DIR}/Geometry.h ${HEADERS_DIR}/UsefulFunctions.h ${HEADERS_DIR}/XFEMOperators.h ${HEADERS_DIR}/Core.h ${HEADERS_DIR}/GetPot ${HEADERS_DIR}/ParserDefinitions.h ${HEADERS_DIR}/ParserSpiritGrammar.h )

TARGET_LINK_LIBRARIES(Test2 ${GETFEM_LIBRARIES} ${QHULL_LIBRARIES} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} )


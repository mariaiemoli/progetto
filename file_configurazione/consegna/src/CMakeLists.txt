INCLUDE_DIRECTORIES( ${PACS_SOURCE_DIR}/include)
SET(HEADERS_DIR ${PACS_SOURCE_DIR}/include)

SET(TARGET_SRC BC.cc BCHandler.cc DarcyFractured.cc Exporter.cc FractureData.cc FractureHandler.cc FractureIntersect.cc FracturesSet.cc IntersectData.cc LevelSetData.cc LevelSetHandler.cc MatrixBifurcationHandler.cc MediumData.cc MeshHandler.cc Parser.cc PointData.cc QuadrilaterData.cc StringUtility.cc TriangleData.cc Geometry.cc UsefulFunctions.cc XFEMOperators.cc  )

SET( TARGET_HANDLER  ${HEADERS_DIR}/BC.h ${HEADERS_DIR}/BCHandler.h ${HEADERS_DIR}/DarcyFractured.h ${HEADERS_DIR}/Exporter.h ${HEADERS_DIR}/FractureData.h ${HEADERS_DIR}/FractureHandler.h ${HEADERS_DIR}/FractureIntersect.h ${HEADERS_DIR}/FracturesSet.h ${HEADERS_DIR}/IntersectData.h ${HEADERS_DIR}/LevelSetData.h ${HEADERS_DIR}/LevelSetHandler.h ${HEADERS_DIR}/MatrixBifurcationHandler.h ${HEADERS_DIR}/MediumData.h ${HEADERS_DIR}/MeshHandler.h ${HEADERS_DIR}/Parser.h ${HEADERS_DIR}/PointData.h ${HEADERS_DIR}/QuadrilaterData.h ${HEADERS_DIR}/StringUtility.h ${HEADERS_DIR}/TriangleData.h ${HEADERS_DIR}/Geometry.h ${HEADERS_DIR}/UsefulFunctions.h ${HEADERS_DIR}/XFEMOperators.h ${HEADERS_DIR}/Core.h ${HEADERS_DIR}/GetPot ${HEADERS_DIR}/ParserDefinitions.h ${HEADERS_DIR}/ParserSpiritGrammar.h)




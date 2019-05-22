IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

FIND_PATH(MKNIX_LINK_DIR HINTS ${PROJECT_SOURCE_DIR}/../MkniX-build/src  NAMES libmknix.so)
FIND_PATH(MKNIX_SOURCE_DIR HINTS ${PROJECT_SOURCE_DIR}/../MkniX)

include_directories (
${MkniX-LamellaRT_SOURCE_DIR}/src-coarse 
${MKNIX_SOURCE_DIR}/src/simulation/
${MKNIX_SOURCE_DIR}/src/core/
${MKNIX_SOURCE_DIR}/src/
${MKNIX_SOURCE_DIR}/LMX-1.0/src/
)

add_executable (alicia ALICIA.cpp)

link_directories ( ${PROJECT_BINARY_DIR})

target_link_libraries ( alicia debug
 ${MKNIX_LINK_DIR}/../../MkniX-debug/src/libmknix.so 
)
target_link_libraries ( alicia optimized
 ${MKNIX_LINK_DIR}/libmknix.so
)

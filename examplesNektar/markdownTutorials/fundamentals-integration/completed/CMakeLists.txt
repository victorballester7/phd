cmake_minimum_required(VERSION 2.8)
PROJECT(Nektar-exercises)

FIND_PACKAGE(Nektar++ PATHS ../../../builds)

INCLUDE_DIRECTORIES(${NEKTAR++_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${NEKTAR++_TP_INCLUDE_DIRS})
LINK_DIRECTORIES(${NEKTAR++_LIBRARY_DIRS})
LINK_DIRECTORIES(${NEKTAR++_TP_LIBRARY_DIRS})
ADD_DEFINITIONS(${NEKTAR++_DEFINITIONS})

SET(LinkLibraries StdRegions)

SET(NekExerciseDir ${PROJECT_SOURCE_DIR})

SET(StdIntegration1DSource  ${NekExerciseDir}/StdIntegration1D.cpp)
ADD_EXECUTABLE(StdIntegration1D ${StdIntegration1DSource})
TARGET_LINK_LIBRARIES(StdIntegration1D ${NEKTAR++_LIBRARIES})

SET(StdIntegration2DSource  ${NekExerciseDir}/StdIntegration2D.cpp)
ADD_EXECUTABLE(StdIntegration2D ${StdIntegration2DSource})
TARGET_LINK_LIBRARIES(StdIntegration2D ${NEKTAR++_LIBRARIES})

SET(LocIntegration2DSource  ${NekExerciseDir}/LocIntegration2D.cpp)
ADD_EXECUTABLE(LocIntegration2D ${LocIntegration2DSource})
TARGET_LINK_LIBRARIES(LocIntegration2D ${NEKTAR++_LIBRARIES})


set(crackVisLT_SOURCES
    vtkCellCenter3D.cpp
    vtkCellCenter3D.h
    CombinationsUtils.cpp
    CellSearch.cpp
    Timeris.cpp
    CombinationsUtils.h
    CellSearch.h
    Timeris.h
    vtkVoronoi3D.cpp
    vtkVoronoi3D.h
    VoroCPPManager.cpp
    VoroCPPManager.h
    vtkThresholdPolyData.cpp
    vtkThresholdPolyData.h

)
add_library(algorithmlibs SHARED ${crackVisLT_SOURCES})




ADD_PARAVIEW_PLUGIN(CrackAlgorithms "1.0"
SERVER_MANAGER_XML CrackAlgorithms.xml
SERVER_MANAGER_SOURCES CrackAlgorithms.cxx )
target_link_libraries(CrackAlgorithms LINK_PRIVATE vtksys  -lstdc++fs algorithmlibs )

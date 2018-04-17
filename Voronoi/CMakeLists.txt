ADD_PARAVIEW_PLUGIN(Voronoi "1.0"
SERVER_MANAGER_XML Voronoi.xml
SERVER_MANAGER_SOURCES Voronoi.cxx)
target_link_libraries(Voronoi LINK_PRIVATE vtksys  -lstdc++fs )

ADD_PARAVIEW_PLUGIN(Porocity "1.0"
SERVER_MANAGER_XML Porocity.xml
SERVER_MANAGER_SOURCES Porocity.cxx)
target_link_libraries(Porocity LINK_PRIVATE vtksys  -lstdc++fs )

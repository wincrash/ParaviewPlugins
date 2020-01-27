ADD_PARAVIEW_PLUGIN(MicroStructure "1.0"
SERVER_MANAGER_XML MicroStructure.xml
SERVER_MANAGER_SOURCES MicroStructure.cxx)
target_link_libraries(MicroStructure LINK_PRIVATE vtksys  -lstdc++fs )

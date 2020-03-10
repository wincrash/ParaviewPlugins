
ADD_PARAVIEW_PLUGIN(HDF5ReaderConvex "1.0"
SERVER_MANAGER_XML HDF5ReaderConvex_Server.xml
SERVER_MANAGER_SOURCES HDF5ReaderConvex.cxx)
target_link_libraries(HDF5ReaderConvex LINK_PRIVATE vtksys  -lstdc++fs )


ADD_PARAVIEW_PLUGIN(HDF5ReaderMultiBlock "1.0"
SERVER_MANAGER_XML HDF5ReaderMultiBlock_Server.xml
SERVER_MANAGER_SOURCES HDF5ReaderMultiBlock.cxx)
target_link_libraries(HDF5ReaderMultiBlock LINK_PRIVATE vtksys  -lstdc++fs )

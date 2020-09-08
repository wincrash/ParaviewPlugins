
ADD_PARAVIEW_PLUGIN(McdemReader "1.0"
SERVER_MANAGER_XML McdemReader_Server.xml
SERVER_MANAGER_SOURCES McdemReader.cxx)
target_link_libraries(McdemReader LINK_PRIVATE vtksys  -lstdc++fs ${HDF5_LIBRARIES} -lhdf5_cpp  )

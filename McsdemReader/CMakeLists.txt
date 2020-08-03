
ADD_PARAVIEW_PLUGIN(McsdemReader "1.0"
SERVER_MANAGER_XML McsdemReader_Server.xml
SERVER_MANAGER_SOURCES McsdemReader.cxx)
target_link_libraries(McsdemReader LINK_PRIVATE vtksys  -lstdc++fs ${HDF5_LIBRARIES} -lhdf5_cpp  )

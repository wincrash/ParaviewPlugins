
ADD_PARAVIEW_PLUGIN(HDF5Reader "1.0"
SERVER_MANAGER_XML HDF5Reader_Server.xml
SERVER_MANAGER_SOURCES HDF5Reader.cxx)
target_link_libraries(HDF5Reader LINK_PRIVATE vtksys  -lstdc++fs )

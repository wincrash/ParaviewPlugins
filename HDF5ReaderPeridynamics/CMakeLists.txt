
ADD_PARAVIEW_PLUGIN(HDF5ReaderPeridynamics "1.0"
SERVER_MANAGER_XML HDF5ReaderPeridynamics_Server.xml
SERVER_MANAGER_SOURCES HDF5ReaderPeridynamics.cxx)
target_link_libraries(HDF5ReaderPeridynamics LINK_PRIVATE vtksys  -lstdc++fs )

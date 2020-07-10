
ADD_PARAVIEW_PLUGIN(HDF5ReaderBonds "1.0"
SERVER_MANAGER_XML HDF5ReaderBonds_Server.xml
SERVER_MANAGER_SOURCES HDF5ReaderBonds.cxx)
target_link_libraries(HDF5ReaderBonds LINK_PRIVATE vtksys  -lstdc++fs ${HDF5_LIBRARIES} -lhdf5_cpp  )

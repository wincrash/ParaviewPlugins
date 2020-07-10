
ADD_PARAVIEW_PLUGIN(HDF5ReaderTetraBonds "1.0"
SERVER_MANAGER_XML HDF5ReaderTetraBonds_Server.xml
SERVER_MANAGER_SOURCES HDF5ReaderTetraBonds.cxx)
target_link_libraries(HDF5ReaderTetraBonds LINK_PRIVATE vtksys  -lstdc++fs ${HDF5_LIBRARIES} -lhdf5_cpp  )

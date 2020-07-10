
ADD_PARAVIEW_PLUGIN(HDF5ReaderTetraConvex "1.0"
SERVER_MANAGER_XML HDF5ReaderTetraConvex_Server.xml
SERVER_MANAGER_SOURCES HDF5ReaderTetraConvex.cxx)
target_link_libraries(HDF5ReaderTetraConvex LINK_PRIVATE vtksys  -lstdc++fs ${HDF5_LIBRARIES} -lhdf5_cpp  )

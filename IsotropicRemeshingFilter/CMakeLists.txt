find_package(CGAL)
if(CGAL_FOUND)
  include( ${CGAL_USE_FILE} )
endif(CGAL_FOUND)

ADD_PARAVIEW_PLUGIN(IsotropicRemeshingFilter "1.0"
    SERVER_MANAGER_XML IsotropicRemeshingFilter.xml
    SERVER_MANAGER_SOURCES vtkIsotropicRemeshingFilter.cxx)
    target_link_libraries(IsotropicRemeshingFilter LINK_PRIVATE
      CGAL::CGAL 
      ${Boost_LIBRARIES})
    

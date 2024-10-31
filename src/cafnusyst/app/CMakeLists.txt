LIST(APPEND TARGETS_TO_BUILD 
UpdateReweight
)

foreach(targ ${TARGETS_TO_BUILD})
  add_executable(${targ} ${targ}.cxx)
  target_link_libraries(${targ} 
    nusyst::systproviders
    ${ROOT_LIBRARIES}
    ROOT::Geom
    ROOT::MathMore
    ${ANAOBJ_LIBS}
    cafnusyst::commondeps
    cafnusyst::interface
    cafnusyst::utility
  )
  target_include_directories(${targ}
    PRIVATE
    $ENV{SRPROXY_INC}
  )

endforeach()

install(TARGETS ${TARGETS_TO_BUILD} DESTINATION bin/)

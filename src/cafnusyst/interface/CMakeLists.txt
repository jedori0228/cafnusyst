SET(IFCE_IMPLFILES
  WeightUpdater.cxx
)

SET(IFCE_HDRFILES
  WeightUpdater.h
)


add_library(cafnusyst_interface SHARED ${IFCE_IMPLFILES})
add_library(cafnusyst::interface ALIAS cafnusyst_interface)

set_target_properties(cafnusyst_interface PROPERTIES 
  PUBLIC_HEADER "${IFCE_HDRFILES}"
  EXPORT_NAME interface )

target_link_libraries(cafnusyst_interface
PUBLIC
  nusyst::systproviders
  ${ROOT_LIBRARIES}
  ROOT::Geom
  ROOT::MathMore
  ${ANAOBJ_LIBS}
  cafnusyst::commondeps
  cafnusyst::utility
)
target_include_directories(cafnusyst_interface
  PRIVATE
  $ENV{SRPROXY_INC}
)


#install(TARGETS ${TARGETS_TO_BUILD} DESTINATION interface/)

install(TARGETS cafnusyst_interface
    EXPORT cafnusyst-targets
    PUBLIC_HEADER DESTINATION include/cafnusyst/interface COMPONENT Development)

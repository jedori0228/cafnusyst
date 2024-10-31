SET(IFCE_IMPLFILES
  Utilities.cxx
)

SET(IFCE_HDRFILES
  Utilities.h
)


add_library(cafnusyst_utility SHARED ${IFCE_IMPLFILES})
add_library(cafnusyst::utility ALIAS cafnusyst_utility)

set_target_properties(cafnusyst_utility PROPERTIES 
  PUBLIC_HEADER "${IFCE_HDRFILES}"
  EXPORT_NAME utility )

target_link_libraries(cafnusyst_utility
PUBLIC
  nusyst::systproviders
  ${ROOT_LIBRARIES}
  ROOT::Geom
  ROOT::MathMore
  ${ANAOBJ_LIBS}
  cafnusyst::commondeps
)
target_include_directories(cafnusyst_utility
  PRIVATE
  $ENV{SRPROXY_INC}
)

install(TARGETS cafnusyst_utility
    EXPORT cafnusyst-targets
    PUBLIC_HEADER DESTINATION include/cafnusyst/utility COMPONENT Development)

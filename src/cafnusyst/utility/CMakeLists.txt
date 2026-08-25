SET(IFCE_IMPLFILES
  Utilities.cxx
  NuTreeHelper.cxx
)

SET(IFCE_HDRFILES
  Utilities.h
  NuTreeHelper.h
  ResourceMonitor.h
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
PRIVATE
  # See the ANAOBJ_PROXY_LIB/ANAOBJ_FLAT_LIB comment in the top-level
  # CMakeLists.txt: PRIVATE so this does NOT propagate to cafnusyst_utility's
  # own consumers (cafnusyst_interface, the app executables) -- they already
  # depend on cafnusyst::utility and pick up these symbols dynamically from
  # it. Linking this static archive into more than one shared
  # library/executable in the same process causes duplicate ROOT dictionary
  # registration and a heap-corruption crash at exit.
  ${ANAOBJ_PROXY_LIB}
)
target_include_directories(cafnusyst_utility
  PRIVATE
  $ENV{SRPROXY_INC}
)

install(TARGETS cafnusyst_utility
    EXPORT cafnusyst-targets
    PUBLIC_HEADER DESTINATION include/cafnusyst/utility COMPONENT Development)

SET(IFCE_IMPLFILES
  WeightUpdater.cxx
  LoopRecord.cxx
)

SET(IFCE_HDRFILES
  WeightUpdater.h
  LoopRecord.h
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
PRIVATE
  # See the ANAOBJ_PROXY_LIB/ANAOBJ_FLAT_LIB comment in the top-level
  # CMakeLists.txt: WeightUpdater.cxx is the only place that references
  # caf::FlatStandardRecord (flat::Flat<caf::StandardRecord>) directly, so
  # this static archive belongs here -- PRIVATE so it isn't also embedded
  # into whatever links cafnusyst::interface (the app executables), which
  # would cause duplicate ROOT dictionary registration and a heap-corruption
  # crash at exit.
  ${ANAOBJ_FLAT_LIB}
)
target_include_directories(cafnusyst_interface
  PRIVATE
  $ENV{SRPROXY_INC}
)


#install(TARGETS ${TARGETS_TO_BUILD} DESTINATION interface/)

install(TARGETS cafnusyst_interface
    EXPORT cafnusyst-targets
    PUBLIC_HEADER DESTINATION include/cafnusyst/interface COMPONENT Development)

cmake_minimum_required (VERSION 3.8)
add_executable(RadialDistributionFunction main.cpp)
if(NOT DONT_USE_CUDA)
  SET_SOURCE_FILES_PROPERTIES(main.cpp PROPERTIES LANGUAGE CUDA)
endif()

set_target_properties(RadialDistributionFunction PROPERTIES OUTPUT_NAME "rdf")

install(TARGETS RadialDistributionFunction  RUNTIME DESTINATION bin)

configure_file ( "defines.h.in" ${CMAKE_CURRENT_SOURCE_DIR}/defines.h  )


find_package(Git)
if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
  execute_process(COMMAND bash  ${PROJECT_SOURCE_DIR}/cmake/gitversion.sh OUTPUT_VARIABLE GITCOMMIT)
else()
  SET(GITCOMMIT "unknown")
endif()

configure_file("gitversion.h.in" "${CMAKE_CURRENT_SOURCE_DIR}/gitversion.h")

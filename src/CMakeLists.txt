# compile everything except the main source file into a statis lib to be linked
# to the main executable
add_library(${LIBRARY_NAME} STATIC
vac_splitting.cpp
qcd_const.cpp
convolution.cpp
	)

set_target_properties(${LIBRARY_NAME} PROPERTIES PREFIX "")

install(TARGETS ${LIBRARY_NAME} DESTINATION lib/${LIBRARY_NAME})
# install the library headers
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/${LIBRARY_NAME})

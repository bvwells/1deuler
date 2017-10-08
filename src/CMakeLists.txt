########################################
# Set up how to compile the source files
########################################

# Define the executable name
SET(1DEULEREXE 1deuler.exe)

# Add the source files
SET(1DEULER_Sources ${SRC}/main.f90
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${1DEULEREXE} ${1DEULER_Sources})

# This links foo to the bar library
#TARGET_LINK_LIBRARIES(${1DEULEREXE} ${SOME_LIBRARY})

#####################################
# Tell how to install this executable
#####################################

set_target_properties( ${1DEULEREXE}
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)

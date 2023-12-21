cmake_minimum_required(VERSION 3.1)
project(example)

# Find Geant4 package
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
else()
    find_package(Geant4 REQUIRED)
endif()

# Find ROOT package
find_package(ROOT REQUIRED COMPONENTS RIO Tree)

# Setup Geant4 include directories and compile definitions
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})

# Locate sources and headers for this project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Add the executable, and link it to the Geant4 and ROOT libraries
add_executable(example example.cc ${sources} ${headers})
target_link_libraries(example ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

# Copy scripts to the build directory
set(EXAMPLE_SCRIPTS
  init_vis.mac
  shoot.mac
  vis.mac
)

foreach(_script ${EXAMPLE_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
install(TARGETS example DESTINATION bin)

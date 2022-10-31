# Motocyte configuration

set(CMAKE_BUILD_TYPE Release) # enable 10x faster vtk visualizer
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(motocyte_sources
		VolumeSpecies.cpp
    Walker.cpp
		HDF5File.cpp
		HDF5Interface.cpp
		MacroSpecies.cpp
		MesoSpecies.cpp
		MicroSpace.cpp
		MesoSpace.cpp
		MicroSpecies.cpp
    Model.cpp
    MotileSpecies.cpp
    Motocyte.cpp
    PointSpecies.cpp
		Populator.cpp
    Process.cpp
    pugixml.cpp
    Reaction.cpp
    Space.cpp
    SpaceCompartment.cpp
    Species.cpp
    Stepper.cpp
    XMLModelReader.cpp
    XMLModelWriter.cpp
    VisualWriter.cpp
    NumberWriter.cpp
    TrackWriter.cpp
    TrackReader.cpp
    TrackSpecies.cpp
    TrackIDSpecies.cpp
    XDMFFile.cpp)

add_executable(motocyte ${motocyte_sources})

target_compile_options(motocyte PUBLIC -std=c++17 -O3 -g -Wall)

target_include_directories(motocyte PUBLIC .)

find_package(HDF5 REQUIRED COMPONENTS CXX)

find_package(Boost COMPONENTS filesystem system REQUIRED)

find_package(GSL REQUIRED) 

target_include_directories(motocyte PUBLIC
  ${Boost_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIR})

target_link_libraries(motocyte PUBLIC
  ${Boost_LIBRARIES}
  ${HDF5_LIBRARIES}
  ${GSL_LIBRARIES})


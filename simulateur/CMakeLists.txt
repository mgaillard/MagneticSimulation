add_executable(MagneticSimulation)

message(STATUS "Creating target 'MagneticSimulation'")

# Files for the project

set(HEADER_FILES
    constants.h
    consoleapplication.h
    export_vtk.h
    interfacexml.h
    renduimage.h
    scene.h
    shape.h
    shapecircle.h
    shapepolygon.h
    simulateur.h
    solver.h
    solvercolpivhouseholderqr.h
    solverdense.h
    solversparse.h
    solversparselu.h
)

set(SRC_FILES
    consoleapplication.cpp
    export_vtk.cpp
    interfacexml.cpp
    main.cpp
    renduimage.cpp
    scene.cpp
    shape.cpp
    shapecircle.cpp
    shapepolygon.cpp
    simulateur.cpp
    simulateur.cpp
    solver.cpp
    solvercolpivhouseholderqr.cpp
    solverdense.cpp
    solversparse.cpp
    solversparselu.cpp
)

# Setup filters in Visual Studio
source_group("Header Files" FILES ${HEADER_FILES})
source_group("Source Files" FILES ${SRC_FILES})

# Copy the ReSharper config file
configure_file(
    "${CMAKE_SOURCE_DIR}/MagneticSimulation.sln.DotSettings"
    ${CMAKE_BINARY_DIR} COPYONLY
)

target_sources(MagneticSimulation
    PUBLIC
    ${HEADER_FILES}
    PRIVATE
    ${SRC_FILES}
)

target_link_libraries(MagneticSimulation 
    PRIVATE
    OpenMP::OpenMP_CXX
    Qt5::Core
    Qt5::Gui
    Qt5::Xml
    Eigen3::Eigen
)

if(TARGET Qt5::windeployqt)
    add_custom_command(TARGET MagneticSimulation
        POST_BUILD
        COMMAND set PATH=%PATH%$<SEMICOLON>${qt5_install_prefix}/bin
        COMMAND Qt5::windeployqt "$<TARGET_FILE_DIR:MagneticSimulation>/$<TARGET_FILE_NAME:MagneticSimulation>"
    )
endif()

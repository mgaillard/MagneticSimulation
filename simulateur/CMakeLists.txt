add_executable(MagneticSimulation)

message(STATUS "Creating target 'MagneticSimulation'")

# Files for the project

set(HEADER_FILES
    constants.h
    consoleapplication.h
    contour.h
    export_image.h
    export_vtk.h
    interfacexml.h
    magneticsimulation.h
    scene.h
    shape.h
    shapecircle.h
    shapepolygon.h
    solver.h
    solvercolpivhouseholderqr.h
    solverdense.h
    solversparse.h
    solversparselu.h
    streamlines.h
)

set(SRC_FILES
    consoleapplication.cpp
    contour.cpp
    export_image.cpp
    export_vtk.cpp
    interfacexml.cpp
    main.cpp
    magneticsimulation.cpp
    scene.cpp
    shape.cpp
    shapecircle.cpp
    shapepolygon.cpp
    solver.cpp
    solvercolpivhouseholderqr.cpp
    solverdense.cpp
    solversparse.cpp
    solversparselu.cpp
    streamlines.cpp
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
    spdlog::spdlog
    OpenMP::OpenMP_CXX
    Eigen3::Eigen
    Qt5::Core
    Qt5::Gui
    Qt5::Xml
)

if(TARGET spdlog::spdlog)
    add_custom_command(TARGET MagneticSimulation
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:spdlog>" "$<TARGET_FILE_DIR:MagneticSimulation>"
    )
endif()

if(TARGET Qt5::windeployqt)
    add_custom_command(TARGET MagneticSimulation
        POST_BUILD
        COMMAND set PATH=%PATH%$<SEMICOLON>${qt5_install_prefix}/bin
        COMMAND Qt5::windeployqt "$<TARGET_FILE_DIR:MagneticSimulation>/$<TARGET_FILE_NAME:MagneticSimulation>"
    )
endif()

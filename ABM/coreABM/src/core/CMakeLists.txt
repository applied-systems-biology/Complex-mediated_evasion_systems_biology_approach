add_subdirectory(analyser)
add_subdirectory(basic)
add_subdirectory(simulation)
add_subdirectory(utils)
add_subdirectory(visualisation)

add_library(core SHARED
        $<TARGET_OBJECTS:analyser>
        $<TARGET_OBJECTS:basic>
        $<TARGET_OBJECTS:simulation>
        $<TARGET_OBJECTS:utils>
        $<TARGET_OBJECTS:visualisation>
        )

add_library(abm::core ALIAS core)
target_include_directories(core PRIVATE ${PROJECT_SOURCE_DIR}/src)

target_link_libraries(analyser PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
target_link_libraries(basic PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
target_link_libraries(simulation PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
target_link_libraries(utils PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
target_link_libraries(visualisation PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)


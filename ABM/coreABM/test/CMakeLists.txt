add_executable(test_configurations  src/testConfigurations.cpp)
target_include_directories(test_configurations PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(test_configurations PRIVATE
        project_options
        simulation
        analyser
        basic
        utils
        visualisation
        simulatorExample
        simulatorCME
        Boost::filesystem)

add_executable(test_units  src/testUnits.cpp)
target_include_directories(test_units PRIVATE ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(test_units PRIVATE
        project_options
        analyser
        simulation
        basic
        utils
        visualisation
        simulatorExample
        simulatorCME
        Boost::filesystem)

add_dependencies(test_units simulation)
add_dependencies(test_units analyser)
add_dependencies(test_units basic)
add_dependencies(test_units utils)
add_dependencies(test_units visualisation)

add_test(NAME configurations_functions_tests COMMAND test_configurations)
add_test(NAME analyser_functions_tests COMMAND test_units)
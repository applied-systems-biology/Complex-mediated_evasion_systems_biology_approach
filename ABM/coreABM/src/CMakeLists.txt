add_subdirectory(core)
add_subdirectory(apps/example)
add_subdirectory(apps/CME)

target_link_libraries(core PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
set_target_properties(core PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(simulatorExample PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)
target_link_libraries(simulatorCME PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)


add_executable(coreABM main.cpp)
target_link_libraries(coreABM PUBLIC
project_options
project_warnings
abm::core
simulatorCME
simulatorExample
Boost::filesystem
OpenMP::OpenMP_CXX)

target_include_directories(coreABM PRIVATE ${PROJECT_SOURCE_DIR}/src)

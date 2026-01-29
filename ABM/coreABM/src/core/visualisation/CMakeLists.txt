add_library(visualisation SHARED
        PovFile.cpp
        PovRayObject.cpp
        Visualizer.cpp)
target_include_directories(visualisation PRIVATE  ${PROJECT_SOURCE_DIR}/src)
add_library(analyser SHARED
        histogram_measurement.cpp
        pair_measurment.cpp
        Analyser.cpp
        InSituMeasurements.cpp)
target_include_directories(analyser PRIVATE ${PROJECT_SOURCE_DIR}/src)
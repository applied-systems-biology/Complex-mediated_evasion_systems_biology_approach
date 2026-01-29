add_library(simulatorExample SHARED
        SimulatorExample.cpp
        CuboidSiteExample.cpp
        io_utils_example.cpp
        FungalCellExample.cpp
        CellFactoryExample.cpp
        AnalyserExample.cpp
        InSituMeasurementsExample.cpp
        )

target_include_directories(simulatorExample PRIVATE ../..)
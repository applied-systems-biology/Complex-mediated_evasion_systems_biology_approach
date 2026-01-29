add_library(simulatorCME SHARED
        SimulatorCME.cpp
        CuboidSiteCME.cpp
        io_utilsCME.cpp
        factories/CellFactoryCME.cpp
        factories/InteractionStateFactoryCME.cpp
        factories/InteractionFactoryCME.cpp
        interaction/ComplexFormationInteraction.cpp
        interaction/ComplexExtensionInteraction.cpp
        interaction/UptakeInteraction.cpp
        interactiontypes/ComplexFormation.cpp
        interactiontypes/ComplexExtension.cpp
        interactiontypes/Uptake.cpp
        cells/Amp.cpp
        cells/Complex.cpp
        cells/Defensive.cpp
        cells/Pathogenic_cell.cpp
        InSituMeasurementsCME.cpp
        AnalyserCME.cpp
        cells/Drug.cpp
        cells/Drug.h
        factories/RateFactoryCME.cpp
        factories/RateFactoryCME.h
        Rate/ConcentrationDependentRate.cpp
        Rate/ConcentrationDependentRate.h
)

target_include_directories(simulatorCME PRIVATE ../..)
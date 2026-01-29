#!/bin/bash

#
# Copyright by Yann Bachelot
#
# Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
# https://www.leibniz-hki.de/en/applied-systems-biology.html
# HKI-Center for Systems Biology of Infection
# Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Institute (HKI)
# Adolf-Reichwein-Straße 23, 07745 Jena, Germany
#
# The project code is licensed under BSD2-Clause.
# See the LICENSE file provided with the code for the full license.
#

# first program argument is how your your new simulator is called
# you still have to add the simulator to the main.cpp and the test/src/testConfigurations.cpp

cp -r apps/example/ apps/$1

replace_word=$1

replace_word="$(tr '[:lower:]' '[:upper:]' <<< ${1:0:1})${1:1}"

find ./apps/$1 -depth -name "*Example*" -execdir bash -c 'mv "$1" "${1//Example/'"$replace_word"'}"' _ {} \;

find ./apps/$1 -depth -name "*example*" -execdir bash -c 'mv "$1" "${1//example/'"$1"'}"' _ {} \;

find ./apps/$1 -type f -print0 | xargs -0 sed -i "s/Example/$replace_word/g"

find ./apps/$1 -type f -print0 | xargs -0 sed -i "s/example/$1/g"

sed -i "/add_subdirectory(apps\/example)/a add_subdirectory(apps\/$1)" CMakeLists.txt

sed -i "/target_link_libraries(simulatorExample PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)/a target_link_libraries(simulator$replace_word PRIVATE project_options project_warnings OpenMP::OpenMP_CXX)" CMakeLists.txt

sed -i "/abm::core/a simulator$replace_word" CMakeLists.txt

file(REMOVE_RECURSE
  "libcore.pdb"
  "libcore.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/core.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()

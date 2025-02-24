# README.md

mkdir build
cd build
conan install .. --output-folder=. -g CMakeDeps -g CMakeToolchain --build=missing 
cmake .. -DCMAKE_POLICY_DEFAULT_CMP0091="NEW" -DCMAKE_TOOLCHAIN_FILE="build/generators/conan_toolchain.cmake" -G "Visual Studio 17 2022" 
cmake --build . --config Release
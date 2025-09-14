# Simple operations for modular arithmetics

This project is made to demonstrate:
1. Several elementary operations for modular arithmetic, frequently used for developing NTT
2. Template project which uses CMake for automation of work with C++ code.
3. Simple way to integrate Google Benchmark into project and evaluate performance of the mathematical operations 

# How to compile and run on Linux

```
mkdir build
cd build
cmake ..
make 
./bin/main
```

## How to compile and run on Windows

Download and install MSYS2 (https://www.msys2.org/) using default settings. Start the MSYS2 MINGW 64-bit shell and execute the following command
```
pacman -Syu
```
Run the following commands to install all pre-requisites
```
pacman -S mingw-w64-x86_64-gcc
pacman -S mingw-w64-x86_64-cmake
pacman -S make
pacman -S git
```
Clone the repo. Create a directory where the binaries will be built. The typical choice is a subfolder build. In this case, the commands are:
```
mkdir build
cd build
cmake .. -G"Unix Makefiles"
```
Build the program using the command and run an example
```
make
./bin/main.exe
```
## How to compile and run on MacOS

Install the Mac terminal command line functions if needed (type git at the command line to trigger the install). Then install home-brew if not already present:
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
```
Install pre-requisites cmake library using Homebrew:
```
brew install cmake
```
Compile and run using following commands (from the project directory)
```
mkdir build
cd build
cmake ..
make 
./bin/main
```

# How to run benchmakrs

Benchmarks are based on the Google Benchmark library. After compilation, you can run benchmarks with command from the root directory of the project:

```
./build/bin/benchmark/ntmath-bench
```

Note, that you can verify correctness operations vs naive approach by running this file:
```
./build/bin/examples/test
```

# How to run unit tests

Configure the cmake with parameter `-DBUILD_TESTING=ON` to turn the unit tests on. Unit tests are based on Google Test

```
cmake .. -DBUILD_TESTING=ON
make 
ctest
```
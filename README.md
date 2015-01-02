BCC Poisson Surface Reconstruction
=======

Requirements
----
This project should build in any *NIX environment, but requires cmake, ViennaCL, Eigen and FFTW to be installed. On OSX, these dependecies can be installed via MacPorts with the command

`sudo port install cmake eigen3 viennacl fftw3`

I would also recommend to use GCC 4.x where x >= 9, as I've only tested the code with GCC (OSX uses the Clang compiler by default, which also does not support OpenMP which I've used in some of the code to get a cheap speed increase.) To install GCC 4.9 on OSX with MacPorts, use

`sudo port install gcc49`

Build instructions
----
In the root of the project, use the commands

`mkdir build`

`cd build`

`cmake ..`

`make poisson`

If you want to build everything, including a minimalistic test suite, just run `make`. If you want to compile with GCC 4.9, the command for generating the makefile would be 

`CC=gcc-4.9 CXX=g++-4.9 cmake ..`

Running
---

Run
./poisson
to see a list of options. 

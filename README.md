BCC Poisson Surface Reconstruction
=======

Requirements
----
This project should build in any *NIX environment, but requires cmake, ViennaCL, Eigen and FFTW to be install. On OSX, these can be installed via MacPorts with the command

sudo port install cmake eigen3 viennacl fftw3


Build instructions
----

cd build
cmake ..
make poisson

Running
---

Run
./poisson
to see a list of options.

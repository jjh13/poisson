#define CATCH_CONFIG_MAIN 
#include <fstream>
#include <vector>
#include <string>

#include <catch/catch.hpp>

#include <sisl/lattice/bcc_odd.hpp>
#include <sisl/basis/quintic.hpp>
#include <poisson/ppm_writer.hpp>

#include <iostream>

using namespace sisl;
using namespace std;


TEST_CASE("PPM Writer test", "ppm test") {
	ppm_writer<> image(256, 256);

	for(int i = 0; i < 256; i++) 
		for(int j = 0; j < 256; j++)
			image.at(i,j) = vector3<unsigned char>(i,j,256-i);

	image.write("out.ppm");

}
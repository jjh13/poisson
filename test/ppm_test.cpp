#define CATCH_CONFIG_MAIN 
#include <fstream>
#include <vector>
#include <string>

#include <catch/catch.hpp>

#include <sisl/lattice/cartesian2_odd.hpp>
#include <sisl/basis/tp2linear.hpp>
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

TEST_CASE("Test cartesian odd PPM", "2d") {
	cartesian2_odd<tp2linear<double,double>, double, double> lat(1./(32.-1.));
	ppm_writer<double, double> image(256, 256);
	

	for(int x = 15; x < 25; x++) {
		for(int y = 15; y < 25; y++) {
 			lat.SV(x,y, 100.);
 		}
	}

	printf("%f", lat.f(0.5,0.5));
	lat.frequencyFilter(tp2linear<double,double>::divgradFilter);
	lat.frequencyFilter(tp2linear<double,double>::poissonFilter);

	int i=0,j=0;
	for(double x = 0; x < 1.; x+= 1./256., i++) {
		j = 0;
		for(double y = 0; y < 1.; y+= 1./256., j++) {
			unsigned char v = lat.f(x,y);

			image.at2(i,j) = vector3<double>(v,v,v);

		}
	}


	image.write2("imgtst.ppm");
}
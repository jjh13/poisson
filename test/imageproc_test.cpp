#define CATCH_CONFIG_MAIN 
#include <fstream>
#include <vector>
#include <string>

#include <catch/catch.hpp>

#include <sisl/lattice/cartesian2_odd.hpp>
#include <sisl/basis/tp2linear.hpp>
#include <sisl/basis/tp2cubic.hpp>
#include <sisl/basis/tp2quadratic.hpp>

#include <sisl/basis/dir4bs.hpp>
#include <sisl/basis/dir8bs.hpp>
#include <sisl/basis/zpelement.hpp>

#include <poisson/ppm_writer.hpp>

#include <string>
#include <iostream>
#include <tuple>

using namespace sisl;
using namespace std;

template <class GF>
void lobb_test(int res, vector<tuple<int, int, double>> qFilter, double fm) {

	// Setup the reconstruction spaces
	cartesian2_odd<GF, double, double> lat(1./(double(res)-1.));
	cartesian2_odd<GF, double, double> flat(1./(double(res)-1.));

	// Setup the output images
	ppm_writer<double, double> image(1024, 1024);
	ppm_writer<double, double> imagef(1024, 1024);

	printf("Sampling function ... \n");
	for(int x = 1; x < res - 1; x++) {
		for(int y = 1; y < res - 1; y++) {
			vector2<double> p = lat.getSitePosition(x,y);
			double xx = (p.i-0.5), yy = (p.j-0.5);
			double M = cos(2.*fm*M_PI*cos(sqrt(xx*xx + yy*yy)));

 			lat.SV(x,y, M);
 		}
	}

	for(int x = 1; x < res - 1; x++) {
		for(int y = 1; y < res - 1; y++) {
			double filtered = 0;

			for(auto it = qFilter.begin(); it != qFilter.end(); ++it) {
				tuple<int, int, double> t = *it;

				int dx = get<0>(t);
				int dy = get<1>(t);
				double w = get<2>(t);
				filtered += lat.GV(x+dx,y+dy) * w;

			}
 			flat.SV(x,y, filtered);
 		}
	}


	int i=0,j=0;
	for(double x = 0; x < 1.; x+= 1./1024., i++) {
		j = 0;
		for(double y = 0; y < 1.; y+= 1./1024., j++) {
			double v = lat.f(x,y);
			double vf = flat.f(x,y);
			image.at(i,j) = vector3<double>(v,v,v);
			imagef.at(i,j) = vector3<double>(vf,vf,vf);
		}
	}

	image.write((GF::getBasisName() + string("_") + to_string(res) + string(".ppm")));
	imagef.write((GF::getBasisName() + string("_filtered_") + to_string(res) + string(".ppm")));
}


TEST_CASE("Test 4dir boxspline  PPM", "2d") {

	double freq = 200.;

	vector<tuple<int, int, double>> qFilter4dbs = {
		make_tuple(-1,1, 1./8.),
		make_tuple(0,1, -3./8.),

		make_tuple(-1, 0, -1./4.),
		make_tuple( 0, 0, 2),

		make_tuple(1, 0, -1./4.),
		make_tuple(0, -1, -3./8.),
		make_tuple(1,-1, 1./8.),
	};


	vector<tuple<int, int, double>> qFilterLin = {
		make_tuple(0, 0, 1),
	};

	vector<tuple<int, int, double>> qFilterQuad = {
		make_tuple(0, 1, -1./2.),
		make_tuple(0,-1, -1./2.),
		make_tuple(0, 0, 4),
		make_tuple(-1, 0, -1./2.),
		make_tuple(-1, 0, -1./2.),
	};

	vector<tuple<int, int, double>> qFilterCubic = {
		make_tuple(0, 1, -1./3.),
		make_tuple(0,-1, -1./3.),
		make_tuple(0, 0, 7./3.),
		make_tuple(-1, 0, -1./3.),
		make_tuple(-1, 0, -1./3.),
	};

	vector<tuple<int, int, double>> qFilterZP= {
		make_tuple(0, 1, -1./4.),
		make_tuple(0,-1, -1./4.),
		make_tuple(0, 0, 2.),
		make_tuple(-1, 0, -1./4.),
		make_tuple(-1, 0, -1./4.),
	};

	vector<tuple<int, int, double>> qFilter6dir = {
		make_tuple(-2,-2,-1483./1800.), make_tuple(-2,-1,2653./1200.),
		make_tuple(-2,0,-1./20.), make_tuple(-2,1,2117./3600.),
		make_tuple(-1,-2,2653./1200.), make_tuple(-1,-1,-1717./360.),
		make_tuple(-1,0,-2277./400.), make_tuple(-1,2,-179./900.),
		make_tuple(0,-1,-2297./400.), make_tuple(0,0,6659./240.),
		make_tuple(0,1,-12197./1200.), make_tuple(0,2,2653./1200.),
		make_tuple(1,-2,1937./3600.), make_tuple(1,0,-12137./1200.),
		make_tuple(1,1,2027./1800.), make_tuple(2,-1,-67./450.),
		make_tuple(2,0,2593./1200.), make_tuple(2,2,-313./3600.)
	};

	for(int res = 64; res <= 512; res *= 2) {

		// Done filters
		lobb_test<dir8bs<double, double>>(res, qFilter6dir, freq); 
		lobb_test<zp_element<double, double>>(res, qFilterZP, freq); 
		lobb_test<tp2quadratic<double, double>>(res, qFilterQuad, freq); 
		lobb_test<tp2linear<double, double>>(res, qFilterLin, freq); 
		lobb_test<dir4bs<double, double>>(res, qFilter4dbs, freq); 
		lobb_test<tp2cubic<double, double>>(res, qFilterCubic, freq); 
	}
}

#define CATCH_CONFIG_MAIN 
#include <fstream>
#include <vector>
#include <string>

#include <catch/catch.hpp>

#include <sisl/lattice/cartesian2_odd.hpp>
#include <sisl/basis/tp2linear.hpp>
#include <sisl/basis/tp2cubic.hpp>

#include <sisl/basis/dir4bs.hpp>
#include <sisl/basis/zpelement.hpp>

#include <poisson/ppm_writer.hpp>

#include <iostream>
#include <tuple>

using namespace sisl;
using namespace std;


TEST_CASE("PPM Writer test", "ppm test") {
	ppm_writer<double, double> image(256, 256);

	for(int i = 0; i < 256; i++) 
		for(int j = 0; j < 256; j++)
			image.at(i,j) = vector3<double>(i,j,256-i);

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
			double v = lat.f(x,y);
			image.at(i,j) = vector3<double>(v,v,v);
		}
	}
	image.write("imgtst.ppm");
}

TEST_CASE("Test cartesian odd PPM2", "2d") {
	cartesian2_odd<zp_element<double,double>, double, double> lat(1./(32.-1.));
	ppm_writer<double, double> image(256, 256);
	

	for(int x = 15; x < 25; x++) {
		for(int y = 15; y < 25; y++) {
 			lat.SV(x,y, 100.);
 		}
	}

	int i=0,j=0;
	for(double x = 0; x < 1.; x+= 1./256., i++) {
		j = 0;
		for(double y = 0; y < 1.; y+= 1./256., j++) {
			double v = lat.f(x,y);
			image.at(i,j) = vector3<double>(v,v,v);
		}
	}
	image.write("imgtst_zp.ppm");
}

TEST_CASE("Test splatting", "") {
	typedef tp2cubic<double,double> gType;
	typedef cartesian2_odd<gType, double, double> lType;
	ppm_writer<double, double> image(256, 256);

	lType *fspace_x = new lType(1./(128.-1.));
	lType *fspace_y = new lType(1./(128.-1.));
	lType *ind  = new lType(1./(128.-1.));


	vector<vector2<double>> inputs;
	vector<double> samples_x, samples_y;

	double theta = 0.1;
	for(int i = 0; i < 63; i++) {
		double x = cos(theta*i) * 0.25 + 0.5;
		double y = sin(theta*i) * 0.25 + 0.5;

		inputs.push_back(vector2<double>(x,y));
		samples_x.push_back(cos(theta*i));
		samples_y.push_back(sin(theta*i));
	}

	//
	utility::splatSamplesOntoLattice2D<lType, gType, double, double>(
				inputs,
				samples_x,
				fspace_x,
				vector2<double>(0,0));
	
	utility::splatSamplesOntoLattice2D<lType, gType, double, double>(
				inputs,
				samples_y,
				fspace_y,
				vector2<double>(0,0));


	ind->forEachLatticeSite([&](const int & i, const int &j) { 
		double sx = 1./(fspace_x->getScale()*12.);

		double dx = (
		 -1.*fspace_x->GV(i-2, j) + 
		  8.*fspace_x->GV(i-1, j) +
		 -8.*fspace_x->GV(i+1, j) + 
			 fspace_x->GV(i+2, j)) * sx;
		
		double dy = (
		 (-1.*fspace_y->GV(i, j-2)) + 
		 ( 8.*fspace_y->GV(i, j-1)) + 
		 (-8.*fspace_y->GV(i, j+1)) + 
		 (	fspace_y->GV(i, j+2))) * sx;

		return dx+dy;
	});

	ind->frequencyFilter(tp2linear<double,double>::poissonFilter);

	int i=0,j=0;
	for(double x = 0; x < 1.; x+= 1./256., i++, j = 0) {
		for(double y = 0; y < 1.; y+= 1./256., j++) {
			double v = ind->f(x,y);
			image.at(i,j) = vector3<double>(v,v,v);
		}
	}
	image.write("splatt.ppm");

}


TEST_CASE("Test pointfit2", "") {
	typedef tp2cubic<double,double> gType;
	typedef cartesian2_odd<gType, double, double> lType;
	ppm_writer<double, double> image(256, 256);

	lType *fspace_x = new lType(1./(128.-1.));
	lType *fspace_y = new lType(1./(128.-1.));
	lType *ind  = new lType(1./(128.-1.));


	vector<vector2<double>> inputs;
	vector<double> samples_x, samples_y;

	double theta = 0.01;
	for(int i = 0; i < 2000; i++) {
		double x = cos(theta*i) * 0.25 + 0.5;
		double y = sin(theta*i) * 0.25 + 0.5;

		inputs.push_back(vector2<double>(x,y));
		samples_x.push_back(cos(theta*i));
		samples_y.push_back(sin(theta*i));
	}
	
	auto optimW = vector<tuple<int,int,double>>();
	auto bl2 = gType::getBeppoLevi2Norm();
	auto ac = gType::autoCorrelation();
	
	double lambda1 = 0.1;
	double lambda2 = 1;

	for(unsigned int i = 0; i < ac.size(); i++){
		auto b = bl2[i];
		auto a = ac[i];

		if(std::get<0>(a) != std::get<0>(b) ||
			std::get<1>(a) != std::get<1>(b))
			throw "Index mis-match while creating the optimization weight vector!";

		optimW.push_back(
			std::make_tuple(std::get<0>(a), std::get<1>(a), 
				(std::get<2>(a)*lambda1 + std::get<2>(b)*lambda2)));
	}

	utility::fitPoints2D<lType, gType, double, double> (
				inputs,
				samples_x,
				optimW,
				fspace_x,
				vector2<double>(0,0)
			);
	utility::fitPoints2D<lType, gType, double, double> (
				inputs,
				samples_y,
				optimW,
				fspace_y,
				vector2<double>(0,0)
			);

	ind->forEachLatticeSite([&](const int & i, const int &j) { 
		double sx = 1./(fspace_x->getScale()*12.);

		double dx = (
		 -1.*fspace_x->GV(i-2, j) + 
		  8.*fspace_x->GV(i-1, j) +
		 -8.*fspace_x->GV(i+1, j) + 
			 fspace_x->GV(i+2, j)) * sx;
		
		double dy = (
		 (-1.*fspace_y->GV(i, j-2)) + 
		 ( 8.*fspace_y->GV(i, j-1)) + 
		 (-8.*fspace_y->GV(i, j+1)) + 
		 (	fspace_y->GV(i, j+2))) * sx;

		return dx+dy;
	});
 
	ind->frequencyFilter(tp2linear<double,double>::poissonFilter);

	int i=0,j=0;
	for(double x = 0; x < 1.; x+= 1./256., i++, j = 0) {
		for(double y = 0; y < 1.; y+= 1./256., j++) {
			double v = ind->f(x,y);
			image.at(i,j) = vector3<double>(v,v,v);
		}
	}
	image.write("pfit.ppm");

}

TEST_CASE("Test 4dir boxspline  PPM", "2d") {
#define RES 128
	cartesian2_odd<zp_element<double,double>, double, double> lat(1./(double(RES)-1.));
	cartesian2_odd<zp_element<double,double>, double, double> flat(1./(double(RES)-1.));
	ppm_writer<double, double> image(1024, 1024);
	ppm_writer<double, double> imagef(1024, 1024);

	vector<tuple<int, int, double>> qFilter = {
		make_tuple(-1,1, 1./8.),
		make_tuple(0,1, -3./8.),

		make_tuple(-1, 0, -1./4.),
		make_tuple( 0, 0, 2),

		make_tuple(1, 0, -1./4.),
		make_tuple(0, -1, -3./8.),
		make_tuple(1,-1, 1./8.),
	};

	
	printf("Sampling function ... \n");

	for(int x = 1; x < RES - 1; x++) {
		for(int y = 1; y < RES - 1; y++) {
			vector2<double> p = lat.getSitePosition(x,y);
			double xx = (p.i-0.5), yy = (p.j-0.5);
			double M = cos(160.*M_PI*cos(sqrt(xx*xx + yy*yy)));

 			lat.SV(x,y, M);
 		}
	}

	printf("Filtering function ... \n");

	for(int x = 1; x < RES - 1; x++) {
		for(int y = 1; y < RES - 1; y++) {
			double filtered =  lat.GV(x,y);

			// for(auto it = qFilter.begin(); it != qFilter.end(); ++it) {
			// 	tuple<int, int, double> t = *it;

			// 	int dx = get<0>(t);
			// 	int dy = get<1>(t);
			// 	double w = get<2>(t);
			// 	filtered += lat.GV(x+dx,y+dy) * w;

			// }


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
	image.write("img4dir.ppm");
	imagef.write("img4dir_filtered.ppm");
}

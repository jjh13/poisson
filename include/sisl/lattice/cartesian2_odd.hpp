
#ifndef _CARTESIAN2_ODD_H
#define _CARTESIAN2_ODD_H

#include <sisl/array.hpp>
#include <sisl/sparse_array.hpp>
#include <sisl/lattice.hpp>
#include <sisl/basis_function.hpp>
#include <sisl/fftwalloc.hpp>
#include <fftw3.h>

#include <iostream>

namespace sisl{
template <class GF, class O = float, class I = float>
class cartesian2_odd : public shift_invariant_space2<GF,O,I> {
public:
	cartesian2_odd() : cartesian2_odd(1) {} // should initialize with h = 1
	cartesian2_odd(const I &h) : 
		dh(h), 
		res((1./h)+1),
		NA(res-2,res-2,vector2<O>(0,0)),
		L(res-2,res-2) { };
	~cartesian2_odd(){}

	virtual O GV(const int &x, const int &y) const {
		int xr = MODP(x, res-1);
		int yr = MODP(y, res-1);

		if(xr == 0 || yr == 0) 
			return 0;
		 return L(xr-1, yr-1) *
		 			(x < 0 || x > (res - 1) ? -1 : 1) * 
		 			(y < 0 || y > (res - 1) ? -1 : 1);
	}

	virtual void SV(const int &x, const int &y, const O &value) {
		int xr = MODP(x, (res - 1));
		int yr = MODP(y, (res - 1));

		if(xr == 0 || yr == 0) 
			return;

		L(xr-1, yr-1) = value;
	}
	virtual vector2<O> GN(const int &x, const int &y) {
		int xr = MODP(x, res-1);
		int yr = MODP(y, res-1);
		
		if(xr == 0 || yr == 0) 
			return vector2<O>(0,0);

		//if(!NA.hasValue(xr-1,yr-1))
		// 	NA(xr-1,yr-1,zr-1) = GF::approximateGradient(x,y,z, this);

		// return NA(xr-1, yr-1, zr-1) *
		// 			(x < 0 || x > (res - 1) ? -1 : 1) * 
		// 			(y < 0 || y > (res - 1) ? -1 : 1) * 
		// 			(z < 0 || z > (res - 1) ? -1 : 1);
		return vector2<O>(0,0);
	}

	virtual void SN(const int &x, const int &y, const vector2<O> &value) {
		int xr = MODP(x, (res - 1));
		int yr = MODP(y, (res - 1));
		
		if(xr == 0 || yr == 0) 
			return;

		NA(xr-1, yr-1) = value;
	}
	virtual const int lIndex(const int &x, const int &y) const { 
		int xr = MODP(x, (res - 1));
		int yr = MODP(y, (res - 1));

		if(xr == 0 || yr == 0 || x >=  (res - 1) || y >=  (res - 1)  || x < 0 || y < 0) 
			return -1;

		return L.lIndex(xr-1, yr-1);
	}
	virtual const int numberOfLatticeSites() const {
		return L.size();
	}

	virtual I getScale() const { 
		return dh;
	}
	virtual const int getResolution() const {
		return res-2;
	}

	virtual vector2<I> getSitePosition(const int &x, const int &y) const { 
		return vector2<I>(I(x)*dh,I(y)*dh); 
	}
	virtual vector2<int> getNearestIndex(const vector2<I> &pos) const {
		return vector2<int>((int)floor(pos.i/dh), (int)floor(pos.j/dh));
	}

	// Evaluates the function  
	virtual const O f(const I &x, const I &y) { return GF::convolutionSum(vector2<I>(x,y), this); }
	virtual vector2<O> grad_f(const I &x, const I &y) { return vector2<O>(0,0); } //return GF::convolutionSumNormal(vector3<I>(x,y,z), this); }

	// Evaluates the function  
	virtual const O f(const vector2<I>& p ) { return  GF::convolutionSum(p, this); }
	virtual vector2<O> grad_f(const vector2<I>& p) { return vector2<O>(0,0); } //GF::convolutionSumNormal(p, this); }

	//
	virtual void forEachLatticeSite(std::function<O(const int &, const int &)> lambda) {
		#pragma omp parallel for
		for(int i = 1; i < this->res - 1; i++)
			for(int j = 1; j < this->res - 1; j++){
					SV(i,j, lambda(i,j));
				}
	}
	
	//
	virtual void spatialFilter(std::function<O(const int &, const int &, std::function<O(const int &, const int & )>)> lambda) {
		// cartesian_odd <GF, O, I> src = *this; // Copy the lattice

		// std::function<O(const int &, const int &, const int &)> accessFunc = 
		// 	[&](const int &x, const int &y, const int &z) ->  O {
		// 		return (O) src.GV(x,y,z); 
		// 	};
		
		// #pragma omp parallel for
		// for(int i = 1; i < this->res - 1; i++)
		// 	for(int j = 1; i < this->res - 1; j++)
		// 		for(int k = 1; i < this->res - 1; k++) {
		// 			SV(i,j,k, lambda(i,j,k, accessFunc));
		// 		}
	} 
	//
	virtual void frequencyFilter(std::function<O(const I &u, const I &v, const I &h)> lambda) {
		I scale =  (0.25*dh*dh);
		I invdh = 1./(dh*dh);
		O *lattice = L.getArray();

		// Do the forward transform
		fftw_plan p = fftw_plan_r2r_2d(res-2, res-2,
					lattice, lattice, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);
		fftw_execute(p);
		fftw_destroy_plan(p);

		// // Apply the filter(s)
		#pragma omp parallel for
		for(int i = 0; i<res-2; i++) 
			for(int j = 0; j<res-2; j++) {
					I u = M_PI*I(i+1)*dh;
					I v = M_PI*I(j+1)*dh;
					L(i,j) *= lambda(u,v,dh) * scale;
				}
		
		// Do the inverse transform
		p = fftw_plan_r2r_2d(res-2, res-2, 
					lattice, lattice, FFTW_RODFT00, FFTW_RODFT00, FFTW_ESTIMATE);

		fftw_execute(p);
		fftw_destroy_plan(p);
	} 

	// Convolution of two functions
	//inline cartesian_odd<GF, O, I> operator*(const cartesian_odd<GF, O, I> &v) {throw;};

	// Copy a function
	inline cartesian2_odd<GF, O, I> operator=(const cartesian2_odd<GF, O, I> &v) {throw;};

	virtual std::string getLatticeName() const {
		return std::string("cartesian2");
	}
	virtual std::string getBasisName() const {
		return GF::getBasisName();
	}
private:
	I dh;
	int res;
	array2<O, _fftwalloc<O> > L;
	sparse_array2<vector2<O> > NA;
};
};

#endif // _CARTESIAN_H

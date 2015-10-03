/**
 * Basic layout for a generating function on 
 * a lattice
 * 
 * @author Joshua Horacsek
 */

#ifndef _BASIS_FUNCTION_H_
#define _BASIS_FUNCTION_H_

#include <vector>
#include <string>
#include <tuple>
#include <sisl/primitives.hpp>
#include <sisl/lattice.hpp>

namespace sisl{


template <class O, class I>
class basis_function {
public:
	static const O M(const I &x, const I &y, const I &z) { return (O)0; };
	static const O M(const vector3<I> &p) { return (O)0; };

	// This function should return the intersection of the closure of the support of
	// the generator and the lattice.
	static std::vector<std::tuple<int,int,int>> getSupport() {throw "basis_function()::getSupport() - Not Implemented!";};

	// Gets the lattices sites that actually contribute to point p
	// p should be in the vornoi cell of the 0 element of the lattice,
	// if not, the behaviour of this function is undefined.
	static std::vector<std::tuple<int,int,int>> getEffectiveSupport(const vector3<I> &p) {throw "basis_function()::getEffectiveSupport() - Not Implemented!";};

	// Gets the lattices sites that actually contribute to point p
	// p should be in the vornoi cell of the 0 element of the lattice,
	// if not, the behaviour of this function is undefined.
	static std::vector<std::tuple<int,int,int,O>> getBeppoLevi2Norm(){throw "basis_function()::getBeppoLevi2Norm() - Not Implemented!";};
	static std::vector<std::tuple<int,int,int,O>> getBeppoLevi1Norm(){throw "basis_function()::getBeppoLevi1Norm() - Not Implemented!";};
	static std::vector<std::tuple<int,int,int,O>> autoCorrelation(){throw "basis_function()::autoCorrelation() - Not Implemented!";};

	// This does the actual semi-descrete convolution sum
	// the idea behind this, is to allow basis functions
	// to provide a potentially optimized version of
	// the convolution sum for particular lattices.
	static const O convolutionSum(const vector3<I> &p, const shift_invariant_space3< basis_function, O, I> *lattice) {throw "basis_function()::convolutionSum() - Not Implemented!";};

	//
	static std::string getBasisName(){ throw "basis_function::getBasisName not implemented!";}
};

template <class O, class I>
class basis_function2 {
public:
	static const O M(const I &x, const I &y) { return (O)0; };
	static const O M(const vector2<I> &p) { return (O)0; };

	// This function should return the intersection of the closure of the support of
	// the generator and the lattice.
	static std::vector<std::tuple<int,int>> getSupport() {throw "basis_function()::getSupport() - Not Implemented!";};

	// Gets the lattices sites that actually contribute to point p
	// p should be in the vornoi cell of the 0 element of the lattice,
	// if not, the behaviour of this function is undefined.
	static std::vector<std::tuple<int,int>> getEffectiveSupport(const vector3<I> &p) {throw "basis_function()::getEffectiveSupport() - Not Implemented!";};

	// Gets the lattices sites that actually contribute to point p
	// p should be in the vornoi cell of the 0 element of the lattice,
	// if not, the behaviour of this function is undefined.
	static std::vector<std::tuple<int,int,O>> getBeppoLevi2Norm(){throw "basis_function()::getBeppoLevi2Norm() - Not Implemented!";};
	static std::vector<std::tuple<int,int,O>> getBeppoLevi1Norm(){throw "basis_function()::getBeppoLevi1Norm() - Not Implemented!";};
	static std::vector<std::tuple<int,int,O>> autoCorrelation(){throw "basis_function()::autoCorrelation() - Not Implemented!";};

	// This does the actual semi-descrete convolution sum
	// the idea behind this, is to allow basis functions
	// to provide a potentially optimized version of
	// the convolution sum for particular lattices.
	static const O convolutionSum(const vector2<I> &p, const shift_invariant_space2< basis_function2, O, I> *lattice) {throw "basis_function()::convolutionSum() - Not Implemented!";};

	//
	static std::string getBasisName(){ throw "basis_function::getBasisName not implemented!";}
};

};
#endif // _BASIS_FUNCTION_H_

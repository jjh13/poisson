
#include <math.h>
#include <sisl/basis_function.hpp>
#include <sisl/primitives.hpp>
#include <sisl/lattice.hpp>
#include <vector>

#ifndef _TP_2D_8DBS_H_
#define _TP_2D_8DBS_H_

namespace sisl{
template <class O, class I>
class dir8bs : public basis_function <O,I> {


public:
    static std::string getBasisName(){
        return std::string("dir8bs");
    }
    static const O M(const I &x, const I &y) { return (O) box_spline(x,y); };
    static const O M(const vector2<I> &p) { return (O) box_spline(p.i, p.j); };

    // This function should return the intersection of the closure of the support of
    // the generator and the lattice.
    static std::vector<std::tuple<int,int>> getSupport() {
        using namespace std;
        std::vector<std::tuple<int,int>> support;
        for(int i = -4; i <= 4; i++)
            for(int j = -4; j <= 4; j++) {
                support.push_back(make_tuple(i,j));
            }
        return support;
    };

    static std::vector<std::tuple<int,int>> getEffectiveSupport(const vector2<I> &p) {
        using namespace std;
        std::vector<std::tuple<int,int>> support;
        for(int i = -4; i <= 4; i++)
            for(int j = -4; j <= 4; j++) {
                    support.push_back(make_tuple(i,j));
                }
        return support;
    }

    static std::vector<std::tuple<int,int,O>> getBeppoLevi2Norm() {throw "basis_function()::getBeppoLevi2Norm() - Not Implemented!";};
    static std::vector<std::tuple<int,int,O>> getBeppoLevi1Norm(){throw "basis_function()::getBeppoLevi1Norm() - Not Implemented!";};
    static std::vector<std::tuple<int,int,O>> autoCorrelation(){throw "basis_function()::autoCorrelation() - Not Implemented!";};

    inline static const O convolutionSum(const vector2<I> &p, const shift_invariant_space2< dir8bs, O, I> *lattice) {
        O sum = 0;
        I dh = lattice->getScale();
        vector2<I> vox(p.i/dh, p.j/dh);     

        int vx = (int)floor(vox.i), 
            vy = (int)floor(vox.j);

        for(int i = -4; i <= 4; i++)
            for(int j = -4; j <= 4; j++) {
                sum += lattice->GV(vx + i, vy + j) * M(vox.i - I(vx + i), vox.j - I(vy + j));
            }
        return sum;     
    }

private: 
    static inline O __pp_r1__(const I &x_0, const I &x_1) {
       O __pp_r1___result;
       __pp_r1___result = (1.0L/960.0L)*pow(x_0, 6) - 1.0L/240.0L*pow(x_0, 5)*x_1 - 23.0L/960.0L*pow(x_0, 5) + (7.0L/96.0L)*pow(x_0, 4)*x_1 + (175.0L/768.0L)*pow(x_0, 4) - 49.0L/96.0L*pow(x_0, 3)*x_1 - 147.0L/128.0L*pow(x_0, 3) + (343.0L/192.0L)*pow(x_0, 2)*x_1 + (9947.0L/3072.0L)*pow(x_0, 2) - 2401.0L/768.0L*x_0*x_1 - 74431.0L/15360.0L*x_0 + (16807.0L/7680.0L)*x_1 + 184877.0L/61440.0L;
       return __pp_r1___result;
    }
    static inline O __pp_r2__(const I &x_0, const I &x_1) {
       O __pp_r2___result;
       __pp_r2___result = (1.0L/960.0L)*pow(x_0, 6) - 1.0L/240.0L*pow(x_0, 5)*x_1 - 23.0L/960.0L*pow(x_0, 5) + (7.0L/96.0L)*pow(x_0, 4)*x_1 + (175.0L/768.0L)*pow(x_0, 4) - 49.0L/96.0L*pow(x_0, 3)*x_1 - 147.0L/128.0L*pow(x_0, 3) + (343.0L/192.0L)*pow(x_0, 2)*x_1 + (9947.0L/3072.0L)*pow(x_0, 2) - 2401.0L/768.0L*x_0*x_1 - 74431.0L/15360.0L*x_0 + (16807.0L/7680.0L)*x_1 + 184877.0L/61440.0L;
       return __pp_r2___result;
    }
    static inline O __pp_r3__(const I &x_0, const I &x_1) {
       O __pp_r3___result;
       __pp_r3___result = (11.0L/8640.0L)*pow(x_0, 6) - 1.0L/360.0L*pow(x_0, 5)*x_1 - 9.0L/320.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/96.0L)*pow(x_0, 4)*x_1 + (199.0L/768.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) - 37.0L/96.0L*pow(x_0, 3)*x_1 - 163.0L/128.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (3.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (271.0L/192.0L)*pow(x_0, 2)*x_1 + (10811.0L/3072.0L)*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 1.0L/48.0L*x_0*pow(x_1, 4) + (1.0L/8.0L)*x_0*pow(x_1, 3) - 3.0L/8.0L*x_0*pow(x_1, 2) - 1969.0L/768.0L*x_0*x_1 - 15923.0L/3072.0L*x_0 + (1.0L/4320.0L)*pow(x_1, 6) - 1.0L/240.0L*pow(x_1, 5) + (1.0L/32.0L)*pow(x_1, 4) - 1.0L/8.0L*pow(x_1, 3) + (9.0L/32.0L)*pow(x_1, 2) + (2843.0L/1536.0L)*x_1 + 39049.0L/12288.0L;
       return __pp_r3___result;
    }
    static inline O __pp_r4__(const I &x_0, const I &x_1) {
       O __pp_r4___result;
       __pp_r4___result = -83.0L/43200.0L*pow(x_0, 6) + (17.0L/3600.0L)*pow(x_0, 5)*x_1 + (487.0L/14400.0L)*pow(x_0, 5) - 1.0L/90.0L*pow(x_0, 4)*pow(x_1, 2) - 103.0L/1440.0L*pow(x_0, 4)*x_1 - 2783.0L/11520.0L*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) + (13.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (617.0L/1440.0L)*pow(x_0, 3)*x_1 + (15307.0L/17280.0L)*pow(x_0, 3) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 4) - 13.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 169.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) - 3643.0L/2880.0L*pow(x_0, 2)*x_1 - 79283.0L/46080.0L*pow(x_0, 2) + (1.0L/1800.0L)*x_0*pow(x_1, 5) + (13.0L/720.0L)*x_0*pow(x_1, 4) + (169.0L/720.0L)*x_0*pow(x_1, 3) + (2197.0L/1440.0L)*x_0*pow(x_1, 2) + (21107.0L/11520.0L)*x_0*x_1 + (368707.0L/230400.0L)*x_0 - 1.0L/21600.0L*pow(x_1, 6) - 13.0L/7200.0L*pow(x_1, 5) - 169.0L/5760.0L*pow(x_1, 4) - 2197.0L/8640.0L*pow(x_1, 3) - 28561.0L/23040.0L*pow(x_1, 2) - 29797.0L/28800.0L*x_1 - 1334153.0L/2764800.0L;
       return __pp_r4___result;
    }
    static inline O __pp_r5__(const I &x_0, const I &x_1) {
       O __pp_r5___result;
       __pp_r5___result = -83.0L/43200.0L*pow(x_0, 6) + (17.0L/3600.0L)*pow(x_0, 5)*x_1 + (487.0L/14400.0L)*pow(x_0, 5) - 1.0L/90.0L*pow(x_0, 4)*pow(x_1, 2) - 103.0L/1440.0L*pow(x_0, 4)*x_1 - 2783.0L/11520.0L*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) + (13.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (617.0L/1440.0L)*pow(x_0, 3)*x_1 + (15307.0L/17280.0L)*pow(x_0, 3) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 4) - 13.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 169.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) - 3643.0L/2880.0L*pow(x_0, 2)*x_1 - 79283.0L/46080.0L*pow(x_0, 2) + (1.0L/1800.0L)*x_0*pow(x_1, 5) + (13.0L/720.0L)*x_0*pow(x_1, 4) + (169.0L/720.0L)*x_0*pow(x_1, 3) + (2197.0L/1440.0L)*x_0*pow(x_1, 2) + (21107.0L/11520.0L)*x_0*x_1 + (368707.0L/230400.0L)*x_0 - 1.0L/21600.0L*pow(x_1, 6) - 13.0L/7200.0L*pow(x_1, 5) - 169.0L/5760.0L*pow(x_1, 4) - 2197.0L/8640.0L*pow(x_1, 3) - 28561.0L/23040.0L*pow(x_1, 2) - 29797.0L/28800.0L*x_1 - 1334153.0L/2764800.0L;
       return __pp_r5___result;
    }
    static inline O __pp_r6__(const I &x_0, const I &x_1) {
       O __pp_r6___result;
       __pp_r6___result = -73.0L/43200.0L*pow(x_0, 6) + (11.0L/1800.0L)*pow(x_0, 5)*x_1 + (427.0L/14400.0L)*pow(x_0, 5) - 11.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 133.0L/1440.0L*pow(x_0, 4)*x_1 - 2423.0L/11520.0L*pow(x_0, 4) + (13.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (37.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (797.0L/1440.0L)*pow(x_0, 3)*x_1 + (13147.0L/17280.0L)*pow(x_0, 3) + (1.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 41.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 31.0L/60.0L*pow(x_0, 2)*pow(x_1, 2) - 4723.0L/2880.0L*pow(x_0, 2)*x_1 - 66323.0L/46080.0L*pow(x_0, 2) + (7.0L/3600.0L)*x_0*pow(x_1, 5) - 1.0L/360.0L*x_0*pow(x_1, 4) + (259.0L/720.0L)*x_0*pow(x_1, 3) + (1657.0L/1440.0L)*x_0*pow(x_1, 2) + (27587.0L/11520.0L)*x_0*x_1 + (290947.0L/230400.0L)*x_0 + (1.0L/5400.0L)*pow(x_1, 6) - 43.0L/7200.0L*pow(x_1, 5) + (11.0L/5760.0L)*pow(x_1, 4) - 3277.0L/8640.0L*pow(x_1, 3) - 22081.0L/23040.0L*pow(x_1, 2) - 39517.0L/28800.0L*x_1 - 867593.0L/2764800.0L;
       return __pp_r6___result;
    }
    static inline O __pp_r7__(const I &x_0, const I &x_1) {
       O __pp_r7___result;
       __pp_r7___result = -7.0L/4800.0L*pow(x_0, 6) + (7.0L/3600.0L)*pow(x_0, 5)*x_1 + (367.0L/14400.0L)*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) - 43.0L/1440.0L*pow(x_0, 4)*x_1 - 2063.0L/11520.0L*pow(x_0, 4) - 1.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) + (11.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (257.0L/1440.0L)*pow(x_0, 3)*x_1 + (10987.0L/17280.0L)*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) - 79.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) - 1483.0L/2880.0L*pow(x_0, 2)*x_1 - 53363.0L/46080.0L*pow(x_0, 2) - 1.0L/450.0L*x_0*pow(x_1, 5) - 17.0L/720.0L*x_0*pow(x_1, 4) - 11.0L/720.0L*x_0*pow(x_1, 3) + (1117.0L/1440.0L)*x_0*pow(x_1, 2) + (8147.0L/11520.0L)*x_0*x_1 + (213187.0L/230400.0L)*x_0 + (1.0L/2400.0L)*pow(x_1, 6) + (47.0L/7200.0L)*pow(x_1, 5) + (191.0L/5760.0L)*pow(x_1, 4) - 37.0L/8640.0L*pow(x_1, 3) - 15601.0L/23040.0L*pow(x_1, 2) - 10357.0L/28800.0L*x_1 - 401033.0L/2764800.0L;
       return __pp_r7___result;
    }
    static inline O __pp_r8__(const I &x_0, const I &x_1) {
       O __pp_r8___result;
       __pp_r8___result = -7.0L/4800.0L*pow(x_0, 6) + (7.0L/3600.0L)*pow(x_0, 5)*x_1 + (367.0L/14400.0L)*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) - 43.0L/1440.0L*pow(x_0, 4)*x_1 - 2063.0L/11520.0L*pow(x_0, 4) - 1.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) + (11.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (257.0L/1440.0L)*pow(x_0, 3)*x_1 + (10987.0L/17280.0L)*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) - 79.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) - 1483.0L/2880.0L*pow(x_0, 2)*x_1 - 53363.0L/46080.0L*pow(x_0, 2) - 1.0L/450.0L*x_0*pow(x_1, 5) - 17.0L/720.0L*x_0*pow(x_1, 4) - 11.0L/720.0L*x_0*pow(x_1, 3) + (1117.0L/1440.0L)*x_0*pow(x_1, 2) + (8147.0L/11520.0L)*x_0*x_1 + (213187.0L/230400.0L)*x_0 + (1.0L/2400.0L)*pow(x_1, 6) + (47.0L/7200.0L)*pow(x_1, 5) + (191.0L/5760.0L)*pow(x_1, 4) - 37.0L/8640.0L*pow(x_1, 3) - 15601.0L/23040.0L*pow(x_1, 2) - 10357.0L/28800.0L*x_1 - 401033.0L/2764800.0L;
       return __pp_r8___result;
    }
    static inline O __pp_r9__(const I &x_0, const I &x_1) {
       O __pp_r9___result;
       __pp_r9___result = -53.0L/43200.0L*pow(x_0, 6) + (1.0L/300.0L)*pow(x_0, 5)*x_1 + (307.0L/14400.0L)*pow(x_0, 5) - 1.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 73.0L/1440.0L*pow(x_0, 4)*x_1 - 1703.0L/11520.0L*pow(x_0, 4) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (437.0L/1440.0L)*pow(x_0, 3)*x_1 + (8827.0L/17280.0L)*pow(x_0, 3) + (11.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 11.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 17.0L/120.0L*pow(x_0, 2)*pow(x_1, 2) - 2563.0L/2880.0L*pow(x_0, 2)*x_1 - 40403.0L/46080.0L*pow(x_0, 2) - 1.0L/1200.0L*x_0*pow(x_1, 5) - 2.0L/45.0L*x_0*pow(x_1, 4) + (79.0L/720.0L)*x_0*pow(x_1, 3) + (577.0L/1440.0L)*x_0*pow(x_1, 2) + (14627.0L/11520.0L)*x_0*x_1 + (135427.0L/230400.0L)*x_0 + (7.0L/10800.0L)*pow(x_1, 6) + (17.0L/7200.0L)*pow(x_1, 5) + (371.0L/5760.0L)*pow(x_1, 4) - 1117.0L/8640.0L*pow(x_1, 3) - 9121.0L/23040.0L*pow(x_1, 2) - 20077.0L/28800.0L*x_1 + 65527.0L/2764800.0L;
       return __pp_r9___result;
    }
    static inline O __pp_r10__(const I &x_0, const I &x_1) {
       O __pp_r10___result;
       __pp_r10___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 + (1.0L/960.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/32.0L)*pow(x_0, 4)*x_1 - 127.0L/11520.0L*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) - 227.0L/1440.0L*pow(x_0, 3)*x_1 + (113.0L/5760.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) - 37.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (131.0L/320.0L)*pow(x_0, 2)*x_1 + (5201.0L/46080.0L)*pow(x_0, 2) - 1.0L/400.0L*x_0*pow(x_1, 5) - 1.0L/32.0L*x_0*pow(x_1, 4) - 143.0L/1440.0L*x_0*pow(x_1, 3) + (301.0L/960.0L)*x_0*pow(x_1, 2) - 3247.0L/5760.0L*x_0*x_1 - 7261.0L/15360.0L*x_0 + (19.0L/43200.0L)*pow(x_1, 6) + (7.0L/960.0L)*pow(x_1, 5) + (503.0L/11520.0L)*pow(x_1, 4) + (419.0L/5760.0L)*pow(x_1, 3) - 16561.0L/46080.0L*pow(x_1, 2) + (5213.0L/15360.0L)*x_1 + 42829.0L/86400.0L;
       return __pp_r10___result;
    }
    static inline O __pp_r11__(const I &x_0, const I &x_1) {
       O __pp_r11___result;
       __pp_r11___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 + (1.0L/960.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/32.0L)*pow(x_0, 4)*x_1 - 127.0L/11520.0L*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) - 227.0L/1440.0L*pow(x_0, 3)*x_1 + (113.0L/5760.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) - 37.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (131.0L/320.0L)*pow(x_0, 2)*x_1 + (5201.0L/46080.0L)*pow(x_0, 2) - 1.0L/400.0L*x_0*pow(x_1, 5) - 1.0L/32.0L*x_0*pow(x_1, 4) - 143.0L/1440.0L*x_0*pow(x_1, 3) + (301.0L/960.0L)*x_0*pow(x_1, 2) - 3247.0L/5760.0L*x_0*x_1 - 7261.0L/15360.0L*x_0 + (19.0L/43200.0L)*pow(x_1, 6) + (7.0L/960.0L)*pow(x_1, 5) + (503.0L/11520.0L)*pow(x_1, 4) + (419.0L/5760.0L)*pow(x_1, 3) - 16561.0L/46080.0L*pow(x_1, 2) + (5213.0L/15360.0L)*x_1 + 42829.0L/86400.0L;
       return __pp_r11___result;
    }
    static inline O __pp_r12__(const I &x_0, const I &x_1) {
       O __pp_r12___result;
       __pp_r12___result = (11.0L/43200.0L)*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 1.0L/320.0L*pow(x_0, 5) + (7.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 4)*x_1 + (233.0L/11520.0L)*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) - 47.0L/1440.0L*pow(x_0, 3)*x_1 - 607.0L/5760.0L*pow(x_0, 3) + (13.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (53.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (11.0L/320.0L)*pow(x_0, 2)*x_1 + (18161.0L/46080.0L)*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) - 5.0L/96.0L*x_0*pow(x_1, 4) + (37.0L/1440.0L)*x_0*pow(x_1, 3) - 59.0L/960.0L*x_0*pow(x_1, 2) - 7.0L/5760.0L*x_0*x_1 - 2489.0L/3072.0L*x_0 + (29.0L/43200.0L)*pow(x_1, 6) + (1.0L/320.0L)*pow(x_1, 5) + (863.0L/11520.0L)*pow(x_1, 4) - 301.0L/5760.0L*pow(x_1, 3) - 3601.0L/46080.0L*pow(x_1, 2) + (29.0L/15360.0L)*x_1 + 57409.0L/86400.0L;
       return __pp_r12___result;
    }
    static inline O __pp_r13__(const I &x_0, const I &x_1) {
       O __pp_r13___result;
       __pp_r13___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 + (1.0L/64.0L)*pow(x_0, 5) - 7.0L/96.0L*pow(x_0, 4)*x_1 - 21.0L/256.0L*pow(x_0, 4) + (49.0L/96.0L)*pow(x_0, 3)*x_1 + (49.0L/384.0L)*pow(x_0, 3) - 343.0L/192.0L*pow(x_0, 2)*x_1 + (343.0L/1024.0L)*pow(x_0, 2) + (2401.0L/768.0L)*x_0*x_1 - 7203.0L/5120.0L*x_0 - 16807.0L/7680.0L*x_1 + 16807.0L/12288.0L;
       return __pp_r13___result;
    }
    static inline O __pp_r14__(const I &x_0, const I &x_1) {
       O __pp_r14___result;
       __pp_r14___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 + (1.0L/64.0L)*pow(x_0, 5) - 7.0L/96.0L*pow(x_0, 4)*x_1 - 21.0L/256.0L*pow(x_0, 4) + (49.0L/96.0L)*pow(x_0, 3)*x_1 + (49.0L/384.0L)*pow(x_0, 3) - 343.0L/192.0L*pow(x_0, 2)*x_1 + (343.0L/1024.0L)*pow(x_0, 2) + (2401.0L/768.0L)*x_0*x_1 - 7203.0L/5120.0L*x_0 - 16807.0L/7680.0L*x_1 + 16807.0L/12288.0L;
       return __pp_r14___result;
    }
    static inline O __pp_r15__(const I &x_0, const I &x_1) {
       O __pp_r15___result;
       __pp_r15___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 + (77.0L/2880.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 5.0L/288.0L*pow(x_0, 4)*x_1 - 445.0L/2304.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/9.0L)*pow(x_0, 3)*pow(x_1, 2) + (19.0L/288.0L)*pow(x_0, 3)*x_1 + (2489.0L/3456.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/9.0L)*pow(x_0, 2)*pow(x_1, 3) - 2.0L/3.0L*pow(x_0, 2)*pow(x_1, 2) - 5.0L/576.0L*pow(x_0, 2)*x_1 - 13297.0L/9216.0L*pow(x_0, 2) - 1.0L/360.0L*x_0*pow(x_1, 5) + (1.0L/18.0L)*x_0*pow(x_1, 4) - 4.0L/9.0L*x_0*pow(x_1, 3) + (16.0L/9.0L)*x_0*pow(x_1, 2) - 989.0L/2304.0L*x_0*x_1 + (13249.0L/9216.0L)*x_0 - 1.0L/2160.0L*pow(x_1, 6) + (1.0L/90.0L)*pow(x_1, 5) - 1.0L/9.0L*pow(x_1, 4) + (16.0L/27.0L)*pow(x_1, 3) - 16.0L/9.0L*pow(x_1, 2) + (3023.0L/4608.0L)*x_1 - 292261.0L/552960.0L;
       return __pp_r15___result;
    }
    static inline O __pp_r16__(const I &x_0, const I &x_1) {
       O __pp_r16___result;
       __pp_r16___result = -7.0L/4800.0L*pow(x_0, 6) + (7.0L/3600.0L)*pow(x_0, 5)*x_1 + (367.0L/14400.0L)*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) - 43.0L/1440.0L*pow(x_0, 4)*x_1 - 2063.0L/11520.0L*pow(x_0, 4) - 1.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) + (11.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (257.0L/1440.0L)*pow(x_0, 3)*x_1 + (10987.0L/17280.0L)*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) - 79.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) - 1483.0L/2880.0L*pow(x_0, 2)*x_1 - 53363.0L/46080.0L*pow(x_0, 2) + (11.0L/1800.0L)*x_0*pow(x_1, 5) - 2.0L/45.0L*x_0*pow(x_1, 4) + (1.0L/180.0L)*x_0*pow(x_1, 3) + (551.0L/720.0L)*x_0*pow(x_1, 2) + (8177.0L/11520.0L)*x_0*x_1 + (213127.0L/230400.0L)*x_0 + (1.0L/400.0L)*pow(x_1, 6) - 13.0L/450.0L*pow(x_1, 5) + (41.0L/360.0L)*pow(x_1, 4) - 89.0L/1080.0L*pow(x_1, 3) - 3679.0L/5760.0L*pow(x_1, 2) - 42523.0L/115200.0L*x_1 - 398423.0L/2764800.0L;
       return __pp_r16___result;
    }
    static inline O __pp_r17__(const I &x_0, const I &x_1) {
       O __pp_r17___result;
       __pp_r17___result = (19.0L/43200.0L)*pow(x_0, 6) - 1.0L/3600.0L*pow(x_0, 5)*x_1 - 127.0L/14400.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 17.0L/1440.0L*pow(x_0, 4)*x_1 + (991.0L/11520.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 11.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (251.0L/1440.0L)*pow(x_0, 3)*x_1 - 8443.0L/17280.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (11.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (121.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 2483.0L/2880.0L*pow(x_0, 2)*x_1 + (73999.0L/46080.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 11.0L/1440.0L*x_0*pow(x_1, 4) - 121.0L/1440.0L*x_0*pow(x_1, 3) - 1331.0L/2880.0L*x_0*pow(x_1, 2) + (10687.0L/5760.0L)*x_0*x_1 - 646237.0L/230400.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (11.0L/14400.0L)*pow(x_1, 5) + (121.0L/11520.0L)*pow(x_1, 4) + (1331.0L/17280.0L)*pow(x_1, 3) + (14641.0L/46080.0L)*pow(x_1, 2) - 343159.0L/230400.0L*x_1 + 347071.0L/172800.0L;
       return __pp_r17___result;
    }
    static inline O __pp_r18__(const I &x_0, const I &x_1) {
       O __pp_r18___result;
       __pp_r18___result = (19.0L/43200.0L)*pow(x_0, 6) - 1.0L/3600.0L*pow(x_0, 5)*x_1 - 127.0L/14400.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 17.0L/1440.0L*pow(x_0, 4)*x_1 + (991.0L/11520.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 11.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (251.0L/1440.0L)*pow(x_0, 3)*x_1 - 8443.0L/17280.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (11.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (121.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 2483.0L/2880.0L*pow(x_0, 2)*x_1 + (73999.0L/46080.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 11.0L/1440.0L*x_0*pow(x_1, 4) - 121.0L/1440.0L*x_0*pow(x_1, 3) - 1331.0L/2880.0L*x_0*pow(x_1, 2) + (10687.0L/5760.0L)*x_0*x_1 - 646237.0L/230400.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (11.0L/14400.0L)*pow(x_1, 5) + (121.0L/11520.0L)*pow(x_1, 4) + (1331.0L/17280.0L)*pow(x_1, 3) + (14641.0L/46080.0L)*pow(x_1, 2) - 343159.0L/230400.0L*x_1 + 347071.0L/172800.0L;
       return __pp_r18___result;
    }
    static inline O __pp_r19__(const I &x_0, const I &x_1) {
       O __pp_r19___result;
       __pp_r19___result = -1.0L/43200.0L*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 + (11.0L/4800.0L)*pow(x_0, 5) - 1.0L/720.0L*pow(x_0, 4)*pow(x_1, 2) + (7.0L/160.0L)*pow(x_0, 4)*x_1 - 289.0L/11520.0L*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/20.0L)*pow(x_0, 3)*pow(x_1, 2) - 389.0L/1440.0L*pow(x_0, 3)*x_1 + (599.0L/5760.0L)*pow(x_0, 3) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 4) + (17.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) - 199.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (293.0L/320.0L)*pow(x_0, 2)*x_1 - 7921.0L/46080.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (23.0L/480.0L)*x_0*pow(x_1, 4) - 761.0L/1440.0L*x_0*pow(x_1, 3) + (421.0L/320.0L)*x_0*pow(x_1, 2) - 9793.0L/5760.0L*x_0*x_1 + (3041.0L/76800.0L)*x_0 - 19.0L/43200.0L*pow(x_1, 6) + (19.0L/1600.0L)*pow(x_1, 5) - 1159.0L/11520.0L*pow(x_1, 4) + (3857.0L/5760.0L)*pow(x_1, 3) - 67279.0L/46080.0L*pow(x_1, 2) + (34689.0L/25600.0L)*x_1 + 19391.0L/172800.0L;
       return __pp_r19___result;
    }
    static inline O __pp_r20__(const I &x_0, const I &x_1) {
       O __pp_r20___result;
       __pp_r20___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 + (1.0L/960.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/32.0L)*pow(x_0, 4)*x_1 - 127.0L/11520.0L*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) - 227.0L/1440.0L*pow(x_0, 3)*x_1 + (113.0L/5760.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) - 37.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (131.0L/320.0L)*pow(x_0, 2)*x_1 + (5201.0L/46080.0L)*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) - 5.0L/96.0L*x_0*pow(x_1, 4) - 113.0L/1440.0L*x_0*pow(x_1, 3) + (97.0L/320.0L)*x_0*pow(x_1, 2) - 101.0L/180.0L*x_0*x_1 - 1453.0L/3072.0L*x_0 + (109.0L/43200.0L)*pow(x_1, 6) - 9.0L/320.0L*pow(x_1, 5) + (1433.0L/11520.0L)*pow(x_1, 4) - 31.0L/5760.0L*pow(x_1, 3) - 14791.0L/46080.0L*pow(x_1, 2) + (1689.0L/5120.0L)*x_1 + 686569.0L/1382400.0L;
       return __pp_r20___result;
    }
    static inline O __pp_r21__(const I &x_0, const I &x_1) {
       O __pp_r21___result;
       __pp_r21___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 - 1.0L/30.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/12.0L)*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/12.0L*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 - 25.0L/16.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/32.0L)*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/96.0L*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) - 75.0L/64.0L*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 - 3375.0L/512.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (1.0L/960.0L)*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) + (25.0L/128.0L)*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) + (3375.0L/1024.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r21___result;
    }
    static inline O __pp_r22__(const I &x_0, const I &x_1) {
       O __pp_r22___result;
       __pp_r22___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 - 1.0L/30.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/12.0L)*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/12.0L*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 - 25.0L/16.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/32.0L)*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/96.0L*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) - 75.0L/64.0L*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 - 3375.0L/512.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (1.0L/960.0L)*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) + (25.0L/128.0L)*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) + (3375.0L/1024.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r22___result;
    }
    static inline O __pp_r23__(const I &x_0, const I &x_1) {
       O __pp_r23___result;
       __pp_r23___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 - 1.0L/36.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/18.0L)*pow(x_0, 4)*x_1 + (37.0L/144.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/36.0L*pow(x_0, 3)*pow(x_1, 2) - 29.0L/72.0L*pow(x_0, 3)*x_1 - 547.0L/432.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (419.0L/288.0L)*pow(x_0, 2)*x_1 + (8077.0L/2304.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (5.0L/288.0L)*x_0*pow(x_1, 4) + (19.0L/288.0L)*x_0*pow(x_1, 3) - 163.0L/576.0L*x_0*pow(x_1, 2) - 6029.0L/2304.0L*x_0*x_1 - 119107.0L/23040.0L*x_0 - 1.0L/4800.0L*pow(x_1, 6) - 13.0L/2880.0L*pow(x_1, 5) - 83.0L/2304.0L*pow(x_1, 4) - 349.0L/3456.0L*pow(x_1, 3) + (1933.0L/9216.0L)*pow(x_1, 2) + (86339.0L/46080.0L)*x_1 + 1753837.0L/552960.0L;
       return __pp_r23___result;
    }
    static inline O __pp_r24__(const I &x_0, const I &x_1) {
       O __pp_r24___result;
       __pp_r24___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 - 1.0L/36.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/18.0L)*pow(x_0, 4)*x_1 + (37.0L/144.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/36.0L*pow(x_0, 3)*pow(x_1, 2) - 29.0L/72.0L*pow(x_0, 3)*x_1 - 547.0L/432.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (419.0L/288.0L)*pow(x_0, 2)*x_1 + (8077.0L/2304.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (5.0L/288.0L)*x_0*pow(x_1, 4) + (19.0L/288.0L)*x_0*pow(x_1, 3) - 163.0L/576.0L*x_0*pow(x_1, 2) - 6029.0L/2304.0L*x_0*x_1 - 119107.0L/23040.0L*x_0 - 1.0L/4800.0L*pow(x_1, 6) - 13.0L/2880.0L*pow(x_1, 5) - 83.0L/2304.0L*pow(x_1, 4) - 349.0L/3456.0L*pow(x_1, 3) + (1933.0L/9216.0L)*pow(x_1, 2) + (86339.0L/46080.0L)*x_1 + 1753837.0L/552960.0L;
       return __pp_r24___result;
    }
    static inline O __pp_r25__(const I &x_0, const I &x_1) {
       O __pp_r25___result;
       __pp_r25___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 - 1.0L/36.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/18.0L)*pow(x_0, 4)*x_1 + (37.0L/144.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/36.0L*pow(x_0, 3)*pow(x_1, 2) - 29.0L/72.0L*pow(x_0, 3)*x_1 - 547.0L/432.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (419.0L/288.0L)*pow(x_0, 2)*x_1 + (8077.0L/2304.0L)*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/144.0L)*x_0*pow(x_1, 4) + (1.0L/18.0L)*x_0*pow(x_1, 3) - 83.0L/288.0L*x_0*pow(x_1, 2) - 377.0L/144.0L*x_0*x_1 - 11911.0L/2304.0L*x_0 - 1.0L/800.0L*pow(x_1, 6) + (1.0L/144.0L)*pow(x_1, 5) - 1.0L/288.0L*pow(x_1, 4) - 29.0L/432.0L*pow(x_1, 3) + (523.0L/2304.0L)*pow(x_1, 2) + (4327.0L/2304.0L)*x_1 + 10963.0L/3456.0L;
       return __pp_r25___result;
    }
    static inline O __pp_r26__(const I &x_0, const I &x_1) {
       O __pp_r26___result;
       __pp_r26___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 + (3.0L/100.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 4.0L/45.0L*pow(x_0, 4)*x_1 - 17.0L/80.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 + (553.0L/720.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 31.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 2299.0L/1440.0L*pow(x_0, 2)*x_1 - 5579.0L/3840.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (17.0L/480.0L)*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) + (1193.0L/960.0L)*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 + (16339.0L/12800.0L)*x_0 - 11.0L/43200.0L*pow(x_1, 6) - 91.0L/14400.0L*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) - 6139.0L/17280.0L*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) - 310891.0L/230400.0L*x_1 - 294811.0L/921600.0L;
       return __pp_r26___result;
    }
    static inline O __pp_r27__(const I &x_0, const I &x_1) {
       O __pp_r27___result;
       __pp_r27___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 + (3.0L/100.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 4.0L/45.0L*pow(x_0, 4)*x_1 - 17.0L/80.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 + (553.0L/720.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 31.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 2299.0L/1440.0L*pow(x_0, 2)*x_1 - 5579.0L/3840.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (17.0L/480.0L)*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) + (1193.0L/960.0L)*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 + (16339.0L/12800.0L)*x_0 - 11.0L/43200.0L*pow(x_1, 6) - 91.0L/14400.0L*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) - 6139.0L/17280.0L*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) - 310891.0L/230400.0L*x_1 - 294811.0L/921600.0L;
       return __pp_r27___result;
    }
    static inline O __pp_r28__(const I &x_0, const I &x_1) {
       O __pp_r28___result;
       __pp_r28___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 + (3.0L/100.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 4.0L/45.0L*pow(x_0, 4)*x_1 - 17.0L/80.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 + (553.0L/720.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 31.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 2299.0L/1440.0L*pow(x_0, 2)*x_1 - 5579.0L/3840.0L*pow(x_0, 2) - 1.0L/400.0L*x_0*pow(x_1, 5) + (1.0L/40.0L)*x_0*pow(x_1, 4) + (209.0L/720.0L)*x_0*pow(x_1, 3) + (99.0L/80.0L)*x_0*pow(x_1, 2) + (13481.0L/5760.0L)*x_0*x_1 + (12253.0L/9600.0L)*x_0 - 7.0L/5400.0L*pow(x_1, 6) + (37.0L/7200.0L)*pow(x_1, 5) - 21.0L/640.0L*pow(x_1, 4) - 2777.0L/8640.0L*pow(x_1, 3) - 7777.0L/7680.0L*pow(x_1, 2) - 154943.0L/115200.0L*x_1 - 147203.0L/460800.0L;
       return __pp_r28___result;
    }
    static inline O __pp_r29__(const I &x_0, const I &x_1) {
       O __pp_r29___result;
       __pp_r29___result = -1.0L/800.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (13.0L/600.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 17.0L/360.0L*pow(x_0, 4)*x_1 - 3.0L/20.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) + (103.0L/360.0L)*pow(x_0, 3)*x_1 + (373.0L/720.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 31.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 1219.0L/1440.0L*pow(x_0, 2)*x_1 - 3419.0L/3840.0L*pow(x_0, 2) - 19.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/60.0L*x_0*pow(x_1, 4) + (29.0L/720.0L)*x_0*pow(x_1, 3) + (39.0L/80.0L)*x_0*pow(x_1, 2) + (7001.0L/5760.0L)*x_0*x_1 + (5773.0L/9600.0L)*x_0 - 1.0L/1200.0L*pow(x_1, 6) + (97.0L/7200.0L)*pow(x_1, 5) + (19.0L/640.0L)*pow(x_1, 4) - 617.0L/8640.0L*pow(x_1, 3) - 3457.0L/7680.0L*pow(x_1, 2) - 77183.0L/115200.0L*x_1 + 8317.0L/460800.0L;
       return __pp_r29___result;
    }
    static inline O __pp_r30__(const I &x_0, const I &x_1) {
       O __pp_r30___result;
       __pp_r30___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 1.0L/360.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/72.0L)*pow(x_0, 4)*x_1 + (13.0L/720.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/36.0L*pow(x_0, 3)*pow(x_1, 2) - 1.0L/20.0L*pow(x_0, 3)*x_1 - 53.0L/540.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/36.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) + (7.0L/90.0L)*pow(x_0, 2)*x_1 + (137.0L/360.0L)*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) - 7.0L/288.0L*x_0*pow(x_1, 4) - 7.0L/160.0L*x_0*pow(x_1, 3) + (73.0L/2880.0L)*x_0*pow(x_1, 2) - 71.0L/1280.0L*x_0*x_1 - 3671.0L/4608.0L*x_0 - 7.0L/8640.0L*pow(x_1, 6) + (41.0L/2880.0L)*pow(x_1, 5) + (463.0L/11520.0L)*pow(x_1, 4) + (97.0L/17280.0L)*pow(x_1, 3) - 6101.0L/46080.0L*pow(x_1, 2) + (1337.0L/46080.0L)*x_1 + 1821463.0L/2764800.0L;
       return __pp_r30___result;
    }
    static inline O __pp_r31__(const I &x_0, const I &x_1) {
       O __pp_r31___result;
       __pp_r31___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 - 1.0L/30.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/12.0L)*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/12.0L*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 - 25.0L/16.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/32.0L)*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/96.0L*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) - 75.0L/64.0L*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 - 3375.0L/512.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (1.0L/960.0L)*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) + (25.0L/128.0L)*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) + (3375.0L/1024.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r31___result;
    }
    static inline O __pp_r32__(const I &x_0, const I &x_1) {
       O __pp_r32___result;
       __pp_r32___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 - 1.0L/30.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/12.0L)*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/12.0L*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 - 25.0L/16.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/32.0L)*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/96.0L*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) - 75.0L/64.0L*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 - 3375.0L/512.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (1.0L/960.0L)*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) + (25.0L/128.0L)*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) + (3375.0L/1024.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r32___result;
    }
    static inline O __pp_r33__(const I &x_0, const I &x_1) {
       O __pp_r33___result;
       __pp_r33___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 - 1.0L/30.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/12.0L)*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/12.0L*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 - 25.0L/16.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/32.0L)*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/96.0L*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) - 75.0L/64.0L*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 - 3375.0L/512.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (1.0L/960.0L)*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) + (25.0L/128.0L)*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) + (3375.0L/1024.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r33___result;
    }
    static inline O __pp_r34__(const I &x_0, const I &x_1) {
       O __pp_r34___result;
       __pp_r34___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 - 1.0L/30.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/12.0L)*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/12.0L*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 - 25.0L/16.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/32.0L)*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/96.0L*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) - 75.0L/64.0L*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 - 3375.0L/512.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (1.0L/960.0L)*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) + (25.0L/128.0L)*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) + (3375.0L/1024.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r34___result;
    }
    static inline O __pp_r35__(const I &x_0, const I &x_1) {
       O __pp_r35___result;
       __pp_r35___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 - 1.0L/30.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/12.0L)*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/12.0L*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 - 25.0L/16.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/32.0L)*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/96.0L*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) - 75.0L/64.0L*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 - 3375.0L/512.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (1.0L/960.0L)*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) + (25.0L/128.0L)*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) + (3375.0L/1024.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r35___result;
    }
    static inline O __pp_r36__(const I &x_0, const I &x_1) {
       O __pp_r36___result;
       __pp_r36___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 - 1.0L/36.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/18.0L)*pow(x_0, 4)*x_1 + (37.0L/144.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/36.0L*pow(x_0, 3)*pow(x_1, 2) - 29.0L/72.0L*pow(x_0, 3)*x_1 - 547.0L/432.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (419.0L/288.0L)*pow(x_0, 2)*x_1 + (8077.0L/2304.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (5.0L/288.0L)*x_0*pow(x_1, 4) + (19.0L/288.0L)*x_0*pow(x_1, 3) - 163.0L/576.0L*x_0*pow(x_1, 2) - 6029.0L/2304.0L*x_0*x_1 - 119107.0L/23040.0L*x_0 - 1.0L/4800.0L*pow(x_1, 6) - 13.0L/2880.0L*pow(x_1, 5) - 83.0L/2304.0L*pow(x_1, 4) - 349.0L/3456.0L*pow(x_1, 3) + (1933.0L/9216.0L)*pow(x_1, 2) + (86339.0L/46080.0L)*x_1 + 1753837.0L/552960.0L;
       return __pp_r36___result;
    }
    static inline O __pp_r37__(const I &x_0, const I &x_1) {
       O __pp_r37___result;
       __pp_r37___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 + (3.0L/100.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 4.0L/45.0L*pow(x_0, 4)*x_1 - 17.0L/80.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 + (553.0L/720.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 31.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 2299.0L/1440.0L*pow(x_0, 2)*x_1 - 5579.0L/3840.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (17.0L/480.0L)*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) + (1193.0L/960.0L)*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 + (16339.0L/12800.0L)*x_0 - 11.0L/43200.0L*pow(x_1, 6) - 91.0L/14400.0L*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) - 6139.0L/17280.0L*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) - 310891.0L/230400.0L*x_1 - 294811.0L/921600.0L;
       return __pp_r37___result;
    }
    static inline O __pp_r38__(const I &x_0, const I &x_1) {
       O __pp_r38___result;
       __pp_r38___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 - 1.0L/30.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/12.0L)*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/12.0L*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 - 25.0L/16.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/32.0L)*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/96.0L*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) - 75.0L/64.0L*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 - 3375.0L/512.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (1.0L/960.0L)*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) + (25.0L/128.0L)*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) + (3375.0L/1024.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r38___result;
    }
    static inline O __pp_r39__(const I &x_0, const I &x_1) {
       O __pp_r39___result;
       __pp_r39___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 - 1.0L/36.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/18.0L)*pow(x_0, 4)*x_1 + (37.0L/144.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/36.0L*pow(x_0, 3)*pow(x_1, 2) - 29.0L/72.0L*pow(x_0, 3)*x_1 - 547.0L/432.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (419.0L/288.0L)*pow(x_0, 2)*x_1 + (8077.0L/2304.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (5.0L/288.0L)*x_0*pow(x_1, 4) + (19.0L/288.0L)*x_0*pow(x_1, 3) - 163.0L/576.0L*x_0*pow(x_1, 2) - 6029.0L/2304.0L*x_0*x_1 - 119107.0L/23040.0L*x_0 - 1.0L/4800.0L*pow(x_1, 6) - 13.0L/2880.0L*pow(x_1, 5) - 83.0L/2304.0L*pow(x_1, 4) - 349.0L/3456.0L*pow(x_1, 3) + (1933.0L/9216.0L)*pow(x_1, 2) + (86339.0L/46080.0L)*x_1 + 1753837.0L/552960.0L;
       return __pp_r39___result;
    }
    static inline O __pp_r40__(const I &x_0, const I &x_1) {
       O __pp_r40___result;
       __pp_r40___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 + (3.0L/100.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 4.0L/45.0L*pow(x_0, 4)*x_1 - 17.0L/80.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 + (553.0L/720.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 31.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 2299.0L/1440.0L*pow(x_0, 2)*x_1 - 5579.0L/3840.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (17.0L/480.0L)*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) + (1193.0L/960.0L)*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 + (16339.0L/12800.0L)*x_0 - 11.0L/43200.0L*pow(x_1, 6) - 91.0L/14400.0L*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) - 6139.0L/17280.0L*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) - 310891.0L/230400.0L*x_1 - 294811.0L/921600.0L;
       return __pp_r40___result;
    }
    static inline O __pp_r41__(const I &x_0, const I &x_1) {
       O __pp_r41___result;
       __pp_r41___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 + (41.0L/1200.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 4.0L/45.0L*pow(x_0, 4)*x_1 - 127.0L/480.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 + (1481.0L/1440.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 31.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 2299.0L/1440.0L*pow(x_0, 2)*x_1 - 2693.0L/1280.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (17.0L/480.0L)*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) + (1193.0L/960.0L)*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 + (80267.0L/38400.0L)*x_0 - 11.0L/43200.0L*pow(x_1, 6) - 91.0L/14400.0L*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) - 6139.0L/17280.0L*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) - 310891.0L/230400.0L*x_1 - 669811.0L/921600.0L;
       return __pp_r41___result;
    }
    static inline O __pp_r42__(const I &x_0, const I &x_1) {
       O __pp_r42___result;
       __pp_r42___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 + (41.0L/1200.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 4.0L/45.0L*pow(x_0, 4)*x_1 - 127.0L/480.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 + (1481.0L/1440.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 31.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 2299.0L/1440.0L*pow(x_0, 2)*x_1 - 2693.0L/1280.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (17.0L/480.0L)*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) + (1193.0L/960.0L)*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 + (80267.0L/38400.0L)*x_0 - 11.0L/43200.0L*pow(x_1, 6) - 91.0L/14400.0L*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) - 6139.0L/17280.0L*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) - 310891.0L/230400.0L*x_1 - 669811.0L/921600.0L;
       return __pp_r42___result;
    }
    static inline O __pp_r43__(const I &x_0, const I &x_1) {
       O __pp_r43___result;
       __pp_r43___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 + (41.0L/1200.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 4.0L/45.0L*pow(x_0, 4)*x_1 - 127.0L/480.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 + (1481.0L/1440.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 31.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 2299.0L/1440.0L*pow(x_0, 2)*x_1 - 2693.0L/1280.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (17.0L/480.0L)*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) + (1193.0L/960.0L)*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 + (80267.0L/38400.0L)*x_0 - 11.0L/43200.0L*pow(x_1, 6) - 91.0L/14400.0L*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) - 6139.0L/17280.0L*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) - 310891.0L/230400.0L*x_1 - 669811.0L/921600.0L;
       return __pp_r43___result;
    }
    static inline O __pp_r44__(const I &x_0, const I &x_1) {
       O __pp_r44___result;
       __pp_r44___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 + (41.0L/1200.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 4.0L/45.0L*pow(x_0, 4)*x_1 - 127.0L/480.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 + (1481.0L/1440.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 31.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 2299.0L/1440.0L*pow(x_0, 2)*x_1 - 2693.0L/1280.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (17.0L/480.0L)*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) + (1193.0L/960.0L)*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 + (80267.0L/38400.0L)*x_0 - 11.0L/43200.0L*pow(x_1, 6) - 91.0L/14400.0L*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) - 6139.0L/17280.0L*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) - 310891.0L/230400.0L*x_1 - 669811.0L/921600.0L;
       return __pp_r44___result;
    }
    static inline O __pp_r45__(const I &x_0, const I &x_1) {
       O __pp_r45___result;
       __pp_r45___result = -11.0L/10800.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 + (13.0L/600.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 19.0L/720.0L*pow(x_0, 4)*x_1 - 41.0L/240.0L*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (29.0L/180.0L)*pow(x_0, 3)*x_1 + (941.0L/1440.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (7.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 1.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 679.0L/1440.0L*pow(x_0, 2)*x_1 - 1613.0L/1280.0L*pow(x_0, 2) - 1.0L/400.0L*x_0*pow(x_1, 5) - 13.0L/480.0L*x_0*pow(x_1, 4) - 107.0L/1440.0L*x_0*pow(x_1, 3) + (113.0L/960.0L)*x_0*pow(x_1, 2) + (7537.0L/11520.0L)*x_0*x_1 + (41387.0L/38400.0L)*x_0 + (19.0L/43200.0L)*pow(x_1, 6) + (89.0L/14400.0L)*pow(x_1, 5) + (109.0L/3840.0L)*pow(x_1, 4) + (341.0L/17280.0L)*pow(x_1, 3) - 953.0L/5120.0L*pow(x_1, 2) - 77611.0L/230400.0L*x_1 - 203251.0L/921600.0L;
       return __pp_r45___result;
    }
    static inline O __pp_r46__(const I &x_0, const I &x_1) {
       O __pp_r46___result;
       __pp_r46___result = -11.0L/10800.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 + (13.0L/600.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 19.0L/720.0L*pow(x_0, 4)*x_1 - 41.0L/240.0L*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (29.0L/180.0L)*pow(x_0, 3)*x_1 + (941.0L/1440.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (7.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 1.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 679.0L/1440.0L*pow(x_0, 2)*x_1 - 1613.0L/1280.0L*pow(x_0, 2) - 1.0L/400.0L*x_0*pow(x_1, 5) - 13.0L/480.0L*x_0*pow(x_1, 4) - 107.0L/1440.0L*x_0*pow(x_1, 3) + (113.0L/960.0L)*x_0*pow(x_1, 2) + (7537.0L/11520.0L)*x_0*x_1 + (41387.0L/38400.0L)*x_0 + (19.0L/43200.0L)*pow(x_1, 6) + (89.0L/14400.0L)*pow(x_1, 5) + (109.0L/3840.0L)*pow(x_1, 4) + (341.0L/17280.0L)*pow(x_1, 3) - 953.0L/5120.0L*pow(x_1, 2) - 77611.0L/230400.0L*x_1 - 203251.0L/921600.0L;
       return __pp_r46___result;
    }
    static inline O __pp_r47__(const I &x_0, const I &x_1) {
       O __pp_r47___result;
       __pp_r47___result = (7.0L/3600.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 - 11.0L/600.0L*pow(x_0, 5) + (1.0L/80.0L)*pow(x_0, 4)*pow(x_1, 2) + (53.0L/720.0L)*pow(x_0, 4)*x_1 + (13.0L/240.0L)*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) - 13.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) - 13.0L/45.0L*pow(x_0, 3)*x_1 - 31.0L/1440.0L*pow(x_0, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 4) + (4.0L/45.0L)*pow(x_0, 2)*pow(x_1, 3) + (53.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (779.0L/1440.0L)*pow(x_0, 2)*x_1 - 31.0L/256.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 19.0L/480.0L*x_0*pow(x_1, 4) - 269.0L/1440.0L*x_0*pow(x_1, 3) - 373.0L/960.0L*x_0*pow(x_1, 2) - 1117.0L/2304.0L*x_0*x_1 + (2021.0L/38400.0L)*x_0 + (7.0L/14400.0L)*pow(x_1, 6) + (107.0L/14400.0L)*pow(x_1, 5) + (163.0L/3840.0L)*pow(x_1, 4) + (1799.0L/17280.0L)*pow(x_1, 3) + (101.0L/1024.0L)*pow(x_1, 2) + (40487.0L/230400.0L)*x_1 + 151043.0L/921600.0L;
       return __pp_r47___result;
    }
    static inline O __pp_r48__(const I &x_0, const I &x_1) {
       O __pp_r48___result;
       __pp_r48___result = (7.0L/3600.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 - 11.0L/600.0L*pow(x_0, 5) + (1.0L/80.0L)*pow(x_0, 4)*pow(x_1, 2) + (53.0L/720.0L)*pow(x_0, 4)*x_1 + (13.0L/240.0L)*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) - 13.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) - 13.0L/45.0L*pow(x_0, 3)*x_1 - 31.0L/1440.0L*pow(x_0, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 4) + (4.0L/45.0L)*pow(x_0, 2)*pow(x_1, 3) + (53.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (779.0L/1440.0L)*pow(x_0, 2)*x_1 - 31.0L/256.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 19.0L/480.0L*x_0*pow(x_1, 4) - 269.0L/1440.0L*x_0*pow(x_1, 3) - 373.0L/960.0L*x_0*pow(x_1, 2) - 1117.0L/2304.0L*x_0*x_1 + (2021.0L/38400.0L)*x_0 + (7.0L/14400.0L)*pow(x_1, 6) + (107.0L/14400.0L)*pow(x_1, 5) + (163.0L/3840.0L)*pow(x_1, 4) + (1799.0L/17280.0L)*pow(x_1, 3) + (101.0L/1024.0L)*pow(x_1, 2) + (40487.0L/230400.0L)*x_1 + 151043.0L/921600.0L;
       return __pp_r48___result;
    }
    static inline O __pp_r49__(const I &x_0, const I &x_1) {
       O __pp_r49___result;
       __pp_r49___result = -17.0L/21600.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (3.0L/160.0L)*pow(x_0, 5) + (7.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/24.0L*pow(x_0, 4)*x_1 - 7.0L/45.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (41.0L/180.0L)*pow(x_0, 3)*x_1 + (1759.0L/2880.0L)*pow(x_0, 3) + (13.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (53.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 37.0L/60.0L*pow(x_0, 2)*x_1 - 13741.0L/11520.0L*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) - 5.0L/96.0L*x_0*pow(x_1, 4) + (37.0L/1440.0L)*x_0*pow(x_1, 3) - 59.0L/960.0L*x_0*pow(x_1, 2) + (9361.0L/11520.0L)*x_0*x_1 + (49.0L/48.0L)*x_0 + (29.0L/43200.0L)*pow(x_1, 6) + (1.0L/320.0L)*pow(x_1, 5) + (863.0L/11520.0L)*pow(x_1, 4) - 301.0L/5760.0L*pow(x_1, 3) - 3601.0L/46080.0L*pow(x_1, 2) - 6221.0L/15360.0L*x_1 - 553537.0L/2764800.0L;
       return __pp_r49___result;
    }
    static inline O __pp_r50__(const I &x_0, const I &x_1) {
       O __pp_r50___result;
       __pp_r50___result = -7.0L/8640.0L*pow(x_0, 6) + (1.0L/360.0L)*pow(x_0, 5)*x_1 + (11.0L/576.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/288.0L*pow(x_0, 4)*x_1 - 1817.0L/11520.0L*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/36.0L*pow(x_0, 3)*pow(x_1, 2) + (101.0L/480.0L)*pow(x_0, 3)*x_1 + (10679.0L/17280.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/36.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) - 1651.0L/2880.0L*pow(x_0, 2)*x_1 - 55589.0L/46080.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) - 7.0L/288.0L*x_0*pow(x_1, 4) - 7.0L/160.0L*x_0*pow(x_1, 3) + (73.0L/2880.0L)*x_0*pow(x_1, 2) + (91.0L/120.0L)*x_0*x_1 + (9533.0L/9216.0L)*x_0 - 7.0L/8640.0L*pow(x_1, 6) + (41.0L/2880.0L)*pow(x_1, 5) + (463.0L/11520.0L)*pow(x_1, 4) + (97.0L/17280.0L)*pow(x_1, 3) - 6101.0L/46080.0L*pow(x_1, 2) - 17413.0L/46080.0L*x_1 - 284581.0L/1382400.0L;
       return __pp_r50___result;
    }
    static inline O __pp_r51__(const I &x_0, const I &x_1) {
       O __pp_r51___result;
       __pp_r51___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 + (7.0L/320.0L)*pow(x_0, 5) - 7.0L/288.0L*pow(x_0, 4)*x_1 - 659.0L/3840.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) + (223.0L/1440.0L)*pow(x_0, 3)*x_1 + (3773.0L/5760.0L)*pow(x_0, 3) + (1.0L/18.0L)*pow(x_0, 2)*pow(x_1, 3) - 1.0L/40.0L*pow(x_0, 2)*pow(x_1, 2) - 1331.0L/2880.0L*pow(x_0, 2)*x_1 - 6461.0L/5120.0L*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) - 1.0L/96.0L*x_0*pow(x_1, 4) - 143.0L/1440.0L*x_0*pow(x_1, 3) + (131.0L/960.0L)*x_0*pow(x_1, 2) + (233.0L/360.0L)*x_0*x_1 + (16571.0L/15360.0L)*x_0 - 1.0L/960.0L*pow(x_1, 6) + (49.0L/2880.0L)*pow(x_1, 5) + (101.0L/3840.0L)*pow(x_1, 4) + (737.0L/17280.0L)*pow(x_1, 3) - 2887.0L/15360.0L*pow(x_1, 2) - 3073.0L/9216.0L*x_1 - 101687.0L/460800.0L;
       return __pp_r51___result;
    }
    static inline O __pp_r52__(const I &x_0, const I &x_1) {
       O __pp_r52___result;
       __pp_r52___result = -11.0L/10800.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 + (13.0L/600.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 19.0L/720.0L*pow(x_0, 4)*x_1 - 41.0L/240.0L*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (29.0L/180.0L)*pow(x_0, 3)*x_1 + (941.0L/1440.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (7.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 1.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 679.0L/1440.0L*pow(x_0, 2)*x_1 - 1613.0L/1280.0L*pow(x_0, 2) - 1.0L/400.0L*x_0*pow(x_1, 5) - 13.0L/480.0L*x_0*pow(x_1, 4) - 107.0L/1440.0L*x_0*pow(x_1, 3) + (113.0L/960.0L)*x_0*pow(x_1, 2) + (7537.0L/11520.0L)*x_0*x_1 + (41387.0L/38400.0L)*x_0 + (19.0L/43200.0L)*pow(x_1, 6) + (149.0L/14400.0L)*pow(x_1, 5) + (149.0L/3840.0L)*pow(x_1, 4) + (521.0L/17280.0L)*pow(x_1, 3) - 2779.0L/15360.0L*pow(x_1, 2) - 77311.0L/230400.0L*x_1 - 203131.0L/921600.0L;
       return __pp_r52___result;
    }
    static inline O __pp_r53__(const I &x_0, const I &x_1) {
       O __pp_r53___result;
       __pp_r53___result = (47.0L/21600.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 - 17.0L/800.0L*pow(x_0, 5) + (23.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (7.0L/120.0L)*pow(x_0, 4)*x_1 + (5.0L/72.0L)*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 17.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) - 2.0L/9.0L*pow(x_0, 3)*x_1 - 37.0L/576.0L*pow(x_0, 3) + (17.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/20.0L)*pow(x_0, 2)*pow(x_1, 3) + (43.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (19.0L/48.0L)*pow(x_0, 2)*x_1 - 619.0L/11520.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) - 31.0L/480.0L*x_0*pow(x_1, 4) - 25.0L/288.0L*x_0*pow(x_1, 3) - 109.0L/192.0L*x_0*pow(x_1, 2) - 3761.0L/11520.0L*x_0*x_1 - 83.0L/19200.0L*x_0 + (31.0L/43200.0L)*pow(x_1, 6) + (7.0L/1600.0L)*pow(x_1, 5) + (205.0L/2304.0L)*pow(x_1, 4) + (37.0L/1152.0L)*pow(x_1, 3) + (9521.0L/46080.0L)*pow(x_1, 2) + (8261.0L/76800.0L)*x_1 + 101869.0L/552960.0L;
       return __pp_r53___result;
    }
    static inline O __pp_r54__(const I &x_0, const I &x_1) {
       O __pp_r54___result;
       __pp_r54___result = (31.0L/14400.0L)*pow(x_0, 6) - 11.0L/1800.0L*pow(x_0, 5)*x_1 - 301.0L/14400.0L*pow(x_0, 5) + (7.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (89.0L/1440.0L)*pow(x_0, 4)*x_1 + (155.0L/2304.0L)*pow(x_0, 4) - 13.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 23.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) - 23.0L/96.0L*pow(x_0, 3)*x_1 - 197.0L/3456.0L*pow(x_0, 3) + (1.0L/160.0L)*pow(x_0, 2)*pow(x_1, 4) + (7.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (19.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (253.0L/576.0L)*pow(x_0, 2)*x_1 - 3101.0L/46080.0L*pow(x_0, 2) - 11.0L/1800.0L*x_0*pow(x_1, 5) - 53.0L/1440.0L*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) - 277.0L/576.0L*x_0*pow(x_1, 2) - 731.0L/1920.0L*x_0*x_1 + (2129.0L/230400.0L)*x_0 - 11.0L/14400.0L*pow(x_1, 6) + (223.0L/14400.0L)*pow(x_1, 5) + (125.0L/2304.0L)*pow(x_1, 4) + (311.0L/3456.0L)*pow(x_1, 3) + (7021.0L/46080.0L)*pow(x_1, 2) + (31033.0L/230400.0L)*x_1 + 12343.0L/69120.0L;
       return __pp_r54___result;
    }
    static inline O __pp_r55__(const I &x_0, const I &x_1) {
       O __pp_r55___result;
       __pp_r55___result = (83.0L/43200.0L)*pow(x_0, 6) - 3.0L/400.0L*pow(x_0, 5)*x_1 - 29.0L/1600.0L*pow(x_0, 5) + (1.0L/90.0L)*pow(x_0, 4)*pow(x_1, 2) + (109.0L/1440.0L)*pow(x_0, 4)*x_1 + (41.0L/768.0L)*pow(x_0, 4) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/10.0L*pow(x_0, 3)*pow(x_1, 2) - 85.0L/288.0L*pow(x_0, 3)*x_1 - 23.0L/1152.0L*pow(x_0, 3) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 4) + (19.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (5.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (317.0L/576.0L)*pow(x_0, 2)*x_1 - 629.0L/5120.0L*pow(x_0, 2) - 3.0L/400.0L*x_0*pow(x_1, 5) - 11.0L/480.0L*x_0*pow(x_1, 4) - 61.0L/288.0L*x_0*pow(x_1, 3) - 71.0L/192.0L*x_0*pow(x_1, 2) - 2833.0L/5760.0L*x_0*x_1 + (4123.0L/76800.0L)*x_0 - 43.0L/43200.0L*pow(x_1, 6) + (263.0L/14400.0L)*pow(x_1, 5) + (31.0L/768.0L)*pow(x_1, 4) + (439.0L/3456.0L)*pow(x_1, 3) + (1487.0L/15360.0L)*pow(x_1, 2) + (41273.0L/230400.0L)*x_1 + 3773.0L/23040.0L;
       return __pp_r55___result;
    }
    static inline O __pp_r56__(const I &x_0, const I &x_1) {
       O __pp_r56___result;
       __pp_r56___result = (7.0L/3600.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 - 11.0L/600.0L*pow(x_0, 5) + (1.0L/80.0L)*pow(x_0, 4)*pow(x_1, 2) + (53.0L/720.0L)*pow(x_0, 4)*x_1 + (13.0L/240.0L)*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) - 13.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) - 13.0L/45.0L*pow(x_0, 3)*x_1 - 31.0L/1440.0L*pow(x_0, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 4) + (4.0L/45.0L)*pow(x_0, 2)*pow(x_1, 3) + (53.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (779.0L/1440.0L)*pow(x_0, 2)*x_1 - 31.0L/256.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 19.0L/480.0L*x_0*pow(x_1, 4) - 269.0L/1440.0L*x_0*pow(x_1, 3) - 373.0L/960.0L*x_0*pow(x_1, 2) - 1117.0L/2304.0L*x_0*x_1 + (2021.0L/38400.0L)*x_0 + (7.0L/14400.0L)*pow(x_1, 6) + (167.0L/14400.0L)*pow(x_1, 5) + (203.0L/3840.0L)*pow(x_1, 4) + (1979.0L/17280.0L)*pow(x_1, 3) + (319.0L/3072.0L)*pow(x_1, 2) + (40787.0L/230400.0L)*x_1 + 151163.0L/921600.0L;
       return __pp_r56___result;
    }
    static inline O __pp_r57__(const I &x_0, const I &x_1) {
       O __pp_r57___result;
       __pp_r57___result = (19.0L/43200.0L)*pow(x_0, 6) - 1.0L/3600.0L*pow(x_0, 5)*x_1 - 7.0L/960.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/32.0L)*pow(x_0, 4)*x_1 + (65.0L/768.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/12.0L*pow(x_0, 3)*pow(x_1, 2) - 35.0L/96.0L*pow(x_0, 3)*x_1 - 75.0L/128.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (325.0L/192.0L)*pow(x_0, 2)*x_1 + (6625.0L/3072.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/96.0L*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) - 75.0L/64.0L*x_0*pow(x_1, 2) - 1375.0L/384.0L*x_0*x_1 - 12125.0L/3072.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (1.0L/960.0L)*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) + (25.0L/128.0L)*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) + (8875.0L/3072.0L)*x_1 + 4375.0L/1536.0L;
       return __pp_r57___result;
    }
    static inline O __pp_r58__(const I &x_0, const I &x_1) {
       O __pp_r58___result;
       __pp_r58___result = (19.0L/43200.0L)*pow(x_0, 6) - 1.0L/3600.0L*pow(x_0, 5)*x_1 - 7.0L/960.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/32.0L)*pow(x_0, 4)*x_1 + (65.0L/768.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/12.0L*pow(x_0, 3)*pow(x_1, 2) - 35.0L/96.0L*pow(x_0, 3)*x_1 - 75.0L/128.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (325.0L/192.0L)*pow(x_0, 2)*x_1 + (6625.0L/3072.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/96.0L*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) - 75.0L/64.0L*x_0*pow(x_1, 2) - 1375.0L/384.0L*x_0*x_1 - 12125.0L/3072.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (1.0L/960.0L)*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) + (25.0L/128.0L)*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) + (8875.0L/3072.0L)*x_1 + 4375.0L/1536.0L;
       return __pp_r58___result;
    }
    static inline O __pp_r59__(const I &x_0, const I &x_1) {
       O __pp_r59___result;
       __pp_r59___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 + (311.0L/14400.0L)*pow(x_0, 5) - 59.0L/1440.0L*pow(x_0, 4)*x_1 - 1729.0L/11520.0L*pow(x_0, 4) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (151.0L/1440.0L)*pow(x_0, 3)*x_1 + (7451.0L/17280.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) + (481.0L/2880.0L)*pow(x_0, 2)*x_1 - 14869.0L/46080.0L*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 4) - 7.0L/180.0L*x_0*pow(x_1, 3) - 589.0L/1440.0L*x_0*pow(x_1, 2) - 12689.0L/11520.0L*x_0*x_1 - 166789.0L/230400.0L*x_0 + (1.0L/7200.0L)*pow(x_1, 5) + (7.0L/1440.0L)*pow(x_1, 4) + (589.0L/8640.0L)*pow(x_1, 3) + (1379.0L/2880.0L)*pow(x_1, 2) + (73583.0L/57600.0L)*x_1 + 3048191.0L/2764800.0L;
       return __pp_r59___result;
    }
    static inline O __pp_r60__(const I &x_0, const I &x_1) {
       O __pp_r60___result;
       __pp_r60___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 + (311.0L/14400.0L)*pow(x_0, 5) - 59.0L/1440.0L*pow(x_0, 4)*x_1 - 1729.0L/11520.0L*pow(x_0, 4) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (151.0L/1440.0L)*pow(x_0, 3)*x_1 + (7451.0L/17280.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) + (481.0L/2880.0L)*pow(x_0, 2)*x_1 - 14869.0L/46080.0L*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 4) - 7.0L/180.0L*x_0*pow(x_1, 3) - 589.0L/1440.0L*x_0*pow(x_1, 2) - 12689.0L/11520.0L*x_0*x_1 - 166789.0L/230400.0L*x_0 + (1.0L/7200.0L)*pow(x_1, 5) + (7.0L/1440.0L)*pow(x_1, 4) + (589.0L/8640.0L)*pow(x_1, 3) + (1379.0L/2880.0L)*pow(x_1, 2) + (73583.0L/57600.0L)*x_1 + 3048191.0L/2764800.0L;
       return __pp_r60___result;
    }
    static inline O __pp_r61__(const I &x_0, const I &x_1) {
       O __pp_r61___result;
       __pp_r61___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 + (157.0L/4800.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 139.0L/1440.0L*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/10.0L)*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 + (5897.0L/5760.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 19.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 4639.0L/2880.0L*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) + (13.0L/240.0L)*x_0*pow(x_1, 4) + (73.0L/180.0L)*x_0*pow(x_1, 3) + (219.0L/160.0L)*x_0*pow(x_1, 2) + (28271.0L/11520.0L)*x_0*x_1 + (162857.0L/76800.0L)*x_0 - 1.0L/2160.0L*pow(x_1, 6) - 79.0L/7200.0L*pow(x_1, 5) - 17.0L/160.0L*pow(x_1, 4) - 4531.0L/8640.0L*pow(x_1, 3) - 1247.0L/960.0L*pow(x_1, 2) - 90257.0L/57600.0L*x_1 - 731563.0L/921600.0L;
       return __pp_r61___result;
    }
    static inline O __pp_r62__(const I &x_0, const I &x_1) {
       O __pp_r62___result;
       __pp_r62___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 + (157.0L/4800.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 139.0L/1440.0L*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/10.0L)*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 + (5897.0L/5760.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 19.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 4639.0L/2880.0L*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) + (13.0L/240.0L)*x_0*pow(x_1, 4) + (73.0L/180.0L)*x_0*pow(x_1, 3) + (219.0L/160.0L)*x_0*pow(x_1, 2) + (28271.0L/11520.0L)*x_0*x_1 + (162857.0L/76800.0L)*x_0 - 1.0L/2160.0L*pow(x_1, 6) - 79.0L/7200.0L*pow(x_1, 5) - 17.0L/160.0L*pow(x_1, 4) - 4531.0L/8640.0L*pow(x_1, 3) - 1247.0L/960.0L*pow(x_1, 2) - 90257.0L/57600.0L*x_1 - 731563.0L/921600.0L;
       return __pp_r62___result;
    }
    static inline O __pp_r63__(const I &x_0, const I &x_1) {
       O __pp_r63___result;
       __pp_r63___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 + (157.0L/4800.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 139.0L/1440.0L*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/10.0L)*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 + (5897.0L/5760.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 19.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 4639.0L/2880.0L*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) + (13.0L/240.0L)*x_0*pow(x_1, 4) + (73.0L/180.0L)*x_0*pow(x_1, 3) + (219.0L/160.0L)*x_0*pow(x_1, 2) + (28271.0L/11520.0L)*x_0*x_1 + (162857.0L/76800.0L)*x_0 - 1.0L/2160.0L*pow(x_1, 6) - 79.0L/7200.0L*pow(x_1, 5) - 17.0L/160.0L*pow(x_1, 4) - 4531.0L/8640.0L*pow(x_1, 3) - 1247.0L/960.0L*pow(x_1, 2) - 90257.0L/57600.0L*x_1 - 731563.0L/921600.0L;
       return __pp_r63___result;
    }
    static inline O __pp_r64__(const I &x_0, const I &x_1) {
       O __pp_r64___result;
       __pp_r64___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 + (157.0L/4800.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 139.0L/1440.0L*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/10.0L)*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 + (5897.0L/5760.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 19.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 4639.0L/2880.0L*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) + (13.0L/240.0L)*x_0*pow(x_1, 4) + (73.0L/180.0L)*x_0*pow(x_1, 3) + (219.0L/160.0L)*x_0*pow(x_1, 2) + (28271.0L/11520.0L)*x_0*x_1 + (162857.0L/76800.0L)*x_0 - 1.0L/2160.0L*pow(x_1, 6) - 79.0L/7200.0L*pow(x_1, 5) - 17.0L/160.0L*pow(x_1, 4) - 4531.0L/8640.0L*pow(x_1, 3) - 1247.0L/960.0L*pow(x_1, 2) - 90257.0L/57600.0L*x_1 - 731563.0L/921600.0L;
       return __pp_r64___result;
    }
    static inline O __pp_r65__(const I &x_0, const I &x_1) {
       O __pp_r65___result;
       __pp_r65___result = -1.0L/675.0L*pow(x_0, 6) + (13.0L/1800.0L)*pow(x_0, 5)*x_1 + (59.0L/1800.0L)*pow(x_0, 5) - 1.0L/180.0L*pow(x_0, 4)*pow(x_1, 2) - 23.0L/240.0L*pow(x_0, 4)*x_1 - 47.0L/180.0L*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) + (37.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (11.0L/20.0L)*pow(x_0, 3)*x_1 + (4423.0L/4320.0L)*pow(x_0, 3) - 1.0L/720.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/10.0L*pow(x_0, 2)*pow(x_1, 3) - 263.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 773.0L/480.0L*pow(x_0, 2)*x_1 - 24197.0L/11520.0L*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) + (41.0L/1440.0L)*x_0*pow(x_1, 4) + (151.0L/480.0L)*x_0*pow(x_1, 3) + (3539.0L/2880.0L)*x_0*pow(x_1, 2) + (9019.0L/3840.0L)*x_0*x_1 + (240641.0L/115200.0L)*x_0 - 1.0L/43200.0L*pow(x_1, 6) - 37.0L/4800.0L*pow(x_1, 5) - 713.0L/11520.0L*pow(x_1, 4) - 691.0L/1920.0L*pow(x_1, 3) - 47297.0L/46080.0L*pow(x_1, 2) - 34579.0L/25600.0L*x_1 - 2008793.0L/2764800.0L;
       return __pp_r65___result;
    }
    static inline O __pp_r66__(const I &x_0, const I &x_1) {
       O __pp_r66___result;
       __pp_r66___result = -1.0L/675.0L*pow(x_0, 6) + (13.0L/1800.0L)*pow(x_0, 5)*x_1 + (59.0L/1800.0L)*pow(x_0, 5) - 1.0L/180.0L*pow(x_0, 4)*pow(x_1, 2) - 23.0L/240.0L*pow(x_0, 4)*x_1 - 47.0L/180.0L*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) + (37.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (11.0L/20.0L)*pow(x_0, 3)*x_1 + (4423.0L/4320.0L)*pow(x_0, 3) - 1.0L/720.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/10.0L*pow(x_0, 2)*pow(x_1, 3) - 263.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 773.0L/480.0L*pow(x_0, 2)*x_1 - 24197.0L/11520.0L*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) + (41.0L/1440.0L)*x_0*pow(x_1, 4) + (151.0L/480.0L)*x_0*pow(x_1, 3) + (3539.0L/2880.0L)*x_0*pow(x_1, 2) + (9019.0L/3840.0L)*x_0*x_1 + (240641.0L/115200.0L)*x_0 - 1.0L/43200.0L*pow(x_1, 6) - 37.0L/4800.0L*pow(x_1, 5) - 713.0L/11520.0L*pow(x_1, 4) - 691.0L/1920.0L*pow(x_1, 3) - 47297.0L/46080.0L*pow(x_1, 2) - 34579.0L/25600.0L*x_1 - 2008793.0L/2764800.0L;
       return __pp_r66___result;
    }
    static inline O __pp_r67__(const I &x_0, const I &x_1) {
       O __pp_r67___result;
       __pp_r67___result = -17.0L/21600.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (73.0L/3600.0L)*pow(x_0, 5) + (7.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/30.0L*pow(x_0, 4)*x_1 - 241.0L/1440.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/45.0L*pow(x_0, 3)*pow(x_1, 2) + (7.0L/40.0L)*pow(x_0, 3)*x_1 + (2803.0L/4320.0L)*pow(x_0, 3) + (13.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/40.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 233.0L/480.0L*pow(x_0, 2)*x_1 - 14477.0L/11520.0L*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) - 49.0L/1440.0L*x_0*pow(x_1, 4) - 29.0L/480.0L*x_0*pow(x_1, 3) + (299.0L/2880.0L)*x_0*pow(x_1, 2) + (2539.0L/3840.0L)*x_0*x_1 + (124001.0L/115200.0L)*x_0 + (29.0L/43200.0L)*pow(x_1, 6) + (23.0L/4800.0L)*pow(x_1, 5) + (367.0L/11520.0L)*pow(x_1, 4) + (29.0L/1920.0L)*pow(x_1, 3) - 8417.0L/46080.0L*pow(x_1, 2) - 8659.0L/25600.0L*x_1 - 609113.0L/2764800.0L;
       return __pp_r67___result;
    }
    static inline O __pp_r68__(const I &x_0, const I &x_1) {
       O __pp_r68___result;
       __pp_r68___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 + (157.0L/4800.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 139.0L/1440.0L*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/10.0L)*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 + (5897.0L/5760.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 19.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 4639.0L/2880.0L*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) + (11.0L/480.0L)*x_0*pow(x_1, 4) + (449.0L/1440.0L)*x_0*pow(x_1, 3) + (393.0L/320.0L)*x_0*pow(x_1, 2) + (1691.0L/720.0L)*x_0*x_1 + (160427.0L/76800.0L)*x_0 - 13.0L/8640.0L*pow(x_1, 6) - 143.0L/14400.0L*pow(x_1, 5) - 81.0L/1280.0L*pow(x_1, 4) - 6227.0L/17280.0L*pow(x_1, 3) - 15767.0L/15360.0L*pow(x_1, 2) - 311213.0L/230400.0L*x_1 - 334799.0L/460800.0L;
       return __pp_r68___result;
    }
    static inline O __pp_r69__(const I &x_0, const I &x_1) {
       O __pp_r69___result;
       __pp_r69___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 + (157.0L/4800.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 139.0L/1440.0L*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/10.0L)*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 + (5897.0L/5760.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 19.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 4639.0L/2880.0L*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) + (11.0L/480.0L)*x_0*pow(x_1, 4) + (449.0L/1440.0L)*x_0*pow(x_1, 3) + (393.0L/320.0L)*x_0*pow(x_1, 2) + (1691.0L/720.0L)*x_0*x_1 + (160427.0L/76800.0L)*x_0 - 13.0L/8640.0L*pow(x_1, 6) - 143.0L/14400.0L*pow(x_1, 5) - 81.0L/1280.0L*pow(x_1, 4) - 6227.0L/17280.0L*pow(x_1, 3) - 15767.0L/15360.0L*pow(x_1, 2) - 311213.0L/230400.0L*x_1 - 334799.0L/460800.0L;
       return __pp_r69___result;
    }
    static inline O __pp_r70__(const I &x_0, const I &x_1) {
       O __pp_r70___result;
       __pp_r70___result = -7.0L/8640.0L*pow(x_0, 6) + (1.0L/360.0L)*pow(x_0, 5)*x_1 + (97.0L/4800.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 49.0L/1440.0L*pow(x_0, 4)*x_1 - 643.0L/3840.0L*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/40.0L*pow(x_0, 3)*pow(x_1, 2) + (251.0L/1440.0L)*pow(x_0, 3)*x_1 + (3737.0L/5760.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (7.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) - 1399.0L/2880.0L*pow(x_0, 2)*x_1 - 19303.0L/15360.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) - 19.0L/480.0L*x_0*pow(x_1, 4) - 91.0L/1440.0L*x_0*pow(x_1, 3) + (33.0L/320.0L)*x_0*pow(x_1, 2) + (119.0L/180.0L)*x_0*x_1 + (82667.0L/76800.0L)*x_0 - 7.0L/8640.0L*pow(x_1, 6) + (37.0L/14400.0L)*pow(x_1, 5) + (39.0L/1280.0L)*pow(x_1, 4) + (253.0L/17280.0L)*pow(x_1, 3) - 2807.0L/15360.0L*pow(x_1, 2) - 77933.0L/230400.0L*x_1 - 101519.0L/460800.0L;
       return __pp_r70___result;
    }
    static inline O __pp_r71__(const I &x_0, const I &x_1) {
       O __pp_r71___result;
       __pp_r71___result = (47.0L/21600.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 - 71.0L/3600.0L*pow(x_0, 5) + (23.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/15.0L)*pow(x_0, 4)*x_1 + (83.0L/1440.0L)*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 11.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) - 11.0L/40.0L*pow(x_0, 3)*x_1 - 113.0L/4320.0L*pow(x_0, 3) + (17.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (3.0L/40.0L)*pow(x_0, 2)*pow(x_1, 3) + (169.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (253.0L/480.0L)*pow(x_0, 2)*x_1 - 271.0L/2304.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) - 67.0L/1440.0L*x_0*pow(x_1, 4) - 83.0L/480.0L*x_0*pow(x_1, 3) - 1159.0L/2880.0L*x_0*pow(x_1, 2) - 367.0L/768.0L*x_0*x_1 + (5903.0L/115200.0L)*x_0 + (31.0L/43200.0L)*pow(x_1, 6) + (29.0L/4800.0L)*pow(x_1, 5) + (529.0L/11520.0L)*pow(x_1, 4) + (191.0L/1920.0L)*pow(x_1, 3) + (941.0L/9216.0L)*pow(x_1, 2) + (4463.0L/25600.0L)*x_1 + 453769.0L/2764800.0L;
       return __pp_r71___result;
    }
    static inline O __pp_r72__(const I &x_0, const I &x_1) {
       O __pp_r72___result;
       __pp_r72___result = (139.0L/43200.0L)*pow(x_0, 6) - 1.0L/100.0L*pow(x_0, 5)*x_1 - 629.0L/14400.0L*pow(x_0, 5) + (23.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (47.0L/480.0L)*pow(x_0, 4)*x_1 + (2329.0L/11520.0L)*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 11.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) - 59.0L/160.0L*pow(x_0, 3)*x_1 - 7337.0L/17280.0L*pow(x_0, 3) + (17.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (3.0L/40.0L)*pow(x_0, 2)*pow(x_1, 3) + (169.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (641.0L/960.0L)*pow(x_0, 2)*x_1 + (4181.0L/9216.0L)*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) - 67.0L/1440.0L*x_0*pow(x_1, 4) - 83.0L/480.0L*x_0*pow(x_1, 3) - 1159.0L/2880.0L*x_0*pow(x_1, 2) - 7.0L/12.0L*x_0*x_1 - 84179.0L/230400.0L*x_0 + (31.0L/43200.0L)*pow(x_1, 6) + (29.0L/4800.0L)*pow(x_1, 5) + (529.0L/11520.0L)*pow(x_1, 4) + (191.0L/1920.0L)*pow(x_1, 3) + (941.0L/9216.0L)*pow(x_1, 2) + (5273.0L/25600.0L)*x_1 + 396377.0L/1382400.0L;
       return __pp_r72___result;
    }
    static inline O __pp_r73__(const I &x_0, const I &x_1) {
       O __pp_r73___result;
       __pp_r73___result = (11.0L/43200.0L)*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 181.0L/14400.0L*pow(x_0, 5) + (7.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (29.0L/1440.0L)*pow(x_0, 4)*x_1 + (761.0L/11520.0L)*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 2.0L/45.0L*pow(x_0, 3)*pow(x_1, 2) - 139.0L/1440.0L*pow(x_0, 3)*x_1 - 1849.0L/17280.0L*pow(x_0, 3) + (13.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (13.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (71.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (551.0L/2880.0L)*pow(x_0, 2)*x_1 + (1697.0L/46080.0L)*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) - 53.0L/1440.0L*x_0*pow(x_1, 4) - 151.0L/1440.0L*x_0*pow(x_1, 3) - 473.0L/2880.0L*x_0*pow(x_1, 2) - 959.0L/5760.0L*x_0*x_1 - 16951.0L/230400.0L*x_0 + (29.0L/43200.0L)*pow(x_1, 6) + (73.0L/14400.0L)*pow(x_1, 5) + (431.0L/11520.0L)*pow(x_1, 4) + (1033.0L/17280.0L)*pow(x_1, 3) - 97.0L/46080.0L*pow(x_1, 2) + (13843.0L/230400.0L)*x_1 + 34841.0L/172800.0L;
       return __pp_r73___result;
    }
    static inline O __pp_r74__(const I &x_0, const I &x_1) {
       O __pp_r74___result;
       __pp_r74___result = (31.0L/14400.0L)*pow(x_0, 6) - 11.0L/1800.0L*pow(x_0, 5)*x_1 - 19.0L/960.0L*pow(x_0, 5) + (7.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (19.0L/288.0L)*pow(x_0, 4)*x_1 + (221.0L/3840.0L)*pow(x_0, 4) - 13.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/8.0L*pow(x_0, 3)*pow(x_1, 2) - 397.0L/1440.0L*pow(x_0, 3)*x_1 - 151.0L/5760.0L*pow(x_0, 3) + (1.0L/160.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) + (1517.0L/2880.0L)*pow(x_0, 2)*x_1 - 1807.0L/15360.0L*pow(x_0, 2) - 11.0L/1800.0L*x_0*pow(x_1, 5) - 5.0L/96.0L*x_0*pow(x_1, 4) - 253.0L/1440.0L*x_0*pow(x_1, 3) - 129.0L/320.0L*x_0*pow(x_1, 2) - 2753.0L/5760.0L*x_0*x_1 + (787.0L/15360.0L)*x_0 - 11.0L/14400.0L*pow(x_1, 6) + (11.0L/2880.0L)*pow(x_1, 5) + (57.0L/1280.0L)*pow(x_1, 4) + (1711.0L/17280.0L)*pow(x_1, 3) + (1567.0L/15360.0L)*pow(x_1, 2) + (8033.0L/46080.0L)*x_1 + 18907.0L/115200.0L;
       return __pp_r74___result;
    }
    static inline O __pp_r75__(const I &x_0, const I &x_1) {
       O __pp_r75___result;
       __pp_r75___result = (23.0L/7200.0L)*pow(x_0, 6) - 37.0L/3600.0L*pow(x_0, 5)*x_1 - 7.0L/160.0L*pow(x_0, 5) + (7.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (7.0L/72.0L)*pow(x_0, 4)*x_1 + (97.0L/480.0L)*pow(x_0, 4) - 13.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/8.0L*pow(x_0, 3)*pow(x_1, 2) - 133.0L/360.0L*pow(x_0, 3)*x_1 - 1223.0L/2880.0L*pow(x_0, 3) + (1.0L/160.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) + (961.0L/1440.0L)*pow(x_0, 2)*x_1 + (871.0L/1920.0L)*pow(x_0, 2) - 11.0L/1800.0L*x_0*pow(x_1, 5) - 5.0L/96.0L*x_0*pow(x_1, 4) - 253.0L/1440.0L*x_0*pow(x_1, 3) - 129.0L/320.0L*x_0*pow(x_1, 2) - 6721.0L/11520.0L*x_0*x_1 - 1403.0L/3840.0L*x_0 - 11.0L/14400.0L*pow(x_1, 6) + (11.0L/2880.0L)*pow(x_1, 5) + (57.0L/1280.0L)*pow(x_1, 4) + (1711.0L/17280.0L)*pow(x_1, 3) + (1567.0L/15360.0L)*pow(x_1, 2) + (9491.0L/46080.0L)*x_1 + 264251.0L/921600.0L;
       return __pp_r75___result;
    }
    static inline O __pp_r76__(const I &x_0, const I &x_1) {
       O __pp_r76___result;
       __pp_r76___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 91.0L/7200.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (7.0L/360.0L)*pow(x_0, 4)*x_1 + (19.0L/288.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 17.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 7.0L/72.0L*pow(x_0, 3)*x_1 - 185.0L/1728.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (11.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (55.0L/288.0L)*pow(x_0, 2)*x_1 + (53.0L/1440.0L)*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) - 61.0L/1440.0L*x_0*pow(x_1, 4) - 31.0L/288.0L*x_0*pow(x_1, 3) - 95.0L/576.0L*x_0*pow(x_1, 2) - 1919.0L/11520.0L*x_0*x_1 - 2119.0L/28800.0L*x_0 - 7.0L/8640.0L*pow(x_1, 6) + (41.0L/14400.0L)*pow(x_1, 5) + (83.0L/2304.0L)*pow(x_1, 4) + (205.0L/3456.0L)*pow(x_1, 3) - 101.0L/46080.0L*pow(x_1, 2) + (13841.0L/230400.0L)*x_1 + 111491.0L/552960.0L;
       return __pp_r76___result;
    }
    static inline O __pp_r77__(const I &x_0, const I &x_1) {
       O __pp_r77___result;
       __pp_r77___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 31.0L/2400.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 4)*x_1 + (1.0L/36.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 3)*x_1 + (3.0L/64.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/16.0L)*pow(x_0, 2)*x_1 - 2539.0L/11520.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) - 11.0L/480.0L*x_0*pow(x_1, 4) + (23.0L/288.0L)*x_0*pow(x_1, 3) - 15.0L/64.0L*x_0*pow(x_1, 2) + (79.0L/11520.0L)*x_0*x_1 + (2477.0L/19200.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 19.0L/4800.0L*pow(x_1, 5) + (109.0L/2304.0L)*pow(x_1, 4) - 91.0L/1152.0L*pow(x_1, 3) + (1841.0L/46080.0L)*pow(x_1, 2) - 1979.0L/76800.0L*x_1 + 77293.0L/552960.0L;
       return __pp_r77___result;
    }
    static inline O __pp_r78__(const I &x_0, const I &x_1) {
       O __pp_r78___result;
       __pp_r78___result = (7.0L/4800.0L)*pow(x_0, 6) - 7.0L/3600.0L*pow(x_0, 5)*x_1 - 181.0L/14400.0L*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) + (29.0L/1440.0L)*pow(x_0, 4)*x_1 + (59.0L/2304.0L)*pow(x_0, 4) + (1.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) - 2.0L/45.0L*pow(x_0, 3)*pow(x_1, 2) - 7.0L/96.0L*pow(x_0, 3)*x_1 + (187.0L/3456.0L)*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (61.0L/576.0L)*pow(x_0, 2)*x_1 - 10781.0L/46080.0L*pow(x_0, 2) - 7.0L/3600.0L*x_0*pow(x_1, 5) + (7.0L/1440.0L)*x_0*pow(x_1, 4) + (1.0L/96.0L)*x_0*pow(x_1, 3) - 85.0L/576.0L*x_0*pow(x_1, 2) - 91.0L/1920.0L*x_0*x_1 + (32849.0L/230400.0L)*x_0 - 7.0L/4800.0L*pow(x_1, 6) + (103.0L/14400.0L)*pow(x_1, 5) + (29.0L/2304.0L)*pow(x_1, 4) - 73.0L/3456.0L*pow(x_1, 3) - 659.0L/46080.0L*pow(x_1, 2) + (313.0L/230400.0L)*x_1 + 9271.0L/69120.0L;
       return __pp_r78___result;
    }
    static inline O __pp_r79__(const I &x_0, const I &x_1) {
       O __pp_r79___result;
       __pp_r79___result = (53.0L/43200.0L)*pow(x_0, 6) - 1.0L/300.0L*pow(x_0, 5)*x_1 - 47.0L/4800.0L*pow(x_0, 5) + (1.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (49.0L/1440.0L)*pow(x_0, 4)*x_1 + (3.0L/256.0L)*pow(x_0, 4) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) - 37.0L/288.0L*pow(x_0, 3)*x_1 + (35.0L/384.0L)*pow(x_0, 3) - 11.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/45.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (125.0L/576.0L)*pow(x_0, 2)*x_1 - 4447.0L/15360.0L*pow(x_0, 2) - 1.0L/300.0L*x_0*pow(x_1, 5) + (3.0L/160.0L)*x_0*pow(x_1, 4) - 13.0L/288.0L*x_0*pow(x_1, 3) - 7.0L/192.0L*x_0*pow(x_1, 2) - 913.0L/5760.0L*x_0*x_1 + (14363.0L/76800.0L)*x_0 - 73.0L/43200.0L*pow(x_1, 6) + (143.0L/14400.0L)*pow(x_1, 5) - 1.0L/768.0L*pow(x_1, 4) + (55.0L/3456.0L)*pow(x_1, 3) - 1073.0L/15360.0L*pow(x_1, 2) + (10553.0L/230400.0L)*x_1 + 2749.0L/23040.0L;
       return __pp_r79___result;
    }
    static inline O __pp_r80__(const I &x_0, const I &x_1) {
       O __pp_r80___result;
       __pp_r80___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 - 1.0L/100.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (23.0L/720.0L)*pow(x_0, 4)*x_1 + (1.0L/80.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/40.0L*pow(x_0, 3)*pow(x_1, 2) - 11.0L/90.0L*pow(x_0, 3)*x_1 + (43.0L/480.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (299.0L/1440.0L)*pow(x_0, 2)*x_1 - 221.0L/768.0L*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (1.0L/480.0L)*x_0*pow(x_1, 4) - 29.0L/1440.0L*x_0*pow(x_1, 3) - 53.0L/960.0L*x_0*pow(x_1, 2) - 349.0L/2304.0L*x_0*x_1 + (7141.0L/38400.0L)*x_0 - 1.0L/4800.0L*pow(x_1, 6) + (47.0L/14400.0L)*pow(x_1, 5) + (43.0L/3840.0L)*pow(x_1, 4) + (59.0L/17280.0L)*pow(x_1, 3) - 193.0L/3072.0L*pow(x_1, 2) + (10067.0L/230400.0L)*x_1 + 110203.0L/921600.0L;
       return __pp_r80___result;
    }
    static inline O __pp_r81__(const I &x_0, const I &x_1) {
       O __pp_r81___result;
       __pp_r81___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 31.0L/2400.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 4)*x_1 + (1.0L/36.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 3)*x_1 + (3.0L/64.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/16.0L)*pow(x_0, 2)*x_1 - 2539.0L/11520.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) - 11.0L/480.0L*x_0*pow(x_1, 4) + (23.0L/288.0L)*x_0*pow(x_1, 3) - 15.0L/64.0L*x_0*pow(x_1, 2) + (79.0L/11520.0L)*x_0*x_1 + (2477.0L/19200.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 19.0L/4800.0L*pow(x_1, 5) + (109.0L/2304.0L)*pow(x_1, 4) - 91.0L/1152.0L*pow(x_1, 3) + (1841.0L/46080.0L)*pow(x_1, 2) - 1979.0L/76800.0L*x_1 + 77293.0L/552960.0L;
       return __pp_r81___result;
    }
    static inline O __pp_r82__(const I &x_0, const I &x_1) {
       O __pp_r82___result;
       __pp_r82___result = (7.0L/4800.0L)*pow(x_0, 6) - 7.0L/3600.0L*pow(x_0, 5)*x_1 - 181.0L/14400.0L*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) + (29.0L/1440.0L)*pow(x_0, 4)*x_1 + (59.0L/2304.0L)*pow(x_0, 4) + (1.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) - 2.0L/45.0L*pow(x_0, 3)*pow(x_1, 2) - 7.0L/96.0L*pow(x_0, 3)*x_1 + (187.0L/3456.0L)*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (61.0L/576.0L)*pow(x_0, 2)*x_1 - 10781.0L/46080.0L*pow(x_0, 2) - 7.0L/3600.0L*x_0*pow(x_1, 5) + (7.0L/1440.0L)*x_0*pow(x_1, 4) + (1.0L/96.0L)*x_0*pow(x_1, 3) - 85.0L/576.0L*x_0*pow(x_1, 2) - 91.0L/1920.0L*x_0*x_1 + (32849.0L/230400.0L)*x_0 - 7.0L/4800.0L*pow(x_1, 6) + (103.0L/14400.0L)*pow(x_1, 5) + (29.0L/2304.0L)*pow(x_1, 4) - 73.0L/3456.0L*pow(x_1, 3) - 659.0L/46080.0L*pow(x_1, 2) + (313.0L/230400.0L)*x_1 + 9271.0L/69120.0L;
       return __pp_r82___result;
    }
    static inline O __pp_r83__(const I &x_0, const I &x_1) {
       O __pp_r83___result;
       __pp_r83___result = (53.0L/43200.0L)*pow(x_0, 6) - 1.0L/300.0L*pow(x_0, 5)*x_1 - 47.0L/4800.0L*pow(x_0, 5) + (1.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (49.0L/1440.0L)*pow(x_0, 4)*x_1 + (3.0L/256.0L)*pow(x_0, 4) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) - 37.0L/288.0L*pow(x_0, 3)*x_1 + (35.0L/384.0L)*pow(x_0, 3) - 11.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/45.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (125.0L/576.0L)*pow(x_0, 2)*x_1 - 4447.0L/15360.0L*pow(x_0, 2) - 1.0L/300.0L*x_0*pow(x_1, 5) + (3.0L/160.0L)*x_0*pow(x_1, 4) - 13.0L/288.0L*x_0*pow(x_1, 3) - 7.0L/192.0L*x_0*pow(x_1, 2) - 913.0L/5760.0L*x_0*x_1 + (14363.0L/76800.0L)*x_0 - 73.0L/43200.0L*pow(x_1, 6) + (143.0L/14400.0L)*pow(x_1, 5) - 1.0L/768.0L*pow(x_1, 4) + (55.0L/3456.0L)*pow(x_1, 3) - 1073.0L/15360.0L*pow(x_1, 2) + (10553.0L/230400.0L)*x_1 + 2749.0L/23040.0L;
       return __pp_r83___result;
    }
    static inline O __pp_r84__(const I &x_0, const I &x_1) {
       O __pp_r84___result;
       __pp_r84___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 - 1.0L/100.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (23.0L/720.0L)*pow(x_0, 4)*x_1 + (1.0L/80.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/40.0L*pow(x_0, 3)*pow(x_1, 2) - 11.0L/90.0L*pow(x_0, 3)*x_1 + (43.0L/480.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (299.0L/1440.0L)*pow(x_0, 2)*x_1 - 221.0L/768.0L*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (1.0L/480.0L)*x_0*pow(x_1, 4) - 29.0L/1440.0L*x_0*pow(x_1, 3) - 53.0L/960.0L*x_0*pow(x_1, 2) - 349.0L/2304.0L*x_0*x_1 + (7141.0L/38400.0L)*x_0 - 1.0L/4800.0L*pow(x_1, 6) + (47.0L/14400.0L)*pow(x_1, 5) + (43.0L/3840.0L)*pow(x_1, 4) + (59.0L/17280.0L)*pow(x_1, 3) - 193.0L/3072.0L*pow(x_1, 2) + (10067.0L/230400.0L)*x_1 + 110203.0L/921600.0L;
       return __pp_r84___result;
    }
    static inline O __pp_r85__(const I &x_0, const I &x_1) {
       O __pp_r85___result;
       __pp_r85___result = (1.0L/5400.0L)*pow(x_0, 6) + (1.0L/1200.0L)*pow(x_0, 5)*x_1 - 13.0L/1200.0L*pow(x_0, 5) + (1.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/360.0L)*pow(x_0, 4)*x_1 + (7.0L/128.0L)*pow(x_0, 4) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) - 5.0L/144.0L*pow(x_0, 3)*x_1 - 7.0L/96.0L*pow(x_0, 3) - 11.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/45.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (11.0L/144.0L)*pow(x_0, 2)*x_1 - 131.0L/7680.0L*pow(x_0, 2) - 1.0L/300.0L*x_0*pow(x_1, 5) + (3.0L/160.0L)*x_0*pow(x_1, 4) - 13.0L/288.0L*x_0*pow(x_1, 3) - 7.0L/192.0L*x_0*pow(x_1, 2) - 611.0L/11520.0L*x_0*x_1 - 1121.0L/38400.0L*x_0 - 73.0L/43200.0L*pow(x_1, 6) + (143.0L/14400.0L)*pow(x_1, 5) - 1.0L/768.0L*pow(x_1, 4) + (55.0L/3456.0L)*pow(x_1, 3) - 1073.0L/15360.0L*pow(x_1, 2) + (3263.0L/230400.0L)*x_1 + 6877.0L/36864.0L;
       return __pp_r85___result;
    }
    static inline O __pp_r86__(const I &x_0, const I &x_1) {
       O __pp_r86___result;
       __pp_r86___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 - 53.0L/4800.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/1440.0L)*pow(x_0, 4)*x_1 + (71.0L/1280.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/40.0L*pow(x_0, 3)*pow(x_1, 2) - 41.0L/1440.0L*pow(x_0, 3)*x_1 - 143.0L/1920.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (193.0L/2880.0L)*pow(x_0, 2)*x_1 - 47.0L/3072.0L*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (1.0L/480.0L)*x_0*pow(x_1, 4) - 29.0L/1440.0L*x_0*pow(x_1, 3) - 53.0L/960.0L*x_0*pow(x_1, 2) - 53.0L/1152.0L*x_0*x_1 - 2323.0L/76800.0L*x_0 - 1.0L/4800.0L*pow(x_1, 6) + (47.0L/14400.0L)*pow(x_1, 5) + (43.0L/3840.0L)*pow(x_1, 4) + (59.0L/17280.0L)*pow(x_1, 3) - 193.0L/3072.0L*pow(x_1, 2) + (2777.0L/230400.0L)*x_1 + 21521.0L/115200.0L;
       return __pp_r86___result;
    }
    static inline O __pp_r87__(const I &x_0, const I &x_1) {
       O __pp_r87___result;
       __pp_r87___result = -1.0L/43200.0L*pow(x_0, 6) - 1.0L/3600.0L*pow(x_0, 5)*x_1 - 139.0L/14400.0L*pow(x_0, 5) - 1.0L/720.0L*pow(x_0, 4)*pow(x_1, 2) + (11.0L/1440.0L)*pow(x_0, 4)*x_1 + (599.0L/11520.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) - 61.0L/1440.0L*pow(x_0, 3)*x_1 - 1207.0L/17280.0L*pow(x_0, 3) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (29.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (233.0L/2880.0L)*pow(x_0, 2)*x_1 - 173.0L/9216.0L*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (13.0L/1440.0L)*x_0*pow(x_1, 4) - 49.0L/1440.0L*x_0*pow(x_1, 3) - 119.0L/2880.0L*x_0*pow(x_1, 2) - 61.0L/1152.0L*x_0*x_1 - 6649.0L/230400.0L*x_0 - 19.0L/43200.0L*pow(x_1, 6) + (67.0L/14400.0L)*pow(x_1, 5) + (89.0L/11520.0L)*pow(x_1, 4) + (139.0L/17280.0L)*pow(x_1, 3) - 611.0L/9216.0L*pow(x_1, 2) + (3097.0L/230400.0L)*x_1 + 64483.0L/345600.0L;
       return __pp_r87___result;
    }
    static inline O __pp_r88__(const I &x_0, const I &x_1) {
       O __pp_r88___result;
       __pp_r88___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 47.0L/4800.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/160.0L)*pow(x_0, 4)*x_1 + (601.0L/11520.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) - 59.0L/1440.0L*pow(x_0, 3)*x_1 - 403.0L/5760.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (31.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (77.0L/960.0L)*pow(x_0, 2)*x_1 - 863.0L/46080.0L*pow(x_0, 2) + (1.0L/3600.0L)*x_0*pow(x_1, 5) - 11.0L/480.0L*x_0*pow(x_1, 4) - 71.0L/1440.0L*x_0*pow(x_1, 3) - 17.0L/320.0L*x_0*pow(x_1, 2) - 319.0L/5760.0L*x_0*x_1 - 2237.0L/76800.0L*x_0 + (19.0L/43200.0L)*pow(x_1, 6) + (11.0L/4800.0L)*pow(x_1, 5) + (271.0L/11520.0L)*pow(x_1, 4) + (131.0L/5760.0L)*pow(x_1, 3) - 2657.0L/46080.0L*pow(x_1, 2) + (1201.0L/76800.0L)*x_1 + 32281.0L/172800.0L;
       return __pp_r88___result;
    }
    static inline O __pp_r89__(const I &x_0, const I &x_1) {
       O __pp_r89___result;
       __pp_r89___result = -71.0L/7200.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*x_1 + (5.0L/96.0L)*pow(x_0, 4) - 7.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 1.0L/24.0L*pow(x_0, 3)*x_1 - 121.0L/1728.0L*pow(x_0, 3) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (23.0L/288.0L)*pow(x_0, 2)*x_1 - 3.0L/160.0L*pow(x_0, 2) - 1.0L/240.0L*x_0*pow(x_1, 5) - 41.0L/1440.0L*x_0*pow(x_1, 4) - 5.0L/96.0L*x_0*pow(x_1, 3) - 31.0L/576.0L*x_0*pow(x_1, 2) - 71.0L/1280.0L*x_0*x_1 - 839.0L/28800.0L*x_0 - 1.0L/960.0L*pow(x_1, 6) + (1.0L/14400.0L)*pow(x_1, 5) + (17.0L/768.0L)*pow(x_1, 4) + (77.0L/3456.0L)*pow(x_1, 3) - 887.0L/15360.0L*pow(x_1, 2) + (3601.0L/230400.0L)*x_1 + 34433.0L/184320.0L;
       return __pp_r89___result;
    }
    static inline O __pp_r90__(const I &x_0, const I &x_1) {
       O __pp_r90___result;
       __pp_r90___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 47.0L/4800.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/160.0L)*pow(x_0, 4)*x_1 + (601.0L/11520.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) - 59.0L/1440.0L*pow(x_0, 3)*x_1 - 403.0L/5760.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (31.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (77.0L/960.0L)*pow(x_0, 2)*x_1 - 863.0L/46080.0L*pow(x_0, 2) + (31.0L/3600.0L)*x_0*pow(x_1, 5) - 1.0L/480.0L*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) - 41.0L/960.0L*x_0*pow(x_1, 2) - 19.0L/360.0L*x_0*x_1 - 739.0L/25600.0L*x_0 + (109.0L/43200.0L)*pow(x_1, 6) + (1.0L/4800.0L)*pow(x_1, 5) + (121.0L/11520.0L)*pow(x_1, 4) + (41.0L/5760.0L)*pow(x_1, 3) - 3047.0L/46080.0L*pow(x_1, 2) + (1031.0L/76800.0L)*x_1 + 257933.0L/1382400.0L;
       return __pp_r90___result;
    }
    static inline O __pp_r91__(const I &x_0, const I &x_1) {
       O __pp_r91___result;
       __pp_r91___result = -7.0L/5400.0L*pow(x_0, 6) + (19.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/3600.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/40.0L*pow(x_0, 4)*x_1 + (23.0L/1152.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 3)*x_1 - 13.0L/864.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 2)*x_1 - 1643.0L/23040.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/45.0L)*x_0*pow(x_1, 4) - 1.0L/36.0L*x_0*pow(x_1, 3) + (1.0L/144.0L)*x_0*pow(x_1, 2) + (7.0L/5760.0L)*x_0*x_1 - 119.0L/57600.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (23.0L/2400.0L)*pow(x_1, 5) - 1.0L/288.0L*pow(x_1, 4) + (5.0L/576.0L)*pow(x_1, 3) - 961.0L/11520.0L*pow(x_1, 2) + (23.0L/38400.0L)*x_1 + 10003.0L/55296.0L;
       return __pp_r91___result;
    }
    static inline O __pp_r92__(const I &x_0, const I &x_1) {
       O __pp_r92___result;
       __pp_r92___result = -11.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 + (1.0L/14400.0L)*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) - 13.0L/480.0L*pow(x_0, 4)*x_1 + (239.0L/11520.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (59.0L/1440.0L)*pow(x_0, 3)*x_1 - 287.0L/17280.0L*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/240.0L)*pow(x_0, 2)*pow(x_1, 2) - 19.0L/960.0L*pow(x_0, 2)*x_1 - 641.0L/9216.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (1.0L/180.0L)*x_0*pow(x_1, 4) - 1.0L/360.0L*x_0*pow(x_1, 3) - 17.0L/1440.0L*x_0*pow(x_1, 2) + (19.0L/2304.0L)*x_0*x_1 - 719.0L/230400.0L*x_0 - 1.0L/4320.0L*pow(x_1, 6) + (7.0L/2400.0L)*pow(x_1, 5) + (13.0L/1440.0L)*pow(x_1, 4) - 11.0L/2880.0L*pow(x_1, 3) - 11.0L/144.0L*pow(x_1, 2) - 29.0L/19200.0L*x_1 + 500879.0L/2764800.0L;
       return __pp_r92___result;
    }
    static inline O __pp_r93__(const I &x_0, const I &x_1) {
       O __pp_r93___result;
       __pp_r93___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 + (7.0L/4800.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 29.0L/1440.0L*pow(x_0, 4)*x_1 + (199.0L/11520.0L)*pow(x_0, 4) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) + (13.0L/480.0L)*pow(x_0, 3)*x_1 - 23.0L/1920.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) - 17.0L/2880.0L*pow(x_0, 2)*x_1 - 673.0L/9216.0L*pow(x_0, 2) + (1.0L/80.0L)*x_0*pow(x_1, 4) - 1.0L/60.0L*x_0*pow(x_1, 3) + (1.0L/480.0L)*x_0*pow(x_1, 2) + (1.0L/768.0L)*x_0*x_1 - 133.0L/76800.0L*x_0 - 1.0L/2160.0L*pow(x_1, 6) + (31.0L/7200.0L)*pow(x_1, 5) + (1.0L/180.0L)*pow(x_1, 4) + (7.0L/8640.0L)*pow(x_1, 3) - 23.0L/288.0L*pow(x_1, 2) - 7.0L/57600.0L*x_1 + 500239.0L/2764800.0L;
       return __pp_r93___result;
    }
    static inline O __pp_r94__(const I &x_0, const I &x_1) {
       O __pp_r94___result;
       __pp_r94___result = -7.0L/4800.0L*pow(x_0, 6) + (17.0L/3600.0L)*pow(x_0, 5)*x_1 + (19.0L/14400.0L)*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) - 31.0L/1440.0L*pow(x_0, 4)*x_1 + (67.0L/3840.0L)*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 - 209.0L/17280.0L*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) - 19.0L/2880.0L*pow(x_0, 2)*x_1 - 1121.0L/15360.0L*pow(x_0, 2) + (1.0L/1800.0L)*x_0*pow(x_1, 5) - 7.0L/360.0L*x_0*pow(x_1, 4) - 23.0L/720.0L*x_0*pow(x_1, 3) - 7.0L/720.0L*x_0*pow(x_1, 2) - 13.0L/11520.0L*x_0*x_1 - 461.0L/230400.0L*x_0 + (1.0L/2400.0L)*pow(x_1, 6) + (7.0L/3600.0L)*pow(x_1, 5) + (41.0L/1920.0L)*pow(x_1, 4) + (67.0L/4320.0L)*pow(x_1, 3) - 547.0L/7680.0L*pow(x_1, 2) + (239.0L/115200.0L)*x_1 + 166957.0L/921600.0L;
       return __pp_r94___result;
    }
    static inline O __pp_r95__(const I &x_0, const I &x_1) {
       O __pp_r95___result;
       __pp_r95___result = -1.0L/675.0L*pow(x_0, 6) + (1.0L/225.0L)*pow(x_0, 5)*x_1 + (1.0L/800.0L)*pow(x_0, 5) - 1.0L/180.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/45.0L*pow(x_0, 4)*x_1 + (5.0L/288.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/36.0L)*pow(x_0, 3)*x_1 - 7.0L/576.0L*pow(x_0, 3) - 1.0L/720.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/144.0L*pow(x_0, 2)*x_1 - 841.0L/11520.0L*pow(x_0, 2) - 7.0L/1800.0L*x_0*pow(x_1, 5) - 1.0L/40.0L*x_0*pow(x_1, 4) - 5.0L/144.0L*x_0*pow(x_1, 3) - 1.0L/96.0L*x_0*pow(x_1, 2) - 7.0L/5760.0L*x_0*x_1 - 77.0L/38400.0L*x_0 - 23.0L/21600.0L*pow(x_1, 6) - 1.0L/3600.0L*pow(x_1, 5) + (23.0L/1152.0L)*pow(x_1, 4) + (13.0L/864.0L)*pow(x_1, 3) - 1643.0L/23040.0L*pow(x_1, 2) + (119.0L/57600.0L)*x_1 + 50087.0L/276480.0L;
       return __pp_r95___result;
    }
    static inline O __pp_r96__(const I &x_0, const I &x_1) {
       O __pp_r96___result;
       __pp_r96___result = -7.0L/4800.0L*pow(x_0, 6) + (17.0L/3600.0L)*pow(x_0, 5)*x_1 + (19.0L/14400.0L)*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) - 31.0L/1440.0L*pow(x_0, 4)*x_1 + (67.0L/3840.0L)*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 - 209.0L/17280.0L*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) - 19.0L/2880.0L*pow(x_0, 2)*x_1 - 1121.0L/15360.0L*pow(x_0, 2) + (2.0L/225.0L)*x_0*pow(x_1, 5) + (1.0L/720.0L)*x_0*pow(x_1, 4) - 1.0L/90.0L*x_0*pow(x_1, 3) + (1.0L/1440.0L)*x_0*pow(x_1, 2) + (17.0L/11520.0L)*x_0*x_1 - 401.0L/230400.0L*x_0 + (1.0L/400.0L)*pow(x_1, 6) - 1.0L/7200.0L*pow(x_1, 5) + (1.0L/120.0L)*pow(x_1, 4) - 1.0L/8640.0L*pow(x_1, 3) - 51.0L/640.0L*pow(x_1, 2) - 1.0L/7200.0L*x_1 + 166747.0L/921600.0L;
       return __pp_r96___result;
    }
    static inline O __pp_r97__(const I &x_0, const I &x_1) {
       O __pp_r97___result;
       __pp_r97___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 - 29.0L/4800.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 37.0L/1440.0L*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 - 289.0L/640.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 2803.0L/2880.0L*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (1.0L/160.0L)*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) - 337.0L/960.0L*x_0*pow(x_1, 2) + (11327.0L/5760.0L)*x_0*x_1 - 211999.0L/76800.0L*x_0 - 1.0L/4800.0L*pow(x_1, 6) - 29.0L/14400.0L*pow(x_1, 5) - 13.0L/3840.0L*pow(x_1, 4) + (691.0L/17280.0L)*pow(x_1, 3) + (4027.0L/15360.0L)*pow(x_1, 2) - 353399.0L/230400.0L*x_1 + 38279.0L/19200.0L;
       return __pp_r97___result;
    }
    static inline O __pp_r98__(const I &x_0, const I &x_1) {
       O __pp_r98___result;
       __pp_r98___result = (1.0L/4320.0L)*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 - 1.0L/144.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/144.0L*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (25.0L/72.0L)*pow(x_0, 3)*x_1 - 125.0L/216.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 125.0L/72.0L*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 5.0L/144.0L*x_0*pow(x_1, 4) + (25.0L/72.0L)*x_0*pow(x_1, 3) - 125.0L/72.0L*x_0*pow(x_1, 2) + (625.0L/144.0L)*x_0*x_1 - 625.0L/144.0L*x_0 + (1.0L/4320.0L)*pow(x_1, 6) - 1.0L/144.0L*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) - 125.0L/216.0L*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) - 625.0L/144.0L*x_1 + 3125.0L/864.0L;
       return __pp_r98___result;
    }
    static inline O __pp_r99__(const I &x_0, const I &x_1) {
       O __pp_r99___result;
       __pp_r99___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 - 29.0L/4800.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 37.0L/1440.0L*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 - 289.0L/640.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 2803.0L/2880.0L*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (3.0L/80.0L)*x_0*pow(x_1, 4) - 11.0L/90.0L*x_0*pow(x_1, 3) - 101.0L/480.0L*x_0*pow(x_1, 2) + (21439.0L/11520.0L)*x_0*x_1 - 209569.0L/76800.0L*x_0 - 1.0L/800.0L*pow(x_1, 6) + (79.0L/3600.0L)*pow(x_1, 5) - 71.0L/480.0L*pow(x_1, 4) + (947.0L/2160.0L)*pow(x_1, 3) - 1187.0L/3840.0L*pow(x_1, 2) - 128707.0L/115200.0L*x_1 + 574799.0L/307200.0L;
       return __pp_r99___result;
    }
    static inline O __pp_r100__(const I &x_0, const I &x_1) {
       O __pp_r100___result;
       __pp_r100___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 - 29.0L/4800.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 37.0L/1440.0L*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 - 289.0L/640.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 2803.0L/2880.0L*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (1.0L/160.0L)*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) - 337.0L/960.0L*x_0*pow(x_1, 2) + (11327.0L/5760.0L)*x_0*x_1 - 211999.0L/76800.0L*x_0 - 1.0L/4800.0L*pow(x_1, 6) - 29.0L/14400.0L*pow(x_1, 5) - 13.0L/3840.0L*pow(x_1, 4) + (691.0L/17280.0L)*pow(x_1, 3) + (4027.0L/15360.0L)*pow(x_1, 2) - 353399.0L/230400.0L*x_1 + 38279.0L/19200.0L;
       return __pp_r100___result;
    }
    static inline O __pp_r101__(const I &x_0, const I &x_1) {
       O __pp_r101___result;
       __pp_r101___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (73.0L/14400.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (43.0L/1440.0L)*pow(x_0, 4)*x_1 - 449.0L/11520.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (7.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) - 103.0L/480.0L*pow(x_0, 3)*x_1 + (2437.0L/17280.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (41.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (2317.0L/2880.0L)*pow(x_0, 2)*x_1 - 10481.0L/46080.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) + (89.0L/1440.0L)*x_0*pow(x_1, 4) - 227.0L/480.0L*x_0*pow(x_1, 3) + (4109.0L/2880.0L)*x_0*pow(x_1, 2) - 1017.0L/640.0L*x_0*x_1 + (19363.0L/230400.0L)*x_0 - 29.0L/43200.0L*pow(x_1, 6) + (131.0L/14400.0L)*pow(x_1, 5) - 1319.0L/11520.0L*pow(x_1, 4) + (10931.0L/17280.0L)*pow(x_1, 3) - 69839.0L/46080.0L*pow(x_1, 2) + (301961.0L/230400.0L)*x_1 + 16831.0L/172800.0L;
       return __pp_r101___result;
    }
    static inline O __pp_r102__(const I &x_0, const I &x_1) {
       O __pp_r102___result;
       __pp_r102___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (11.0L/2880.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (5.0L/288.0L)*pow(x_0, 4)*x_1 - 287.0L/11520.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/36.0L)*pow(x_0, 3)*pow(x_1, 2) - 49.0L/480.0L*pow(x_0, 3)*x_1 + (979.0L/17280.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (859.0L/2880.0L)*pow(x_0, 2)*x_1 + (2641.0L/46080.0L)*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) - 11.0L/288.0L*x_0*pow(x_1, 4) - 11.0L/480.0L*x_0*pow(x_1, 3) + (1193.0L/2880.0L)*x_0*pow(x_1, 2) - 9.0L/20.0L*x_0*x_1 - 19747.0L/46080.0L*x_0 + (11.0L/4800.0L)*pow(x_1, 6) - 89.0L/2880.0L*pow(x_1, 5) + (1273.0L/11520.0L)*pow(x_1, 4) - 733.0L/17280.0L*pow(x_1, 3) - 17351.0L/46080.0L*pow(x_1, 2) + (13153.0L/46080.0L)*x_1 + 666089.0L/1382400.0L;
       return __pp_r102___result;
    }
    static inline O __pp_r103__(const I &x_0, const I &x_1) {
       O __pp_r103___result;
       __pp_r103___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 - 29.0L/4800.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 37.0L/1440.0L*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 - 289.0L/640.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 2803.0L/2880.0L*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (1.0L/160.0L)*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) - 337.0L/960.0L*x_0*pow(x_1, 2) + (11327.0L/5760.0L)*x_0*x_1 - 211999.0L/76800.0L*x_0 - 1.0L/4800.0L*pow(x_1, 6) - 29.0L/14400.0L*pow(x_1, 5) - 13.0L/3840.0L*pow(x_1, 4) + (691.0L/17280.0L)*pow(x_1, 3) + (4027.0L/15360.0L)*pow(x_1, 2) - 353399.0L/230400.0L*x_1 + 38279.0L/19200.0L;
       return __pp_r103___result;
    }
    static inline O __pp_r104__(const I &x_0, const I &x_1) {
       O __pp_r104___result;
       __pp_r104___result = (1.0L/4320.0L)*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 - 1.0L/144.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/144.0L*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (25.0L/72.0L)*pow(x_0, 3)*x_1 - 125.0L/216.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 125.0L/72.0L*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 5.0L/144.0L*x_0*pow(x_1, 4) + (25.0L/72.0L)*x_0*pow(x_1, 3) - 125.0L/72.0L*x_0*pow(x_1, 2) + (625.0L/144.0L)*x_0*x_1 - 625.0L/144.0L*x_0 + (1.0L/4320.0L)*pow(x_1, 6) - 1.0L/144.0L*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) - 125.0L/216.0L*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) - 625.0L/144.0L*x_1 + 3125.0L/864.0L;
       return __pp_r104___result;
    }
    static inline O __pp_r105__(const I &x_0, const I &x_1) {
       O __pp_r105___result;
       __pp_r105___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 - 29.0L/4800.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 37.0L/1440.0L*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 - 289.0L/640.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 2803.0L/2880.0L*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (3.0L/80.0L)*x_0*pow(x_1, 4) - 11.0L/90.0L*x_0*pow(x_1, 3) - 101.0L/480.0L*x_0*pow(x_1, 2) + (21439.0L/11520.0L)*x_0*x_1 - 209569.0L/76800.0L*x_0 - 1.0L/800.0L*pow(x_1, 6) + (79.0L/3600.0L)*pow(x_1, 5) - 71.0L/480.0L*pow(x_1, 4) + (947.0L/2160.0L)*pow(x_1, 3) - 1187.0L/3840.0L*pow(x_1, 2) - 128707.0L/115200.0L*x_1 + 574799.0L/307200.0L;
       return __pp_r105___result;
    }
    static inline O __pp_r106__(const I &x_0, const I &x_1) {
       O __pp_r106___result;
       __pp_r106___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 - 29.0L/4800.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 37.0L/1440.0L*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 - 289.0L/640.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 2803.0L/2880.0L*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (1.0L/160.0L)*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) - 337.0L/960.0L*x_0*pow(x_1, 2) + (11327.0L/5760.0L)*x_0*x_1 - 211999.0L/76800.0L*x_0 - 1.0L/4800.0L*pow(x_1, 6) - 29.0L/14400.0L*pow(x_1, 5) - 13.0L/3840.0L*pow(x_1, 4) + (691.0L/17280.0L)*pow(x_1, 3) + (4027.0L/15360.0L)*pow(x_1, 2) - 353399.0L/230400.0L*x_1 + 38279.0L/19200.0L;
       return __pp_r106___result;
    }
    static inline O __pp_r107__(const I &x_0, const I &x_1) {
       O __pp_r107___result;
       __pp_r107___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (73.0L/14400.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (43.0L/1440.0L)*pow(x_0, 4)*x_1 - 449.0L/11520.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (7.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) - 103.0L/480.0L*pow(x_0, 3)*x_1 + (2437.0L/17280.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (41.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (2317.0L/2880.0L)*pow(x_0, 2)*x_1 - 10481.0L/46080.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) + (89.0L/1440.0L)*x_0*pow(x_1, 4) - 227.0L/480.0L*x_0*pow(x_1, 3) + (4109.0L/2880.0L)*x_0*pow(x_1, 2) - 1017.0L/640.0L*x_0*x_1 + (19363.0L/230400.0L)*x_0 - 29.0L/43200.0L*pow(x_1, 6) + (131.0L/14400.0L)*pow(x_1, 5) - 1319.0L/11520.0L*pow(x_1, 4) + (10931.0L/17280.0L)*pow(x_1, 3) - 69839.0L/46080.0L*pow(x_1, 2) + (301961.0L/230400.0L)*x_1 + 16831.0L/172800.0L;
       return __pp_r107___result;
    }
    static inline O __pp_r108__(const I &x_0, const I &x_1) {
       O __pp_r108___result;
       __pp_r108___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (11.0L/2880.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (5.0L/288.0L)*pow(x_0, 4)*x_1 - 287.0L/11520.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/36.0L)*pow(x_0, 3)*pow(x_1, 2) - 49.0L/480.0L*pow(x_0, 3)*x_1 + (979.0L/17280.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (859.0L/2880.0L)*pow(x_0, 2)*x_1 + (2641.0L/46080.0L)*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) - 11.0L/288.0L*x_0*pow(x_1, 4) - 11.0L/480.0L*x_0*pow(x_1, 3) + (1193.0L/2880.0L)*x_0*pow(x_1, 2) - 9.0L/20.0L*x_0*x_1 - 19747.0L/46080.0L*x_0 + (11.0L/4800.0L)*pow(x_1, 6) - 89.0L/2880.0L*pow(x_1, 5) + (1273.0L/11520.0L)*pow(x_1, 4) - 733.0L/17280.0L*pow(x_1, 3) - 17351.0L/46080.0L*pow(x_1, 2) + (13153.0L/46080.0L)*x_1 + 666089.0L/1382400.0L;
       return __pp_r108___result;
    }
    static inline O __pp_r109__(const I &x_0, const I &x_1) {
       O __pp_r109___result;
       __pp_r109___result = (17.0L/21600.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 - 1.0L/7200.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (59.0L/720.0L)*pow(x_0, 4)*x_1 - 103.0L/1440.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (7.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) - 19.0L/40.0L*pow(x_0, 3)*x_1 + (4031.0L/8640.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (41.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (131.0L/90.0L)*pow(x_0, 2)*x_1 - 14339.0L/11520.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) + (89.0L/1440.0L)*x_0*pow(x_1, 4) - 227.0L/480.0L*x_0*pow(x_1, 3) + (4109.0L/2880.0L)*x_0*pow(x_1, 2) - 9227.0L/3840.0L*x_0*x_1 + (10859.0L/7200.0L)*x_0 - 29.0L/43200.0L*pow(x_1, 6) + (131.0L/14400.0L)*pow(x_1, 5) - 1319.0L/11520.0L*pow(x_1, 4) + (10931.0L/17280.0L)*pow(x_1, 3) - 69839.0L/46080.0L*pow(x_1, 2) + (395711.0L/230400.0L)*x_1 - 1840079.0L/2764800.0L;
       return __pp_r109___result;
    }
    static inline O __pp_r110__(const I &x_0, const I &x_1) {
       O __pp_r110___result;
       __pp_r110___result = (1.0L/1200.0L)*pow(x_0, 6) - 19.0L/3600.0L*pow(x_0, 5)*x_1 - 1.0L/720.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (5.0L/72.0L)*pow(x_0, 4)*x_1 - 331.0L/5760.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/36.0L)*pow(x_0, 3)*pow(x_1, 2) - 29.0L/80.0L*pow(x_0, 3)*x_1 + (1651.0L/4320.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (1367.0L/1440.0L)*pow(x_0, 2)*x_1 - 22117.0L/23040.0L*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) - 11.0L/288.0L*x_0*pow(x_1, 4) - 11.0L/480.0L*x_0*pow(x_1, 3) + (1193.0L/2880.0L)*x_0*pow(x_1, 2) - 4853.0L/3840.0L*x_0*x_1 + (22939.0L/23040.0L)*x_0 + (11.0L/4800.0L)*pow(x_1, 6) - 89.0L/2880.0L*pow(x_1, 5) + (1273.0L/11520.0L)*pow(x_1, 4) - 733.0L/17280.0L*pow(x_1, 3) - 17351.0L/46080.0L*pow(x_1, 2) + (31903.0L/46080.0L)*x_1 - 777197.0L/2764800.0L;
       return __pp_r110___result;
    }
    static inline O __pp_r111__(const I &x_0, const I &x_1) {
       O __pp_r111___result;
       __pp_r111___result = (11.0L/7200.0L)*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 1.0L/72.0L*pow(x_0, 5) + (1.0L/120.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/144.0L)*pow(x_0, 4)*x_1 + (209.0L/5760.0L)*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/80.0L)*pow(x_0, 3)*x_1 + (31.0L/4320.0L)*pow(x_0, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/9.0L*pow(x_0, 2)*pow(x_1, 3) + (193.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 253.0L/1440.0L*pow(x_0, 2)*x_1 - 2677.0L/23040.0L*pow(x_0, 2) + (41.0L/3600.0L)*x_0*pow(x_1, 5) - 29.0L/288.0L*x_0*pow(x_1, 4) + (169.0L/480.0L)*x_0*pow(x_1, 3) - 2047.0L/2880.0L*x_0*pow(x_1, 2) + (1627.0L/3840.0L)*x_0*x_1 - 389.0L/23040.0L*x_0 + (43.0L/14400.0L)*pow(x_1, 6) - 25.0L/576.0L*pow(x_1, 5) + (2353.0L/11520.0L)*pow(x_1, 4) - 7213.0L/17280.0L*pow(x_1, 3) + (21529.0L/46080.0L)*pow(x_1, 2) - 14753.0L/46080.0L*x_1 + 622483.0L/2764800.0L;
       return __pp_r111___result;
    }
    static inline O __pp_r112__(const I &x_0, const I &x_1) {
       O __pp_r112___result;
       __pp_r112___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 31.0L/2400.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 4)*x_1 + (1.0L/36.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 3)*x_1 + (3.0L/64.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/16.0L)*pow(x_0, 2)*x_1 - 2539.0L/11520.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) - 11.0L/480.0L*x_0*pow(x_1, 4) + (23.0L/288.0L)*x_0*pow(x_1, 3) - 15.0L/64.0L*x_0*pow(x_1, 2) + (79.0L/11520.0L)*x_0*x_1 + (2477.0L/19200.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 59.0L/4800.0L*pow(x_1, 5) + (157.0L/2304.0L)*pow(x_1, 4) - 115.0L/1152.0L*pow(x_1, 3) + (2321.0L/46080.0L)*pow(x_1, 2) - 2179.0L/76800.0L*x_1 + 77437.0L/552960.0L;
       return __pp_r112___result;
    }
    static inline O __pp_r113__(const I &x_0, const I &x_1) {
       O __pp_r113___result;
       __pp_r113___result = (17.0L/21600.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 - 1.0L/7200.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (59.0L/720.0L)*pow(x_0, 4)*x_1 - 103.0L/1440.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (7.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) - 19.0L/40.0L*pow(x_0, 3)*x_1 + (4031.0L/8640.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (41.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (131.0L/90.0L)*pow(x_0, 2)*x_1 - 14339.0L/11520.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) + (89.0L/1440.0L)*x_0*pow(x_1, 4) - 227.0L/480.0L*x_0*pow(x_1, 3) + (4109.0L/2880.0L)*x_0*pow(x_1, 2) - 9227.0L/3840.0L*x_0*x_1 + (10859.0L/7200.0L)*x_0 - 29.0L/43200.0L*pow(x_1, 6) + (131.0L/14400.0L)*pow(x_1, 5) - 1319.0L/11520.0L*pow(x_1, 4) + (10931.0L/17280.0L)*pow(x_1, 3) - 69839.0L/46080.0L*pow(x_1, 2) + (395711.0L/230400.0L)*x_1 - 1840079.0L/2764800.0L;
       return __pp_r113___result;
    }
    static inline O __pp_r114__(const I &x_0, const I &x_1) {
       O __pp_r114___result;
       __pp_r114___result = (1.0L/1200.0L)*pow(x_0, 6) - 19.0L/3600.0L*pow(x_0, 5)*x_1 - 1.0L/720.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (5.0L/72.0L)*pow(x_0, 4)*x_1 - 331.0L/5760.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/36.0L)*pow(x_0, 3)*pow(x_1, 2) - 29.0L/80.0L*pow(x_0, 3)*x_1 + (1651.0L/4320.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (1367.0L/1440.0L)*pow(x_0, 2)*x_1 - 22117.0L/23040.0L*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) - 11.0L/288.0L*x_0*pow(x_1, 4) - 11.0L/480.0L*x_0*pow(x_1, 3) + (1193.0L/2880.0L)*x_0*pow(x_1, 2) - 4853.0L/3840.0L*x_0*x_1 + (22939.0L/23040.0L)*x_0 + (11.0L/4800.0L)*pow(x_1, 6) - 89.0L/2880.0L*pow(x_1, 5) + (1273.0L/11520.0L)*pow(x_1, 4) - 733.0L/17280.0L*pow(x_1, 3) - 17351.0L/46080.0L*pow(x_1, 2) + (31903.0L/46080.0L)*x_1 - 777197.0L/2764800.0L;
       return __pp_r114___result;
    }
    static inline O __pp_r115__(const I &x_0, const I &x_1) {
       O __pp_r115___result;
       __pp_r115___result = (11.0L/7200.0L)*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 1.0L/72.0L*pow(x_0, 5) + (1.0L/120.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/144.0L)*pow(x_0, 4)*x_1 + (209.0L/5760.0L)*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/80.0L)*pow(x_0, 3)*x_1 + (31.0L/4320.0L)*pow(x_0, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/9.0L*pow(x_0, 2)*pow(x_1, 3) + (193.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 253.0L/1440.0L*pow(x_0, 2)*x_1 - 2677.0L/23040.0L*pow(x_0, 2) + (41.0L/3600.0L)*x_0*pow(x_1, 5) - 29.0L/288.0L*x_0*pow(x_1, 4) + (169.0L/480.0L)*x_0*pow(x_1, 3) - 2047.0L/2880.0L*x_0*pow(x_1, 2) + (1627.0L/3840.0L)*x_0*x_1 - 389.0L/23040.0L*x_0 + (43.0L/14400.0L)*pow(x_1, 6) - 25.0L/576.0L*pow(x_1, 5) + (2353.0L/11520.0L)*pow(x_1, 4) - 7213.0L/17280.0L*pow(x_1, 3) + (21529.0L/46080.0L)*pow(x_1, 2) - 14753.0L/46080.0L*x_1 + 622483.0L/2764800.0L;
       return __pp_r115___result;
    }
    static inline O __pp_r116__(const I &x_0, const I &x_1) {
       O __pp_r116___result;
       __pp_r116___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 31.0L/2400.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 4)*x_1 + (1.0L/36.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 3)*x_1 + (3.0L/64.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/16.0L)*pow(x_0, 2)*x_1 - 2539.0L/11520.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) - 11.0L/480.0L*x_0*pow(x_1, 4) + (23.0L/288.0L)*x_0*pow(x_1, 3) - 15.0L/64.0L*x_0*pow(x_1, 2) + (79.0L/11520.0L)*x_0*x_1 + (2477.0L/19200.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 59.0L/4800.0L*pow(x_1, 5) + (157.0L/2304.0L)*pow(x_1, 4) - 115.0L/1152.0L*pow(x_1, 3) + (2321.0L/46080.0L)*pow(x_1, 2) - 2179.0L/76800.0L*x_1 + 77437.0L/552960.0L;
       return __pp_r116___result;
    }
    static inline O __pp_r117__(const I &x_0, const I &x_1) {
       O __pp_r117___result;
       __pp_r117___result = (1.0L/4320.0L)*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 - 1.0L/144.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/144.0L*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (25.0L/72.0L)*pow(x_0, 3)*x_1 - 125.0L/216.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 125.0L/72.0L*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 5.0L/144.0L*x_0*pow(x_1, 4) + (25.0L/72.0L)*x_0*pow(x_1, 3) - 125.0L/72.0L*x_0*pow(x_1, 2) + (625.0L/144.0L)*x_0*x_1 - 625.0L/144.0L*x_0 + (1.0L/4320.0L)*pow(x_1, 6) - 1.0L/144.0L*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) - 125.0L/216.0L*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) - 625.0L/144.0L*x_1 + 3125.0L/864.0L;
       return __pp_r117___result;
    }
    static inline O __pp_r118__(const I &x_0, const I &x_1) {
       O __pp_r118___result;
       __pp_r118___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 - 29.0L/4800.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 37.0L/1440.0L*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 - 289.0L/640.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 2803.0L/2880.0L*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (3.0L/80.0L)*x_0*pow(x_1, 4) - 11.0L/90.0L*x_0*pow(x_1, 3) - 101.0L/480.0L*x_0*pow(x_1, 2) + (21439.0L/11520.0L)*x_0*x_1 - 209569.0L/76800.0L*x_0 - 1.0L/800.0L*pow(x_1, 6) + (79.0L/3600.0L)*pow(x_1, 5) - 71.0L/480.0L*pow(x_1, 4) + (947.0L/2160.0L)*pow(x_1, 3) - 1187.0L/3840.0L*pow(x_1, 2) - 128707.0L/115200.0L*x_1 + 574799.0L/307200.0L;
       return __pp_r118___result;
    }
    static inline O __pp_r119__(const I &x_0, const I &x_1) {
       O __pp_r119___result;
       __pp_r119___result = (1.0L/4320.0L)*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 - 1.0L/144.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/144.0L*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (25.0L/72.0L)*pow(x_0, 3)*x_1 - 125.0L/216.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 125.0L/72.0L*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 5.0L/144.0L*x_0*pow(x_1, 4) + (25.0L/72.0L)*x_0*pow(x_1, 3) - 125.0L/72.0L*x_0*pow(x_1, 2) + (625.0L/144.0L)*x_0*x_1 - 625.0L/144.0L*x_0 + (1.0L/4320.0L)*pow(x_1, 6) - 1.0L/144.0L*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) - 125.0L/216.0L*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) - 625.0L/144.0L*x_1 + 3125.0L/864.0L;
       return __pp_r119___result;
    }
    static inline O __pp_r120__(const I &x_0, const I &x_1) {
       O __pp_r120___result;
       __pp_r120___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 - 29.0L/4800.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 37.0L/1440.0L*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 - 289.0L/640.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 2803.0L/2880.0L*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (3.0L/80.0L)*x_0*pow(x_1, 4) - 11.0L/90.0L*x_0*pow(x_1, 3) - 101.0L/480.0L*x_0*pow(x_1, 2) + (21439.0L/11520.0L)*x_0*x_1 - 209569.0L/76800.0L*x_0 - 1.0L/800.0L*pow(x_1, 6) + (79.0L/3600.0L)*pow(x_1, 5) - 71.0L/480.0L*pow(x_1, 4) + (947.0L/2160.0L)*pow(x_1, 3) - 1187.0L/3840.0L*pow(x_1, 2) - 128707.0L/115200.0L*x_1 + 574799.0L/307200.0L;
       return __pp_r120___result;
    }
    static inline O __pp_r121__(const I &x_0, const I &x_1) {
       O __pp_r121___result;
       __pp_r121___result = (11.0L/8640.0L)*pow(x_0, 6) - 1.0L/360.0L*pow(x_0, 5)*x_1 - 7.0L/576.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/288.0L)*pow(x_0, 4)*x_1 + (125.0L/2304.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (25.0L/288.0L)*pow(x_0, 3)*x_1 - 875.0L/3456.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 625.0L/576.0L*pow(x_0, 2)*x_1 + (10625.0L/9216.0L)*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 5.0L/144.0L*x_0*pow(x_1, 4) + (25.0L/72.0L)*x_0*pow(x_1, 3) - 125.0L/72.0L*x_0*pow(x_1, 2) + (8125.0L/2304.0L)*x_0*x_1 - 26875.0L/9216.0L*x_0 + (1.0L/4320.0L)*pow(x_1, 6) - 1.0L/144.0L*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) - 125.0L/216.0L*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) - 18125.0L/4608.0L*x_1 + 315625.0L/110592.0L;
       return __pp_r121___result;
    }
    static inline O __pp_r122__(const I &x_0, const I &x_1) {
       O __pp_r122___result;
       __pp_r122___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 - 9.0L/800.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (19.0L/720.0L)*pow(x_0, 4)*x_1 + (19.0L/480.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) - 11.0L/360.0L*pow(x_0, 3)*x_1 - 121.0L/960.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 29.0L/90.0L*pow(x_0, 2)*x_1 + (2047.0L/3840.0L)*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (3.0L/80.0L)*x_0*pow(x_1, 4) - 11.0L/90.0L*x_0*pow(x_1, 3) - 101.0L/480.0L*x_0*pow(x_1, 2) + (377.0L/360.0L)*x_0*x_1 - 16699.0L/12800.0L*x_0 - 1.0L/800.0L*pow(x_1, 6) + (79.0L/3600.0L)*pow(x_1, 5) - 71.0L/480.0L*pow(x_1, 4) + (947.0L/2160.0L)*pow(x_1, 3) - 1187.0L/3840.0L*pow(x_1, 2) - 10229.0L/14400.0L*x_1 + 42553.0L/38400.0L;
       return __pp_r122___result;
    }
    static inline O __pp_r123__(const I &x_0, const I &x_1) {
       O __pp_r123___result;
       __pp_r123___result = (17.0L/21600.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 - 1.0L/7200.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (59.0L/720.0L)*pow(x_0, 4)*x_1 - 103.0L/1440.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (7.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) - 19.0L/40.0L*pow(x_0, 3)*x_1 + (4031.0L/8640.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (41.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (131.0L/90.0L)*pow(x_0, 2)*x_1 - 14339.0L/11520.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) + (67.0L/720.0L)*x_0*pow(x_1, 4) - 17.0L/30.0L*x_0*pow(x_1, 3) + (2257.0L/1440.0L)*x_0*pow(x_1, 2) - 301.0L/120.0L*x_0*x_1 + (177389.0L/115200.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (119.0L/3600.0L)*pow(x_1, 5) - 373.0L/1440.0L*pow(x_1, 4) + (2227.0L/2160.0L)*pow(x_1, 3) - 24041.0L/11520.0L*pow(x_1, 2) + (30731.0L/14400.0L)*x_1 - 272383.0L/345600.0L;
       return __pp_r123___result;
    }
    static inline O __pp_r124__(const I &x_0, const I &x_1) {
       O __pp_r124___result;
       __pp_r124___result = -1.0L/4800.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 - 1.0L/960.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 - 25.0L/128.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/12.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/64.0L*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 1.0L/32.0L*x_0*pow(x_1, 4) + (35.0L/96.0L)*x_0*pow(x_1, 3) - 325.0L/192.0L*x_0*pow(x_1, 2) + (1375.0L/384.0L)*x_0*x_1 - 8875.0L/3072.0L*x_0 + (1.0L/4800.0L)*pow(x_1, 6) - 7.0L/960.0L*pow(x_1, 5) + (65.0L/768.0L)*pow(x_1, 4) - 75.0L/128.0L*pow(x_1, 3) + (6625.0L/3072.0L)*pow(x_1, 2) - 12125.0L/3072.0L*x_1 + 4375.0L/1536.0L;
       return __pp_r124___result;
    }
    static inline O __pp_r125__(const I &x_0, const I &x_1) {
       O __pp_r125___result;
       __pp_r125___result = -1.0L/4320.0L*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 - 1.0L/7200.0L*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/720.0L*pow(x_0, 4)*x_1 + (7.0L/1440.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (7.0L/180.0L)*pow(x_0, 3)*x_1 - 589.0L/8640.0L*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) - 589.0L/1440.0L*pow(x_0, 2)*x_1 + (1379.0L/2880.0L)*pow(x_0, 2) - 1.0L/360.0L*x_0*pow(x_1, 5) + (59.0L/1440.0L)*x_0*pow(x_1, 4) - 151.0L/1440.0L*x_0*pow(x_1, 3) - 481.0L/2880.0L*x_0*pow(x_1, 2) + (12689.0L/11520.0L)*x_0*x_1 - 73583.0L/57600.0L*x_0 - 11.0L/8640.0L*pow(x_1, 6) + (311.0L/14400.0L)*pow(x_1, 5) - 1729.0L/11520.0L*pow(x_1, 4) + (7451.0L/17280.0L)*pow(x_1, 3) - 14869.0L/46080.0L*pow(x_1, 2) - 166789.0L/230400.0L*x_1 + 3048191.0L/2764800.0L;
       return __pp_r125___result;
    }
    static inline O __pp_r126__(const I &x_0, const I &x_1) {
       O __pp_r126___result;
       __pp_r126___result = -1.0L/1440.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (79.0L/7200.0L)*pow(x_0, 5) - 1.0L/96.0L*pow(x_0, 4)*pow(x_1, 2) + (13.0L/240.0L)*pow(x_0, 4)*x_1 - 17.0L/160.0L*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (19.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) - 73.0L/180.0L*pow(x_0, 3)*x_1 + (4531.0L/8640.0L)*pow(x_0, 3) - 1.0L/96.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/10.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (219.0L/160.0L)*pow(x_0, 2)*x_1 - 1247.0L/960.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) + (139.0L/1440.0L)*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) + (4639.0L/2880.0L)*x_0*pow(x_1, 2) - 28271.0L/11520.0L*x_0*x_1 + (90257.0L/57600.0L)*x_0 - 1.0L/576.0L*pow(x_1, 6) + (157.0L/4800.0L)*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) + (5897.0L/5760.0L)*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) + (162857.0L/76800.0L)*x_1 - 731563.0L/921600.0L;
       return __pp_r126___result;
    }
    static inline O __pp_r127__(const I &x_0, const I &x_1) {
       O __pp_r127___result;
       __pp_r127___result = (17.0L/21600.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 - 1.0L/7200.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (59.0L/720.0L)*pow(x_0, 4)*x_1 - 103.0L/1440.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (7.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) - 19.0L/40.0L*pow(x_0, 3)*x_1 + (4031.0L/8640.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (41.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (131.0L/90.0L)*pow(x_0, 2)*x_1 - 14339.0L/11520.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) + (89.0L/1440.0L)*x_0*pow(x_1, 4) - 227.0L/480.0L*x_0*pow(x_1, 3) + (4109.0L/2880.0L)*x_0*pow(x_1, 2) - 9227.0L/3840.0L*x_0*x_1 + (10859.0L/7200.0L)*x_0 - 29.0L/43200.0L*pow(x_1, 6) + (131.0L/14400.0L)*pow(x_1, 5) - 1319.0L/11520.0L*pow(x_1, 4) + (10931.0L/17280.0L)*pow(x_1, 3) - 69839.0L/46080.0L*pow(x_1, 2) + (395711.0L/230400.0L)*x_1 - 1840079.0L/2764800.0L;
       return __pp_r127___result;
    }
    static inline O __pp_r128__(const I &x_0, const I &x_1) {
       O __pp_r128___result;
       __pp_r128___result = (1.0L/1200.0L)*pow(x_0, 6) - 19.0L/3600.0L*pow(x_0, 5)*x_1 - 1.0L/720.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (5.0L/72.0L)*pow(x_0, 4)*x_1 - 331.0L/5760.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/36.0L)*pow(x_0, 3)*pow(x_1, 2) - 29.0L/80.0L*pow(x_0, 3)*x_1 + (1651.0L/4320.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (1367.0L/1440.0L)*pow(x_0, 2)*x_1 - 22117.0L/23040.0L*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) - 11.0L/288.0L*x_0*pow(x_1, 4) - 11.0L/480.0L*x_0*pow(x_1, 3) + (1193.0L/2880.0L)*x_0*pow(x_1, 2) - 4853.0L/3840.0L*x_0*x_1 + (22939.0L/23040.0L)*x_0 + (11.0L/4800.0L)*pow(x_1, 6) - 89.0L/2880.0L*pow(x_1, 5) + (1273.0L/11520.0L)*pow(x_1, 4) - 733.0L/17280.0L*pow(x_1, 3) - 17351.0L/46080.0L*pow(x_1, 2) + (31903.0L/46080.0L)*x_1 - 777197.0L/2764800.0L;
       return __pp_r128___result;
    }
    static inline O __pp_r129__(const I &x_0, const I &x_1) {
       O __pp_r129___result;
       __pp_r129___result = (17.0L/21600.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 - 1.0L/7200.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (59.0L/720.0L)*pow(x_0, 4)*x_1 - 103.0L/1440.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (7.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) - 19.0L/40.0L*pow(x_0, 3)*x_1 + (4031.0L/8640.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (41.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (131.0L/90.0L)*pow(x_0, 2)*x_1 - 14339.0L/11520.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) + (67.0L/720.0L)*x_0*pow(x_1, 4) - 17.0L/30.0L*x_0*pow(x_1, 3) + (2257.0L/1440.0L)*x_0*pow(x_1, 2) - 301.0L/120.0L*x_0*x_1 + (177389.0L/115200.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (119.0L/3600.0L)*pow(x_1, 5) - 373.0L/1440.0L*pow(x_1, 4) + (2227.0L/2160.0L)*pow(x_1, 3) - 24041.0L/11520.0L*pow(x_1, 2) + (30731.0L/14400.0L)*x_1 - 272383.0L/345600.0L;
       return __pp_r129___result;
    }
    static inline O __pp_r130__(const I &x_0, const I &x_1) {
       O __pp_r130___result;
       __pp_r130___result = (11.0L/7200.0L)*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 1.0L/72.0L*pow(x_0, 5) + (1.0L/120.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/144.0L)*pow(x_0, 4)*x_1 + (209.0L/5760.0L)*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/80.0L)*pow(x_0, 3)*x_1 + (31.0L/4320.0L)*pow(x_0, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/9.0L*pow(x_0, 2)*pow(x_1, 3) + (193.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 253.0L/1440.0L*pow(x_0, 2)*x_1 - 2677.0L/23040.0L*pow(x_0, 2) + (41.0L/3600.0L)*x_0*pow(x_1, 5) - 29.0L/288.0L*x_0*pow(x_1, 4) + (169.0L/480.0L)*x_0*pow(x_1, 3) - 2047.0L/2880.0L*x_0*pow(x_1, 2) + (1627.0L/3840.0L)*x_0*x_1 - 389.0L/23040.0L*x_0 + (43.0L/14400.0L)*pow(x_1, 6) - 25.0L/576.0L*pow(x_1, 5) + (2353.0L/11520.0L)*pow(x_1, 4) - 7213.0L/17280.0L*pow(x_1, 3) + (21529.0L/46080.0L)*pow(x_1, 2) - 14753.0L/46080.0L*x_1 + 622483.0L/2764800.0L;
       return __pp_r130___result;
    }
    static inline O __pp_r131__(const I &x_0, const I &x_1) {
       O __pp_r131___result;
       __pp_r131___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 31.0L/2400.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 4)*x_1 + (1.0L/36.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 3)*x_1 + (3.0L/64.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/16.0L)*pow(x_0, 2)*x_1 - 2539.0L/11520.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) - 11.0L/480.0L*x_0*pow(x_1, 4) + (23.0L/288.0L)*x_0*pow(x_1, 3) - 15.0L/64.0L*x_0*pow(x_1, 2) + (79.0L/11520.0L)*x_0*x_1 + (2477.0L/19200.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 59.0L/4800.0L*pow(x_1, 5) + (157.0L/2304.0L)*pow(x_1, 4) - 115.0L/1152.0L*pow(x_1, 3) + (2321.0L/46080.0L)*pow(x_1, 2) - 2179.0L/76800.0L*x_1 + 77437.0L/552960.0L;
       return __pp_r131___result;
    }
    static inline O __pp_r132__(const I &x_0, const I &x_1) {
       O __pp_r132___result;
       __pp_r132___result = -1.0L/1440.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (79.0L/7200.0L)*pow(x_0, 5) - 1.0L/96.0L*pow(x_0, 4)*pow(x_1, 2) + (13.0L/240.0L)*pow(x_0, 4)*x_1 - 17.0L/160.0L*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (19.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) - 73.0L/180.0L*pow(x_0, 3)*x_1 + (4531.0L/8640.0L)*pow(x_0, 3) - 1.0L/96.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/10.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (219.0L/160.0L)*pow(x_0, 2)*x_1 - 1247.0L/960.0L*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) + (47.0L/720.0L)*x_0*pow(x_1, 4) - 41.0L/90.0L*x_0*pow(x_1, 3) + (2117.0L/1440.0L)*x_0*pow(x_1, 2) - 1691.0L/720.0L*x_0*x_1 + (176869.0L/115200.0L)*x_0 - 1.0L/1440.0L*pow(x_1, 6) + (7.0L/800.0L)*pow(x_1, 5) - 7.0L/60.0L*pow(x_1, 4) + (1801.0L/2880.0L)*pow(x_1, 3) - 367.0L/240.0L*pow(x_1, 2) + (65431.0L/38400.0L)*x_1 - 77321.0L/115200.0L;
       return __pp_r132___result;
    }
    static inline O __pp_r133__(const I &x_0, const I &x_1) {
       O __pp_r133___result;
       __pp_r133___result = -7.0L/10800.0L*pow(x_0, 6) - 1.0L/1200.0L*pow(x_0, 5)*x_1 + (7.0L/720.0L)*pow(x_0, 5) - 11.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/24.0L)*pow(x_0, 4)*x_1 - 59.0L/640.0L*pow(x_0, 4) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/18.0L)*pow(x_0, 3)*pow(x_1, 2) - 211.0L/720.0L*pow(x_0, 3)*x_1 + (1901.0L/4320.0L)*pow(x_0, 3) + (1.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 17.0L/80.0L*pow(x_0, 2)*pow(x_1, 2) + (69.0L/80.0L)*pow(x_0, 2)*x_1 - 7789.0L/7680.0L*pow(x_0, 2) + (3.0L/400.0L)*x_0*pow(x_1, 5) - 5.0L/144.0L*x_0*pow(x_1, 4) - 1.0L/180.0L*x_0*pow(x_1, 3) + (659.0L/1440.0L)*x_0*pow(x_1, 2) - 6967.0L/5760.0L*x_0*x_1 + (5891.0L/5760.0L)*x_0 + (49.0L/21600.0L)*pow(x_1, 6) - 1.0L/32.0L*pow(x_1, 5) + (13.0L/120.0L)*pow(x_1, 4) - 143.0L/2880.0L*pow(x_1, 3) - 749.0L/1920.0L*pow(x_1, 2) + (5213.0L/7680.0L)*x_1 - 132137.0L/460800.0L;
       return __pp_r133___result;
    }
    static inline O __pp_r134__(const I &x_0, const I &x_1) {
       O __pp_r134___result;
       __pp_r134___result = -1.0L/1440.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (79.0L/7200.0L)*pow(x_0, 5) - 1.0L/96.0L*pow(x_0, 4)*pow(x_1, 2) + (13.0L/240.0L)*pow(x_0, 4)*x_1 - 17.0L/160.0L*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (19.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) - 73.0L/180.0L*pow(x_0, 3)*x_1 + (4531.0L/8640.0L)*pow(x_0, 3) - 1.0L/96.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/10.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (219.0L/160.0L)*pow(x_0, 2)*x_1 - 1247.0L/960.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) + (139.0L/1440.0L)*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) + (4639.0L/2880.0L)*x_0*pow(x_1, 2) - 28271.0L/11520.0L*x_0*x_1 + (90257.0L/57600.0L)*x_0 - 1.0L/576.0L*pow(x_1, 6) + (157.0L/4800.0L)*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) + (5897.0L/5760.0L)*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) + (162857.0L/76800.0L)*x_1 - 731563.0L/921600.0L;
       return __pp_r134___result;
    }
    static inline O __pp_r135__(const I &x_0, const I &x_1) {
       O __pp_r135___result;
       __pp_r135___result = (1.0L/21600.0L)*pow(x_0, 6) + (1.0L/300.0L)*pow(x_0, 5)*x_1 - 1.0L/360.0L*pow(x_0, 5) + (1.0L/360.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/48.0L*pow(x_0, 4)*x_1 + (1.0L/640.0L)*pow(x_0, 4) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (59.0L/720.0L)*pow(x_0, 3)*x_1 + (281.0L/4320.0L)*pow(x_0, 3) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/8.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) - 21.0L/80.0L*pow(x_0, 2)*x_1 - 1309.0L/7680.0L*pow(x_0, 2) + (7.0L/600.0L)*x_0*pow(x_1, 5) - 7.0L/72.0L*x_0*pow(x_1, 4) + (133.0L/360.0L)*x_0*pow(x_1, 3) - 961.0L/1440.0L*x_0*pow(x_1, 2) + (2753.0L/5760.0L)*x_0*x_1 + (59.0L/5760.0L)*x_0 + (2.0L/675.0L)*pow(x_1, 6) - 7.0L/160.0L*pow(x_1, 5) + (97.0L/480.0L)*pow(x_1, 4) - 1223.0L/2880.0L*pow(x_1, 3) + (871.0L/1920.0L)*pow(x_1, 2) - 2563.0L/7680.0L*x_1 + 101143.0L/460800.0L;
       return __pp_r135___result;
    }
    static inline O __pp_r136__(const I &x_0, const I &x_1) {
       O __pp_r136___result;
       __pp_r136___result = (1.0L/360.0L)*pow(x_0, 5)*x_1 - 13.0L/7200.0L*pow(x_0, 5) - 1.0L/90.0L*pow(x_0, 4)*x_1 - 1.0L/144.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) - 11.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/72.0L)*pow(x_0, 3)*x_1 + (181.0L/1728.0L)*pow(x_0, 3) - 17.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 7.0L/288.0L*pow(x_0, 2)*x_1 - 791.0L/2880.0L*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) - 7.0L/360.0L*x_0*pow(x_1, 4) + (7.0L/72.0L)*x_0*pow(x_1, 3) - 55.0L/288.0L*x_0*pow(x_1, 2) + (11.0L/180.0L)*x_0*x_1 + (17987.0L/115200.0L)*x_0 - 91.0L/7200.0L*pow(x_1, 5) + (19.0L/288.0L)*pow(x_1, 4) - 185.0L/1728.0L*pow(x_1, 3) + (53.0L/1440.0L)*pow(x_1, 2) - 4831.0L/115200.0L*x_1 + 9289.0L/69120.0L;
       return __pp_r136___result;
    }
    static inline O __pp_r137__(const I &x_0, const I &x_1) {
       O __pp_r137___result;
       __pp_r137___result = -43.0L/43200.0L*pow(x_0, 6) + (3.0L/400.0L)*pow(x_0, 5)*x_1 - 11.0L/2880.0L*pow(x_0, 5) + (1.0L/360.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/96.0L*pow(x_0, 4)*x_1 + (57.0L/1280.0L)*pow(x_0, 4) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (253.0L/1440.0L)*pow(x_0, 3)*x_1 - 1711.0L/17280.0L*pow(x_0, 3) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/8.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) - 129.0L/320.0L*pow(x_0, 2)*x_1 + (1567.0L/15360.0L)*pow(x_0, 2) + (7.0L/600.0L)*x_0*pow(x_1, 5) - 7.0L/72.0L*x_0*pow(x_1, 4) + (133.0L/360.0L)*x_0*pow(x_1, 3) - 961.0L/1440.0L*x_0*pow(x_1, 2) + (6721.0L/11520.0L)*x_0*x_1 - 9491.0L/46080.0L*x_0 + (2.0L/675.0L)*pow(x_1, 6) - 7.0L/160.0L*pow(x_1, 5) + (97.0L/480.0L)*pow(x_1, 4) - 1223.0L/2880.0L*pow(x_1, 3) + (871.0L/1920.0L)*pow(x_1, 2) - 1403.0L/3840.0L*x_1 + 264251.0L/921600.0L;
       return __pp_r137___result;
    }
    static inline O __pp_r138__(const I &x_0, const I &x_1) {
       O __pp_r138___result;
       __pp_r138___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 - 41.0L/14400.0L*pow(x_0, 5) - 61.0L/1440.0L*pow(x_0, 4)*x_1 + (83.0L/2304.0L)*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) - 11.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (31.0L/288.0L)*pow(x_0, 3)*x_1 - 205.0L/3456.0L*pow(x_0, 3) - 17.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 95.0L/576.0L*pow(x_0, 2)*x_1 - 101.0L/46080.0L*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) - 7.0L/360.0L*x_0*pow(x_1, 4) + (7.0L/72.0L)*x_0*pow(x_1, 3) - 55.0L/288.0L*x_0*pow(x_1, 2) + (1919.0L/11520.0L)*x_0*x_1 - 13841.0L/230400.0L*x_0 - 91.0L/7200.0L*pow(x_1, 5) + (19.0L/288.0L)*pow(x_1, 4) - 185.0L/1728.0L*pow(x_1, 3) + (53.0L/1440.0L)*pow(x_1, 2) - 2119.0L/28800.0L*x_1 + 111491.0L/552960.0L;
       return __pp_r138___result;
    }
    static inline O __pp_r139__(const I &x_0, const I &x_1) {
       O __pp_r139___result;
       __pp_r139___result = -7.0L/5400.0L*pow(x_0, 6) + (19.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/3600.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/40.0L*pow(x_0, 4)*x_1 + (23.0L/1152.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 3)*x_1 - 13.0L/864.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 2)*x_1 - 1643.0L/23040.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/45.0L)*x_0*pow(x_1, 4) - 1.0L/36.0L*x_0*pow(x_1, 3) + (1.0L/144.0L)*x_0*pow(x_1, 2) + (7.0L/5760.0L)*x_0*x_1 - 119.0L/57600.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (23.0L/2400.0L)*pow(x_1, 5) - 1.0L/288.0L*pow(x_1, 4) + (5.0L/576.0L)*pow(x_1, 3) - 961.0L/11520.0L*pow(x_1, 2) + (23.0L/38400.0L)*x_1 + 10003.0L/55296.0L;
       return __pp_r139___result;
    }
    static inline O __pp_r140__(const I &x_0, const I &x_1) {
       O __pp_r140___result;
       __pp_r140___result = -11.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 - 1.0L/14400.0L*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) - 41.0L/1440.0L*pow(x_0, 4)*x_1 + (17.0L/768.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/96.0L)*pow(x_0, 3)*x_1 - 77.0L/3456.0L*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 31.0L/576.0L*pow(x_0, 2)*x_1 - 887.0L/15360.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 1.0L/180.0L*x_0*pow(x_1, 4) + (1.0L/24.0L)*x_0*pow(x_1, 3) - 23.0L/288.0L*x_0*pow(x_1, 2) + (71.0L/1280.0L)*x_0*x_1 - 3601.0L/230400.0L*x_0 - 1.0L/4320.0L*pow(x_1, 6) - 71.0L/7200.0L*pow(x_1, 5) + (5.0L/96.0L)*pow(x_1, 4) - 121.0L/1728.0L*pow(x_1, 3) - 3.0L/160.0L*pow(x_1, 2) - 839.0L/28800.0L*x_1 + 34433.0L/184320.0L;
       return __pp_r140___result;
    }
    static inline O __pp_r141__(const I &x_0, const I &x_1) {
       O __pp_r141___result;
       __pp_r141___result = -7.0L/5400.0L*pow(x_0, 6) + (19.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/3600.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/40.0L*pow(x_0, 4)*x_1 + (23.0L/1152.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 3)*x_1 - 13.0L/864.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 2)*x_1 - 1643.0L/23040.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/45.0L)*x_0*pow(x_1, 4) - 1.0L/36.0L*x_0*pow(x_1, 3) + (1.0L/144.0L)*x_0*pow(x_1, 2) + (7.0L/5760.0L)*x_0*x_1 - 119.0L/57600.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (1.0L/800.0L)*pow(x_1, 5) + (5.0L/288.0L)*pow(x_1, 4) - 7.0L/576.0L*pow(x_1, 3) - 841.0L/11520.0L*pow(x_1, 2) - 77.0L/38400.0L*x_1 + 50087.0L/276480.0L;
       return __pp_r141___result;
    }
    static inline O __pp_r142__(const I &x_0, const I &x_1) {
       O __pp_r142___result;
       __pp_r142___result = -43.0L/43200.0L*pow(x_0, 6) + (3.0L/400.0L)*pow(x_0, 5)*x_1 - 11.0L/2880.0L*pow(x_0, 5) + (1.0L/360.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/96.0L*pow(x_0, 4)*x_1 + (57.0L/1280.0L)*pow(x_0, 4) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (253.0L/1440.0L)*pow(x_0, 3)*x_1 - 1711.0L/17280.0L*pow(x_0, 3) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/8.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) - 129.0L/320.0L*pow(x_0, 2)*x_1 + (1567.0L/15360.0L)*pow(x_0, 2) + (7.0L/600.0L)*x_0*pow(x_1, 5) - 7.0L/72.0L*x_0*pow(x_1, 4) + (133.0L/360.0L)*x_0*pow(x_1, 3) - 961.0L/1440.0L*x_0*pow(x_1, 2) + (6721.0L/11520.0L)*x_0*x_1 - 9491.0L/46080.0L*x_0 + (2.0L/675.0L)*pow(x_1, 6) - 7.0L/160.0L*pow(x_1, 5) + (97.0L/480.0L)*pow(x_1, 4) - 1223.0L/2880.0L*pow(x_1, 3) + (871.0L/1920.0L)*pow(x_1, 2) - 1403.0L/3840.0L*x_1 + 264251.0L/921600.0L;
       return __pp_r142___result;
    }
    static inline O __pp_r143__(const I &x_0, const I &x_1) {
       O __pp_r143___result;
       __pp_r143___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 - 41.0L/14400.0L*pow(x_0, 5) - 61.0L/1440.0L*pow(x_0, 4)*x_1 + (83.0L/2304.0L)*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) - 11.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (31.0L/288.0L)*pow(x_0, 3)*x_1 - 205.0L/3456.0L*pow(x_0, 3) - 17.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 95.0L/576.0L*pow(x_0, 2)*x_1 - 101.0L/46080.0L*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) - 7.0L/360.0L*x_0*pow(x_1, 4) + (7.0L/72.0L)*x_0*pow(x_1, 3) - 55.0L/288.0L*x_0*pow(x_1, 2) + (1919.0L/11520.0L)*x_0*x_1 - 13841.0L/230400.0L*x_0 - 91.0L/7200.0L*pow(x_1, 5) + (19.0L/288.0L)*pow(x_1, 4) - 185.0L/1728.0L*pow(x_1, 3) + (53.0L/1440.0L)*pow(x_1, 2) - 2119.0L/28800.0L*x_1 + 111491.0L/552960.0L;
       return __pp_r143___result;
    }
    static inline O __pp_r144__(const I &x_0, const I &x_1) {
       O __pp_r144___result;
       __pp_r144___result = -7.0L/5400.0L*pow(x_0, 6) + (19.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/3600.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/40.0L*pow(x_0, 4)*x_1 + (23.0L/1152.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 3)*x_1 - 13.0L/864.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 2)*x_1 - 1643.0L/23040.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/45.0L)*x_0*pow(x_1, 4) - 1.0L/36.0L*x_0*pow(x_1, 3) + (1.0L/144.0L)*x_0*pow(x_1, 2) + (7.0L/5760.0L)*x_0*x_1 - 119.0L/57600.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (23.0L/2400.0L)*pow(x_1, 5) - 1.0L/288.0L*pow(x_1, 4) + (5.0L/576.0L)*pow(x_1, 3) - 961.0L/11520.0L*pow(x_1, 2) + (23.0L/38400.0L)*x_1 + 10003.0L/55296.0L;
       return __pp_r144___result;
    }
    static inline O __pp_r145__(const I &x_0, const I &x_1) {
       O __pp_r145___result;
       __pp_r145___result = -11.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 - 1.0L/14400.0L*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) - 41.0L/1440.0L*pow(x_0, 4)*x_1 + (17.0L/768.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/96.0L)*pow(x_0, 3)*x_1 - 77.0L/3456.0L*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 31.0L/576.0L*pow(x_0, 2)*x_1 - 887.0L/15360.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 1.0L/180.0L*x_0*pow(x_1, 4) + (1.0L/24.0L)*x_0*pow(x_1, 3) - 23.0L/288.0L*x_0*pow(x_1, 2) + (71.0L/1280.0L)*x_0*x_1 - 3601.0L/230400.0L*x_0 - 1.0L/4320.0L*pow(x_1, 6) - 71.0L/7200.0L*pow(x_1, 5) + (5.0L/96.0L)*pow(x_1, 4) - 121.0L/1728.0L*pow(x_1, 3) - 3.0L/160.0L*pow(x_1, 2) - 839.0L/28800.0L*x_1 + 34433.0L/184320.0L;
       return __pp_r145___result;
    }
    static inline O __pp_r146__(const I &x_0, const I &x_1) {
       O __pp_r146___result;
       __pp_r146___result = -7.0L/5400.0L*pow(x_0, 6) + (19.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/3600.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/40.0L*pow(x_0, 4)*x_1 + (23.0L/1152.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 3)*x_1 - 13.0L/864.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 2)*x_1 - 1643.0L/23040.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/45.0L)*x_0*pow(x_1, 4) - 1.0L/36.0L*x_0*pow(x_1, 3) + (1.0L/144.0L)*x_0*pow(x_1, 2) + (7.0L/5760.0L)*x_0*x_1 - 119.0L/57600.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (1.0L/800.0L)*pow(x_1, 5) + (5.0L/288.0L)*pow(x_1, 4) - 7.0L/576.0L*pow(x_1, 3) - 841.0L/11520.0L*pow(x_1, 2) - 77.0L/38400.0L*x_1 + 50087.0L/276480.0L;
       return __pp_r146___result;
    }
    static inline O __pp_r147__(const I &x_0, const I &x_1) {
       O __pp_r147___result;
       __pp_r147___result = -11.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 + (1.0L/14400.0L)*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) - 13.0L/480.0L*pow(x_0, 4)*x_1 + (239.0L/11520.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (59.0L/1440.0L)*pow(x_0, 3)*x_1 - 287.0L/17280.0L*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/240.0L)*pow(x_0, 2)*pow(x_1, 2) - 19.0L/960.0L*pow(x_0, 2)*x_1 - 641.0L/9216.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (1.0L/180.0L)*x_0*pow(x_1, 4) - 1.0L/360.0L*x_0*pow(x_1, 3) - 17.0L/1440.0L*x_0*pow(x_1, 2) + (19.0L/2304.0L)*x_0*x_1 - 719.0L/230400.0L*x_0 - 1.0L/4320.0L*pow(x_1, 6) + (7.0L/2400.0L)*pow(x_1, 5) + (13.0L/1440.0L)*pow(x_1, 4) - 11.0L/2880.0L*pow(x_1, 3) - 11.0L/144.0L*pow(x_1, 2) - 29.0L/19200.0L*x_1 + 500879.0L/2764800.0L;
       return __pp_r147___result;
    }
    static inline O __pp_r148__(const I &x_0, const I &x_1) {
       O __pp_r148___result;
       __pp_r148___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 + (7.0L/4800.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 29.0L/1440.0L*pow(x_0, 4)*x_1 + (199.0L/11520.0L)*pow(x_0, 4) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) + (13.0L/480.0L)*pow(x_0, 3)*x_1 - 23.0L/1920.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) - 17.0L/2880.0L*pow(x_0, 2)*x_1 - 673.0L/9216.0L*pow(x_0, 2) + (1.0L/80.0L)*x_0*pow(x_1, 4) - 1.0L/60.0L*x_0*pow(x_1, 3) + (1.0L/480.0L)*x_0*pow(x_1, 2) + (1.0L/768.0L)*x_0*x_1 - 133.0L/76800.0L*x_0 - 1.0L/2160.0L*pow(x_1, 6) + (31.0L/7200.0L)*pow(x_1, 5) + (1.0L/180.0L)*pow(x_1, 4) + (7.0L/8640.0L)*pow(x_1, 3) - 23.0L/288.0L*pow(x_1, 2) - 7.0L/57600.0L*x_1 + 500239.0L/2764800.0L;
       return __pp_r148___result;
    }
    static inline O __pp_r149__(const I &x_0, const I &x_1) {
       O __pp_r149___result;
       __pp_r149___result = -7.0L/4800.0L*pow(x_0, 6) + (17.0L/3600.0L)*pow(x_0, 5)*x_1 + (19.0L/14400.0L)*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) - 31.0L/1440.0L*pow(x_0, 4)*x_1 + (67.0L/3840.0L)*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 - 209.0L/17280.0L*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) - 19.0L/2880.0L*pow(x_0, 2)*x_1 - 1121.0L/15360.0L*pow(x_0, 2) + (2.0L/225.0L)*x_0*pow(x_1, 5) + (1.0L/720.0L)*x_0*pow(x_1, 4) - 1.0L/90.0L*x_0*pow(x_1, 3) + (1.0L/1440.0L)*x_0*pow(x_1, 2) + (17.0L/11520.0L)*x_0*x_1 - 401.0L/230400.0L*x_0 + (1.0L/400.0L)*pow(x_1, 6) - 1.0L/7200.0L*pow(x_1, 5) + (1.0L/120.0L)*pow(x_1, 4) - 1.0L/8640.0L*pow(x_1, 3) - 51.0L/640.0L*pow(x_1, 2) - 1.0L/7200.0L*x_1 + 166747.0L/921600.0L;
       return __pp_r149___result;
    }
    static inline O __pp_r150__(const I &x_0, const I &x_1) {
       O __pp_r150___result;
       __pp_r150___result = -11.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 + (1.0L/14400.0L)*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) - 13.0L/480.0L*pow(x_0, 4)*x_1 + (239.0L/11520.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (59.0L/1440.0L)*pow(x_0, 3)*x_1 - 287.0L/17280.0L*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/240.0L)*pow(x_0, 2)*pow(x_1, 2) - 19.0L/960.0L*pow(x_0, 2)*x_1 - 641.0L/9216.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (1.0L/180.0L)*x_0*pow(x_1, 4) - 1.0L/360.0L*x_0*pow(x_1, 3) - 17.0L/1440.0L*x_0*pow(x_1, 2) + (19.0L/2304.0L)*x_0*x_1 - 719.0L/230400.0L*x_0 - 1.0L/4320.0L*pow(x_1, 6) + (7.0L/2400.0L)*pow(x_1, 5) + (13.0L/1440.0L)*pow(x_1, 4) - 11.0L/2880.0L*pow(x_1, 3) - 11.0L/144.0L*pow(x_1, 2) - 29.0L/19200.0L*x_1 + 500879.0L/2764800.0L;
       return __pp_r150___result;
    }
    static inline O __pp_r151__(const I &x_0, const I &x_1) {
       O __pp_r151___result;
       __pp_r151___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 + (7.0L/4800.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 29.0L/1440.0L*pow(x_0, 4)*x_1 + (199.0L/11520.0L)*pow(x_0, 4) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) + (13.0L/480.0L)*pow(x_0, 3)*x_1 - 23.0L/1920.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) - 17.0L/2880.0L*pow(x_0, 2)*x_1 - 673.0L/9216.0L*pow(x_0, 2) + (1.0L/80.0L)*x_0*pow(x_1, 4) - 1.0L/60.0L*x_0*pow(x_1, 3) + (1.0L/480.0L)*x_0*pow(x_1, 2) + (1.0L/768.0L)*x_0*x_1 - 133.0L/76800.0L*x_0 - 1.0L/2160.0L*pow(x_1, 6) + (31.0L/7200.0L)*pow(x_1, 5) + (1.0L/180.0L)*pow(x_1, 4) + (7.0L/8640.0L)*pow(x_1, 3) - 23.0L/288.0L*pow(x_1, 2) - 7.0L/57600.0L*x_1 + 500239.0L/2764800.0L;
       return __pp_r151___result;
    }
    static inline O __pp_r152__(const I &x_0, const I &x_1) {
       O __pp_r152___result;
       __pp_r152___result = -7.0L/4800.0L*pow(x_0, 6) + (17.0L/3600.0L)*pow(x_0, 5)*x_1 + (19.0L/14400.0L)*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) - 31.0L/1440.0L*pow(x_0, 4)*x_1 + (67.0L/3840.0L)*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 - 209.0L/17280.0L*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) - 19.0L/2880.0L*pow(x_0, 2)*x_1 - 1121.0L/15360.0L*pow(x_0, 2) + (2.0L/225.0L)*x_0*pow(x_1, 5) + (1.0L/720.0L)*x_0*pow(x_1, 4) - 1.0L/90.0L*x_0*pow(x_1, 3) + (1.0L/1440.0L)*x_0*pow(x_1, 2) + (17.0L/11520.0L)*x_0*x_1 - 401.0L/230400.0L*x_0 + (1.0L/400.0L)*pow(x_1, 6) - 1.0L/7200.0L*pow(x_1, 5) + (1.0L/120.0L)*pow(x_1, 4) - 1.0L/8640.0L*pow(x_1, 3) - 51.0L/640.0L*pow(x_1, 2) - 1.0L/7200.0L*x_1 + 166747.0L/921600.0L;
       return __pp_r152___result;
    }
    static inline O __pp_r153__(const I &x_0, const I &x_1) {
       O __pp_r153___result;
       __pp_r153___result = -1.0L/2160.0L*pow(x_0, 6) + (11.0L/2400.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 7.0L/720.0L*pow(x_0, 4)*x_1 + (1.0L/180.0L)*pow(x_0, 4) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 3)*x_1 + (1.0L/960.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/1440.0L*pow(x_0, 2)*x_1 - 23.0L/288.0L*pow(x_0, 2) + (1.0L/80.0L)*x_0*pow(x_1, 4) - 1.0L/60.0L*x_0*pow(x_1, 3) + (1.0L/480.0L)*x_0*pow(x_1, 2) + (1.0L/38400.0L)*x_0 - 1.0L/2160.0L*pow(x_1, 6) + (31.0L/7200.0L)*pow(x_1, 5) + (1.0L/180.0L)*pow(x_1, 4) + (7.0L/8640.0L)*pow(x_1, 3) - 23.0L/288.0L*pow(x_1, 2) + (1.0L/115200.0L)*x_1 + 15617.0L/86400.0L;
       return __pp_r153___result;
    }
    static inline O __pp_r154__(const I &x_0, const I &x_1) {
       O __pp_r154___result;
       __pp_r154___result = -1.0L/2400.0L*pow(x_0, 6) + (1.0L/1800.0L)*pow(x_0, 5)*x_1 + (1.0L/225.0L)*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/90.0L*pow(x_0, 4)*x_1 + (11.0L/1920.0L)*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (13.0L/720.0L)*pow(x_0, 3)*x_1 + (1.0L/1080.0L)*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/720.0L*pow(x_0, 2)*x_1 - 613.0L/7680.0L*pow(x_0, 2) + (2.0L/225.0L)*x_0*pow(x_1, 5) + (1.0L/720.0L)*x_0*pow(x_1, 4) - 1.0L/90.0L*x_0*pow(x_1, 3) + (1.0L/1440.0L)*x_0*pow(x_1, 2) + (1.0L/5760.0L)*x_0*x_1 + (1.0L/57600.0L)*x_0 + (1.0L/400.0L)*pow(x_1, 6) - 1.0L/7200.0L*pow(x_1, 5) + (1.0L/120.0L)*pow(x_1, 4) - 1.0L/8640.0L*pow(x_1, 3) - 51.0L/640.0L*pow(x_1, 2) - 1.0L/115200.0L*x_1 + 83291.0L/460800.0L;
       return __pp_r154___result;
    }
    static inline O __pp_r155__(const I &x_0, const I &x_1) {
       O __pp_r155___result;
       __pp_r155___result = (1.0L/400.0L)*pow(x_0, 6) - 2.0L/225.0L*pow(x_0, 5)*x_1 + (1.0L/7200.0L)*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/720.0L)*pow(x_0, 4)*x_1 + (1.0L/120.0L)*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/90.0L)*pow(x_0, 3)*x_1 + (1.0L/8640.0L)*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/1440.0L)*pow(x_0, 2)*x_1 - 51.0L/640.0L*pow(x_0, 2) - 1.0L/1800.0L*x_0*pow(x_1, 5) + (1.0L/90.0L)*x_0*pow(x_1, 4) - 13.0L/720.0L*x_0*pow(x_1, 3) + (1.0L/720.0L)*x_0*pow(x_1, 2) - 1.0L/5760.0L*x_0*x_1 + (1.0L/115200.0L)*x_0 - 1.0L/2400.0L*pow(x_1, 6) + (1.0L/225.0L)*pow(x_1, 5) + (11.0L/1920.0L)*pow(x_1, 4) + (1.0L/1080.0L)*pow(x_1, 3) - 613.0L/7680.0L*pow(x_1, 2) + (1.0L/57600.0L)*x_1 + 83291.0L/460800.0L;
       return __pp_r155___result;
    }
    static inline O __pp_r156__(const I &x_0, const I &x_1) {
       O __pp_r156___result;
       __pp_r156___result = (11.0L/4320.0L)*pow(x_0, 6) - 1.0L/120.0L*pow(x_0, 5)*x_1 + (1.0L/144.0L)*pow(x_0, 4)*pow(x_1, 2) + (49.0L/5760.0L)*pow(x_0, 4) + (1.0L/80.0L)*pow(x_0, 3)*x_1 + (1.0L/144.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) - 367.0L/4608.0L*pow(x_0, 2) + (1.0L/120.0L)*x_0*pow(x_1, 5) - 1.0L/80.0L*x_0*pow(x_1, 3) + (11.0L/4320.0L)*pow(x_1, 6) + (49.0L/5760.0L)*pow(x_1, 4) - 367.0L/4608.0L*pow(x_1, 2) + 124937.0L/691200.0L;
       return __pp_r156___result;
    }
    static inline O __pp_r157__(const I &x_0, const I &x_1) {
       O __pp_r157___result;
       __pp_r157___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/960.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 - 25.0L/128.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/12.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/64.0L*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/3600.0L)*x_0*pow(x_1, 5) - 1.0L/32.0L*x_0*pow(x_1, 4) + (35.0L/96.0L)*x_0*pow(x_1, 3) - 325.0L/192.0L*x_0*pow(x_1, 2) + (1375.0L/384.0L)*x_0*x_1 - 8875.0L/3072.0L*x_0 + (19.0L/43200.0L)*pow(x_1, 6) - 7.0L/960.0L*pow(x_1, 5) + (65.0L/768.0L)*pow(x_1, 4) - 75.0L/128.0L*pow(x_1, 3) + (6625.0L/3072.0L)*pow(x_1, 2) - 12125.0L/3072.0L*x_1 + 4375.0L/1536.0L;
       return __pp_r157___result;
    }
    static inline O __pp_r158__(const I &x_0, const I &x_1) {
       O __pp_r158___result;
       __pp_r158___result = -1.0L/7200.0L*pow(x_0, 5) - 1.0L/720.0L*pow(x_0, 4)*x_1 + (7.0L/1440.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (7.0L/180.0L)*pow(x_0, 3)*x_1 - 589.0L/8640.0L*pow(x_0, 3) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) - 589.0L/1440.0L*pow(x_0, 2)*x_1 + (1379.0L/2880.0L)*pow(x_0, 2) - 1.0L/240.0L*x_0*pow(x_1, 5) + (59.0L/1440.0L)*x_0*pow(x_1, 4) - 151.0L/1440.0L*x_0*pow(x_1, 3) - 481.0L/2880.0L*x_0*pow(x_1, 2) + (12689.0L/11520.0L)*x_0*x_1 - 73583.0L/57600.0L*x_0 - 1.0L/960.0L*pow(x_1, 6) + (311.0L/14400.0L)*pow(x_1, 5) - 1729.0L/11520.0L*pow(x_1, 4) + (7451.0L/17280.0L)*pow(x_1, 3) - 14869.0L/46080.0L*pow(x_1, 2) - 166789.0L/230400.0L*x_1 + 3048191.0L/2764800.0L;
       return __pp_r158___result;
    }
    static inline O __pp_r159__(const I &x_0, const I &x_1) {
       O __pp_r159___result;
       __pp_r159___result = -1.0L/2160.0L*pow(x_0, 6) - 1.0L/360.0L*pow(x_0, 5)*x_1 + (79.0L/7200.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (13.0L/240.0L)*pow(x_0, 4)*x_1 - 17.0L/160.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) + (19.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) - 73.0L/180.0L*pow(x_0, 3)*x_1 + (4531.0L/8640.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/10.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (219.0L/160.0L)*pow(x_0, 2)*x_1 - 1247.0L/960.0L*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) + (139.0L/1440.0L)*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) + (4639.0L/2880.0L)*x_0*pow(x_1, 2) - 28271.0L/11520.0L*x_0*x_1 + (90257.0L/57600.0L)*x_0 - 13.0L/8640.0L*pow(x_1, 6) + (157.0L/4800.0L)*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) + (5897.0L/5760.0L)*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) + (162857.0L/76800.0L)*x_1 - 731563.0L/921600.0L;
       return __pp_r159___result;
    }
    static inline O __pp_r160__(const I &x_0, const I &x_1) {
       O __pp_r160___result;
       __pp_r160___result = -1.0L/2160.0L*pow(x_0, 6) - 1.0L/360.0L*pow(x_0, 5)*x_1 + (79.0L/7200.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (13.0L/240.0L)*pow(x_0, 4)*x_1 - 17.0L/160.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) + (19.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) - 73.0L/180.0L*pow(x_0, 3)*x_1 + (4531.0L/8640.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/10.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (219.0L/160.0L)*pow(x_0, 2)*x_1 - 1247.0L/960.0L*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) + (139.0L/1440.0L)*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) + (4639.0L/2880.0L)*x_0*pow(x_1, 2) - 28271.0L/11520.0L*x_0*x_1 + (90257.0L/57600.0L)*x_0 - 13.0L/8640.0L*pow(x_1, 6) + (157.0L/4800.0L)*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) + (5897.0L/5760.0L)*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) + (162857.0L/76800.0L)*x_1 - 731563.0L/921600.0L;
       return __pp_r160___result;
    }
    static inline O __pp_r161__(const I &x_0, const I &x_1) {
       O __pp_r161___result;
       __pp_r161___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/960.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 - 25.0L/128.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/12.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/64.0L*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/3600.0L)*x_0*pow(x_1, 5) - 1.0L/32.0L*x_0*pow(x_1, 4) + (35.0L/96.0L)*x_0*pow(x_1, 3) - 325.0L/192.0L*x_0*pow(x_1, 2) + (1375.0L/384.0L)*x_0*x_1 - 8875.0L/3072.0L*x_0 + (19.0L/43200.0L)*pow(x_1, 6) - 7.0L/960.0L*pow(x_1, 5) + (65.0L/768.0L)*pow(x_1, 4) - 75.0L/128.0L*pow(x_1, 3) + (6625.0L/3072.0L)*pow(x_1, 2) - 12125.0L/3072.0L*x_1 + 4375.0L/1536.0L;
       return __pp_r161___result;
    }
    static inline O __pp_r162__(const I &x_0, const I &x_1) {
       O __pp_r162___result;
       __pp_r162___result = -1.0L/7200.0L*pow(x_0, 5) - 1.0L/720.0L*pow(x_0, 4)*x_1 + (7.0L/1440.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (7.0L/180.0L)*pow(x_0, 3)*x_1 - 589.0L/8640.0L*pow(x_0, 3) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) - 589.0L/1440.0L*pow(x_0, 2)*x_1 + (1379.0L/2880.0L)*pow(x_0, 2) - 1.0L/240.0L*x_0*pow(x_1, 5) + (59.0L/1440.0L)*x_0*pow(x_1, 4) - 151.0L/1440.0L*x_0*pow(x_1, 3) - 481.0L/2880.0L*x_0*pow(x_1, 2) + (12689.0L/11520.0L)*x_0*x_1 - 73583.0L/57600.0L*x_0 - 1.0L/960.0L*pow(x_1, 6) + (311.0L/14400.0L)*pow(x_1, 5) - 1729.0L/11520.0L*pow(x_1, 4) + (7451.0L/17280.0L)*pow(x_1, 3) - 14869.0L/46080.0L*pow(x_1, 2) - 166789.0L/230400.0L*x_1 + 3048191.0L/2764800.0L;
       return __pp_r162___result;
    }
    static inline O __pp_r163__(const I &x_0, const I &x_1) {
       O __pp_r163___result;
       __pp_r163___result = -1.0L/2160.0L*pow(x_0, 6) - 1.0L/360.0L*pow(x_0, 5)*x_1 + (79.0L/7200.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (13.0L/240.0L)*pow(x_0, 4)*x_1 - 17.0L/160.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) + (19.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) - 73.0L/180.0L*pow(x_0, 3)*x_1 + (4531.0L/8640.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/10.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (219.0L/160.0L)*pow(x_0, 2)*x_1 - 1247.0L/960.0L*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) + (139.0L/1440.0L)*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) + (4639.0L/2880.0L)*x_0*pow(x_1, 2) - 28271.0L/11520.0L*x_0*x_1 + (90257.0L/57600.0L)*x_0 - 13.0L/8640.0L*pow(x_1, 6) + (157.0L/4800.0L)*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) + (5897.0L/5760.0L)*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) + (162857.0L/76800.0L)*x_1 - 731563.0L/921600.0L;
       return __pp_r163___result;
    }
    static inline O __pp_r164__(const I &x_0, const I &x_1) {
       O __pp_r164___result;
       __pp_r164___result = -1.0L/2160.0L*pow(x_0, 6) - 1.0L/360.0L*pow(x_0, 5)*x_1 + (79.0L/7200.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (13.0L/240.0L)*pow(x_0, 4)*x_1 - 17.0L/160.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) + (19.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) - 73.0L/180.0L*pow(x_0, 3)*x_1 + (4531.0L/8640.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/10.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (219.0L/160.0L)*pow(x_0, 2)*x_1 - 1247.0L/960.0L*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) + (139.0L/1440.0L)*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) + (4639.0L/2880.0L)*x_0*pow(x_1, 2) - 28271.0L/11520.0L*x_0*x_1 + (90257.0L/57600.0L)*x_0 - 13.0L/8640.0L*pow(x_1, 6) + (157.0L/4800.0L)*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) + (5897.0L/5760.0L)*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) + (162857.0L/76800.0L)*x_1 - 731563.0L/921600.0L;
       return __pp_r164___result;
    }
    static inline O __pp_r165__(const I &x_0, const I &x_1) {
       O __pp_r165___result;
       __pp_r165___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 + (143.0L/14400.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 4)*x_1 - 81.0L/1280.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) + (19.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) - 449.0L/1440.0L*pow(x_0, 3)*x_1 + (6227.0L/17280.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/10.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (393.0L/320.0L)*pow(x_0, 2)*x_1 - 15767.0L/15360.0L*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) + (139.0L/1440.0L)*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) + (4639.0L/2880.0L)*x_0*pow(x_1, 2) - 1691.0L/720.0L*x_0*x_1 + (311213.0L/230400.0L)*x_0 - 13.0L/8640.0L*pow(x_1, 6) + (157.0L/4800.0L)*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) + (5897.0L/5760.0L)*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) + (160427.0L/76800.0L)*x_1 - 334799.0L/460800.0L;
       return __pp_r165___result;
    }
    static inline O __pp_r166__(const I &x_0, const I &x_1) {
       O __pp_r166___result;
       __pp_r166___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 + (143.0L/14400.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 4)*x_1 - 81.0L/1280.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) + (19.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) - 449.0L/1440.0L*pow(x_0, 3)*x_1 + (6227.0L/17280.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/10.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (393.0L/320.0L)*pow(x_0, 2)*x_1 - 15767.0L/15360.0L*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) + (139.0L/1440.0L)*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) + (4639.0L/2880.0L)*x_0*pow(x_1, 2) - 1691.0L/720.0L*x_0*x_1 + (311213.0L/230400.0L)*x_0 - 13.0L/8640.0L*pow(x_1, 6) + (157.0L/4800.0L)*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) + (5897.0L/5760.0L)*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) + (160427.0L/76800.0L)*x_1 - 334799.0L/460800.0L;
       return __pp_r166___result;
    }
    static inline O __pp_r167__(const I &x_0, const I &x_1) {
       O __pp_r167___result;
       __pp_r167___result = -1.0L/43200.0L*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 + (37.0L/4800.0L)*pow(x_0, 5) - 1.0L/720.0L*pow(x_0, 4)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 4)*x_1 - 713.0L/11520.0L*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/10.0L)*pow(x_0, 3)*pow(x_1, 2) - 151.0L/480.0L*pow(x_0, 3)*x_1 + (691.0L/1920.0L)*pow(x_0, 3) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 4) + (37.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 263.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (3539.0L/2880.0L)*pow(x_0, 2)*x_1 - 47297.0L/46080.0L*pow(x_0, 2) - 13.0L/1800.0L*x_0*pow(x_1, 5) + (23.0L/240.0L)*x_0*pow(x_1, 4) - 11.0L/20.0L*x_0*pow(x_1, 3) + (773.0L/480.0L)*x_0*pow(x_1, 2) - 9019.0L/3840.0L*x_0*x_1 + (34579.0L/25600.0L)*x_0 - 1.0L/675.0L*pow(x_1, 6) + (59.0L/1800.0L)*pow(x_1, 5) - 47.0L/180.0L*pow(x_1, 4) + (4423.0L/4320.0L)*pow(x_1, 3) - 24197.0L/11520.0L*pow(x_1, 2) + (240641.0L/115200.0L)*x_1 - 2008793.0L/2764800.0L;
       return __pp_r167___result;
    }
    static inline O __pp_r168__(const I &x_0, const I &x_1) {
       O __pp_r168___result;
       __pp_r168___result = -1.0L/43200.0L*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 + (37.0L/4800.0L)*pow(x_0, 5) - 1.0L/720.0L*pow(x_0, 4)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 4)*x_1 - 713.0L/11520.0L*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/10.0L)*pow(x_0, 3)*pow(x_1, 2) - 151.0L/480.0L*pow(x_0, 3)*x_1 + (691.0L/1920.0L)*pow(x_0, 3) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 4) + (37.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 263.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (3539.0L/2880.0L)*pow(x_0, 2)*x_1 - 47297.0L/46080.0L*pow(x_0, 2) - 13.0L/1800.0L*x_0*pow(x_1, 5) + (23.0L/240.0L)*x_0*pow(x_1, 4) - 11.0L/20.0L*x_0*pow(x_1, 3) + (773.0L/480.0L)*x_0*pow(x_1, 2) - 9019.0L/3840.0L*x_0*x_1 + (34579.0L/25600.0L)*x_0 - 1.0L/675.0L*pow(x_1, 6) + (59.0L/1800.0L)*pow(x_1, 5) - 47.0L/180.0L*pow(x_1, 4) + (4423.0L/4320.0L)*pow(x_1, 3) - 24197.0L/11520.0L*pow(x_1, 2) + (240641.0L/115200.0L)*x_1 - 2008793.0L/2764800.0L;
       return __pp_r168___result;
    }
    static inline O __pp_r169__(const I &x_0, const I &x_1) {
       O __pp_r169___result;
       __pp_r169___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (91.0L/14400.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (17.0L/480.0L)*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (31.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 + (6139.0L/17280.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (1193.0L/960.0L)*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) + (4.0L/45.0L)*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) + (2299.0L/1440.0L)*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 + (310891.0L/230400.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (41.0L/1200.0L)*pow(x_1, 5) - 127.0L/480.0L*pow(x_1, 4) + (1481.0L/1440.0L)*pow(x_1, 3) - 2693.0L/1280.0L*pow(x_1, 2) + (80267.0L/38400.0L)*x_1 - 669811.0L/921600.0L;
       return __pp_r169___result;
    }
    static inline O __pp_r170__(const I &x_0, const I &x_1) {
       O __pp_r170___result;
       __pp_r170___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (91.0L/14400.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (17.0L/480.0L)*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (31.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 + (6139.0L/17280.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (1193.0L/960.0L)*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) + (4.0L/45.0L)*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) + (2299.0L/1440.0L)*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 + (310891.0L/230400.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (41.0L/1200.0L)*pow(x_1, 5) - 127.0L/480.0L*pow(x_1, 4) + (1481.0L/1440.0L)*pow(x_1, 3) - 2693.0L/1280.0L*pow(x_1, 2) + (80267.0L/38400.0L)*x_1 - 669811.0L/921600.0L;
       return __pp_r170___result;
    }
    static inline O __pp_r171__(const I &x_0, const I &x_1) {
       O __pp_r171___result;
       __pp_r171___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (91.0L/14400.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (17.0L/480.0L)*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (31.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 + (6139.0L/17280.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (1193.0L/960.0L)*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) + (4.0L/45.0L)*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) + (2299.0L/1440.0L)*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 + (310891.0L/230400.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (41.0L/1200.0L)*pow(x_1, 5) - 127.0L/480.0L*pow(x_1, 4) + (1481.0L/1440.0L)*pow(x_1, 3) - 2693.0L/1280.0L*pow(x_1, 2) + (80267.0L/38400.0L)*x_1 - 669811.0L/921600.0L;
       return __pp_r171___result;
    }
    static inline O __pp_r172__(const I &x_0, const I &x_1) {
       O __pp_r172___result;
       __pp_r172___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (91.0L/14400.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (17.0L/480.0L)*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (31.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 + (6139.0L/17280.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (1193.0L/960.0L)*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) + (4.0L/45.0L)*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) + (2299.0L/1440.0L)*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 + (310891.0L/230400.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (41.0L/1200.0L)*pow(x_1, 5) - 127.0L/480.0L*pow(x_1, 4) + (1481.0L/1440.0L)*pow(x_1, 3) - 2693.0L/1280.0L*pow(x_1, 2) + (80267.0L/38400.0L)*x_1 - 669811.0L/921600.0L;
       return __pp_r172___result;
    }
    static inline O __pp_r173__(const I &x_0, const I &x_1) {
       O __pp_r173___result;
       __pp_r173___result = -11.0L/14400.0L*pow(x_0, 6) + (11.0L/1800.0L)*pow(x_0, 5)*x_1 - 11.0L/2880.0L*pow(x_0, 5) + (1.0L/160.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/96.0L*pow(x_0, 4)*x_1 + (57.0L/1280.0L)*pow(x_0, 4) + (13.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (253.0L/1440.0L)*pow(x_0, 3)*x_1 - 1711.0L/17280.0L*pow(x_0, 3) + (7.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/8.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) - 129.0L/320.0L*pow(x_0, 2)*x_1 + (1567.0L/15360.0L)*pow(x_0, 2) + (37.0L/3600.0L)*x_0*pow(x_1, 5) - 7.0L/72.0L*x_0*pow(x_1, 4) + (133.0L/360.0L)*x_0*pow(x_1, 3) - 961.0L/1440.0L*x_0*pow(x_1, 2) + (6721.0L/11520.0L)*x_0*x_1 - 9491.0L/46080.0L*x_0 + (23.0L/7200.0L)*pow(x_1, 6) - 7.0L/160.0L*pow(x_1, 5) + (97.0L/480.0L)*pow(x_1, 4) - 1223.0L/2880.0L*pow(x_1, 3) + (871.0L/1920.0L)*pow(x_1, 2) - 1403.0L/3840.0L*x_1 + 264251.0L/921600.0L;
       return __pp_r173___result;
    }
    static inline O __pp_r174__(const I &x_0, const I &x_1) {
       O __pp_r174___result;
       __pp_r174___result = -7.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 - 41.0L/14400.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 61.0L/1440.0L*pow(x_0, 4)*x_1 + (83.0L/2304.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 11.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (31.0L/288.0L)*pow(x_0, 3)*x_1 - 205.0L/3456.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 17.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 95.0L/576.0L*pow(x_0, 2)*x_1 - 101.0L/46080.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 7.0L/360.0L*x_0*pow(x_1, 4) + (7.0L/72.0L)*x_0*pow(x_1, 3) - 55.0L/288.0L*x_0*pow(x_1, 2) + (1919.0L/11520.0L)*x_0*x_1 - 13841.0L/230400.0L*x_0 + (1.0L/4320.0L)*pow(x_1, 6) - 91.0L/7200.0L*pow(x_1, 5) + (19.0L/288.0L)*pow(x_1, 4) - 185.0L/1728.0L*pow(x_1, 3) + (53.0L/1440.0L)*pow(x_1, 2) - 2119.0L/28800.0L*x_1 + 111491.0L/552960.0L;
       return __pp_r174___result;
    }
    static inline O __pp_r175__(const I &x_0, const I &x_1) {
       O __pp_r175___result;
       __pp_r175___result = (31.0L/43200.0L)*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 - 29.0L/4800.0L*pow(x_0, 5) + (17.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 67.0L/1440.0L*pow(x_0, 4)*x_1 + (529.0L/11520.0L)*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 3.0L/40.0L*pow(x_0, 3)*pow(x_1, 2) + (83.0L/480.0L)*pow(x_0, 3)*x_1 - 191.0L/1920.0L*pow(x_0, 3) + (23.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 11.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (169.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 1159.0L/2880.0L*pow(x_0, 2)*x_1 + (941.0L/9216.0L)*pow(x_0, 2) + (1.0L/100.0L)*x_0*pow(x_1, 5) - 47.0L/480.0L*x_0*pow(x_1, 4) + (59.0L/160.0L)*x_0*pow(x_1, 3) - 641.0L/960.0L*x_0*pow(x_1, 2) + (7.0L/12.0L)*x_0*x_1 - 5273.0L/25600.0L*x_0 + (139.0L/43200.0L)*pow(x_1, 6) - 629.0L/14400.0L*pow(x_1, 5) + (2329.0L/11520.0L)*pow(x_1, 4) - 7337.0L/17280.0L*pow(x_1, 3) + (4181.0L/9216.0L)*pow(x_1, 2) - 84179.0L/230400.0L*x_1 + 396377.0L/1382400.0L;
       return __pp_r175___result;
    }
    static inline O __pp_r176__(const I &x_0, const I &x_1) {
       O __pp_r176___result;
       __pp_r176___result = (29.0L/43200.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 - 73.0L/14400.0L*pow(x_0, 5) + (13.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 53.0L/1440.0L*pow(x_0, 4)*x_1 + (431.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 13.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (151.0L/1440.0L)*pow(x_0, 3)*x_1 - 1033.0L/17280.0L*pow(x_0, 3) + (7.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 2.0L/45.0L*pow(x_0, 2)*pow(x_1, 3) + (71.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 473.0L/2880.0L*pow(x_0, 2)*x_1 - 97.0L/46080.0L*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 29.0L/1440.0L*x_0*pow(x_1, 4) + (139.0L/1440.0L)*x_0*pow(x_1, 3) - 551.0L/2880.0L*x_0*pow(x_1, 2) + (959.0L/5760.0L)*x_0*x_1 - 13843.0L/230400.0L*x_0 + (11.0L/43200.0L)*pow(x_1, 6) - 181.0L/14400.0L*pow(x_1, 5) + (761.0L/11520.0L)*pow(x_1, 4) - 1849.0L/17280.0L*pow(x_1, 3) + (1697.0L/46080.0L)*pow(x_1, 2) - 16951.0L/230400.0L*x_1 + 34841.0L/172800.0L;
       return __pp_r176___result;
    }
    static inline O __pp_r177__(const I &x_0, const I &x_1) {
       O __pp_r177___result;
       __pp_r177___result = -7.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 - 37.0L/14400.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 19.0L/480.0L*pow(x_0, 4)*x_1 + (39.0L/1280.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (91.0L/1440.0L)*pow(x_0, 3)*x_1 - 253.0L/17280.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/40.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) + (33.0L/320.0L)*pow(x_0, 2)*x_1 - 2807.0L/15360.0L*pow(x_0, 2) - 1.0L/360.0L*x_0*pow(x_1, 5) + (49.0L/1440.0L)*x_0*pow(x_1, 4) - 251.0L/1440.0L*x_0*pow(x_1, 3) + (1399.0L/2880.0L)*x_0*pow(x_1, 2) - 119.0L/180.0L*x_0*x_1 + (77933.0L/230400.0L)*x_0 - 7.0L/8640.0L*pow(x_1, 6) + (97.0L/4800.0L)*pow(x_1, 5) - 643.0L/3840.0L*pow(x_1, 4) + (3737.0L/5760.0L)*pow(x_1, 3) - 19303.0L/15360.0L*pow(x_1, 2) + (82667.0L/76800.0L)*x_1 - 101519.0L/460800.0L;
       return __pp_r177___result;
    }
    static inline O __pp_r178__(const I &x_0, const I &x_1) {
       O __pp_r178___result;
       __pp_r178___result = (29.0L/43200.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 - 23.0L/4800.0L*pow(x_0, 5) + (13.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 49.0L/1440.0L*pow(x_0, 4)*x_1 + (367.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/40.0L*pow(x_0, 3)*pow(x_1, 2) + (29.0L/480.0L)*pow(x_0, 3)*x_1 - 29.0L/1920.0L*pow(x_0, 3) + (7.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/45.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (299.0L/2880.0L)*pow(x_0, 2)*x_1 - 8417.0L/46080.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/30.0L)*x_0*pow(x_1, 4) - 7.0L/40.0L*x_0*pow(x_1, 3) + (233.0L/480.0L)*x_0*pow(x_1, 2) - 2539.0L/3840.0L*x_0*x_1 + (8659.0L/25600.0L)*x_0 - 17.0L/21600.0L*pow(x_1, 6) + (73.0L/3600.0L)*pow(x_1, 5) - 241.0L/1440.0L*pow(x_1, 4) + (2803.0L/4320.0L)*pow(x_1, 3) - 14477.0L/11520.0L*pow(x_1, 2) + (124001.0L/115200.0L)*x_1 - 609113.0L/2764800.0L;
       return __pp_r178___result;
    }
    static inline O __pp_r179__(const I &x_0, const I &x_1) {
       O __pp_r179___result;
       __pp_r179___result = (19.0L/43200.0L)*pow(x_0, 6) + (1.0L/400.0L)*pow(x_0, 5)*x_1 - 89.0L/14400.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 13.0L/480.0L*pow(x_0, 4)*x_1 + (109.0L/3840.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (107.0L/1440.0L)*pow(x_0, 3)*x_1 - 341.0L/17280.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) - 1.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (113.0L/960.0L)*pow(x_0, 2)*x_1 - 953.0L/5120.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) + (19.0L/720.0L)*x_0*pow(x_1, 4) - 29.0L/180.0L*x_0*pow(x_1, 3) + (679.0L/1440.0L)*x_0*pow(x_1, 2) - 7537.0L/11520.0L*x_0*x_1 + (77611.0L/230400.0L)*x_0 - 11.0L/10800.0L*pow(x_1, 6) + (13.0L/600.0L)*pow(x_1, 5) - 41.0L/240.0L*pow(x_1, 4) + (941.0L/1440.0L)*pow(x_1, 3) - 1613.0L/1280.0L*pow(x_1, 2) + (41387.0L/38400.0L)*x_1 - 203251.0L/921600.0L;
       return __pp_r179___result;
    }
    static inline O __pp_r180__(const I &x_0, const I &x_1) {
       O __pp_r180___result;
       __pp_r180___result = -11.0L/14400.0L*pow(x_0, 6) + (11.0L/1800.0L)*pow(x_0, 5)*x_1 - 11.0L/2880.0L*pow(x_0, 5) + (1.0L/160.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/96.0L*pow(x_0, 4)*x_1 + (57.0L/1280.0L)*pow(x_0, 4) + (13.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (253.0L/1440.0L)*pow(x_0, 3)*x_1 - 1711.0L/17280.0L*pow(x_0, 3) + (7.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/8.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) - 129.0L/320.0L*pow(x_0, 2)*x_1 + (1567.0L/15360.0L)*pow(x_0, 2) + (11.0L/1800.0L)*x_0*pow(x_1, 5) - 19.0L/288.0L*x_0*pow(x_1, 4) + (397.0L/1440.0L)*x_0*pow(x_1, 3) - 1517.0L/2880.0L*x_0*pow(x_1, 2) + (2753.0L/5760.0L)*x_0*x_1 - 8033.0L/46080.0L*x_0 + (31.0L/14400.0L)*pow(x_1, 6) - 19.0L/960.0L*pow(x_1, 5) + (221.0L/3840.0L)*pow(x_1, 4) - 151.0L/5760.0L*pow(x_1, 3) - 1807.0L/15360.0L*pow(x_1, 2) + (787.0L/15360.0L)*x_1 + 18907.0L/115200.0L;
       return __pp_r180___result;
    }
    static inline O __pp_r181__(const I &x_0, const I &x_1) {
       O __pp_r181___result;
       __pp_r181___result = (31.0L/43200.0L)*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 - 29.0L/4800.0L*pow(x_0, 5) + (17.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 67.0L/1440.0L*pow(x_0, 4)*x_1 + (529.0L/11520.0L)*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 3.0L/40.0L*pow(x_0, 3)*pow(x_1, 2) + (83.0L/480.0L)*pow(x_0, 3)*x_1 - 191.0L/1920.0L*pow(x_0, 3) + (23.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 11.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (169.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 1159.0L/2880.0L*pow(x_0, 2)*x_1 + (941.0L/9216.0L)*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) - 1.0L/15.0L*x_0*pow(x_1, 4) + (11.0L/40.0L)*x_0*pow(x_1, 3) - 253.0L/480.0L*x_0*pow(x_1, 2) + (367.0L/768.0L)*x_0*x_1 - 4463.0L/25600.0L*x_0 + (47.0L/21600.0L)*pow(x_1, 6) - 71.0L/3600.0L*pow(x_1, 5) + (83.0L/1440.0L)*pow(x_1, 4) - 113.0L/4320.0L*pow(x_1, 3) - 271.0L/2304.0L*pow(x_1, 2) + (5903.0L/115200.0L)*x_1 + 453769.0L/2764800.0L;
       return __pp_r181___result;
    }
    static inline O __pp_r182__(const I &x_0, const I &x_1) {
       O __pp_r182___result;
       __pp_r182___result = (7.0L/14400.0L)*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 107.0L/14400.0L*pow(x_0, 5) + (1.0L/120.0L)*pow(x_0, 4)*pow(x_1, 2) - 19.0L/480.0L*pow(x_0, 4)*x_1 + (163.0L/3840.0L)*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) - 4.0L/45.0L*pow(x_0, 3)*pow(x_1, 2) + (269.0L/1440.0L)*pow(x_0, 3)*x_1 - 1799.0L/17280.0L*pow(x_0, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 4) - 13.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (53.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 373.0L/960.0L*pow(x_0, 2)*x_1 + (101.0L/1024.0L)*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) - 53.0L/720.0L*x_0*pow(x_1, 4) + (13.0L/45.0L)*x_0*pow(x_1, 3) - 779.0L/1440.0L*x_0*pow(x_1, 2) + (1117.0L/2304.0L)*x_0*x_1 - 40487.0L/230400.0L*x_0 + (7.0L/3600.0L)*pow(x_1, 6) - 11.0L/600.0L*pow(x_1, 5) + (13.0L/240.0L)*pow(x_1, 4) - 31.0L/1440.0L*pow(x_1, 3) - 31.0L/256.0L*pow(x_1, 2) + (2021.0L/38400.0L)*x_1 + 151043.0L/921600.0L;
       return __pp_r182___result;
    }
    static inline O __pp_r183__(const I &x_0, const I &x_1) {
       O __pp_r183___result;
       __pp_r183___result = (19.0L/43200.0L)*pow(x_0, 6) + (1.0L/400.0L)*pow(x_0, 5)*x_1 - 89.0L/14400.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 13.0L/480.0L*pow(x_0, 4)*x_1 + (109.0L/3840.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (107.0L/1440.0L)*pow(x_0, 3)*x_1 - 341.0L/17280.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) - 1.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (113.0L/960.0L)*pow(x_0, 2)*x_1 - 953.0L/5120.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) + (19.0L/720.0L)*x_0*pow(x_1, 4) - 29.0L/180.0L*x_0*pow(x_1, 3) + (679.0L/1440.0L)*x_0*pow(x_1, 2) - 7537.0L/11520.0L*x_0*x_1 + (77611.0L/230400.0L)*x_0 - 11.0L/10800.0L*pow(x_1, 6) + (13.0L/600.0L)*pow(x_1, 5) - 41.0L/240.0L*pow(x_1, 4) + (941.0L/1440.0L)*pow(x_1, 3) - 1613.0L/1280.0L*pow(x_1, 2) + (41387.0L/38400.0L)*x_1 - 203251.0L/921600.0L;
       return __pp_r183___result;
    }
    static inline O __pp_r184__(const I &x_0, const I &x_1) {
       O __pp_r184___result;
       __pp_r184___result = (19.0L/43200.0L)*pow(x_0, 6) + (1.0L/400.0L)*pow(x_0, 5)*x_1 - 149.0L/14400.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 13.0L/480.0L*pow(x_0, 4)*x_1 + (149.0L/3840.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (107.0L/1440.0L)*pow(x_0, 3)*x_1 - 521.0L/17280.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) - 1.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (113.0L/960.0L)*pow(x_0, 2)*x_1 - 2779.0L/15360.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) + (19.0L/720.0L)*x_0*pow(x_1, 4) - 29.0L/180.0L*x_0*pow(x_1, 3) + (679.0L/1440.0L)*x_0*pow(x_1, 2) - 7537.0L/11520.0L*x_0*x_1 + (77311.0L/230400.0L)*x_0 - 11.0L/10800.0L*pow(x_1, 6) + (13.0L/600.0L)*pow(x_1, 5) - 41.0L/240.0L*pow(x_1, 4) + (941.0L/1440.0L)*pow(x_1, 3) - 1613.0L/1280.0L*pow(x_1, 2) + (41387.0L/38400.0L)*x_1 - 203131.0L/921600.0L;
       return __pp_r184___result;
    }
    static inline O __pp_r185__(const I &x_0, const I &x_1) {
       O __pp_r185___result;
       __pp_r185___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 - 49.0L/2880.0L*pow(x_0, 5) - 1.0L/96.0L*pow(x_0, 4)*x_1 + (101.0L/3840.0L)*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/18.0L*pow(x_0, 3)*pow(x_1, 2) + (143.0L/1440.0L)*pow(x_0, 3)*x_1 - 737.0L/17280.0L*pow(x_0, 3) - 1.0L/40.0L*pow(x_0, 2)*pow(x_1, 2) + (131.0L/960.0L)*pow(x_0, 2)*x_1 - 2887.0L/15360.0L*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) + (7.0L/288.0L)*x_0*pow(x_1, 4) - 223.0L/1440.0L*x_0*pow(x_1, 3) + (1331.0L/2880.0L)*x_0*pow(x_1, 2) - 233.0L/360.0L*x_0*x_1 + (3073.0L/9216.0L)*x_0 - 1.0L/960.0L*pow(x_1, 6) + (7.0L/320.0L)*pow(x_1, 5) - 659.0L/3840.0L*pow(x_1, 4) + (3773.0L/5760.0L)*pow(x_1, 3) - 6461.0L/5120.0L*pow(x_1, 2) + (16571.0L/15360.0L)*x_1 - 101687.0L/460800.0L;
       return __pp_r185___result;
    }
    static inline O __pp_r186__(const I &x_0, const I &x_1) {
       O __pp_r186___result;
       __pp_r186___result = -7.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 - 41.0L/2880.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 7.0L/288.0L*pow(x_0, 4)*x_1 + (463.0L/11520.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/36.0L*pow(x_0, 3)*pow(x_1, 2) + (7.0L/160.0L)*pow(x_0, 3)*x_1 - 97.0L/17280.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/36.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) + (73.0L/2880.0L)*pow(x_0, 2)*x_1 - 6101.0L/46080.0L*pow(x_0, 2) - 1.0L/360.0L*x_0*pow(x_1, 5) + (11.0L/288.0L)*x_0*pow(x_1, 4) - 101.0L/480.0L*x_0*pow(x_1, 3) + (1651.0L/2880.0L)*x_0*pow(x_1, 2) - 91.0L/120.0L*x_0*x_1 + (17413.0L/46080.0L)*x_0 - 7.0L/8640.0L*pow(x_1, 6) + (11.0L/576.0L)*pow(x_1, 5) - 1817.0L/11520.0L*pow(x_1, 4) + (10679.0L/17280.0L)*pow(x_1, 3) - 55589.0L/46080.0L*pow(x_1, 2) + (9533.0L/9216.0L)*x_1 - 284581.0L/1382400.0L;
       return __pp_r186___result;
    }
    static inline O __pp_r187__(const I &x_0, const I &x_1) {
       O __pp_r187___result;
       __pp_r187___result = (29.0L/43200.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 - 1.0L/320.0L*pow(x_0, 5) + (13.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/96.0L*pow(x_0, 4)*x_1 + (863.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 37.0L/1440.0L*pow(x_0, 3)*x_1 + (301.0L/5760.0L)*pow(x_0, 3) + (7.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (53.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 59.0L/960.0L*pow(x_0, 2)*x_1 - 3601.0L/46080.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/24.0L)*x_0*pow(x_1, 4) - 41.0L/180.0L*x_0*pow(x_1, 3) + (37.0L/60.0L)*x_0*pow(x_1, 2) - 9361.0L/11520.0L*x_0*x_1 + (6221.0L/15360.0L)*x_0 - 17.0L/21600.0L*pow(x_1, 6) + (3.0L/160.0L)*pow(x_1, 5) - 7.0L/45.0L*pow(x_1, 4) + (1759.0L/2880.0L)*pow(x_1, 3) - 13741.0L/11520.0L*pow(x_1, 2) + (49.0L/48.0L)*x_1 - 553537.0L/2764800.0L;
       return __pp_r187___result;
    }
    static inline O __pp_r188__(const I &x_0, const I &x_1) {
       O __pp_r188___result;
       __pp_r188___result = (7.0L/14400.0L)*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 107.0L/14400.0L*pow(x_0, 5) + (1.0L/120.0L)*pow(x_0, 4)*pow(x_1, 2) - 19.0L/480.0L*pow(x_0, 4)*x_1 + (163.0L/3840.0L)*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) - 4.0L/45.0L*pow(x_0, 3)*pow(x_1, 2) + (269.0L/1440.0L)*pow(x_0, 3)*x_1 - 1799.0L/17280.0L*pow(x_0, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 4) - 13.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (53.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 373.0L/960.0L*pow(x_0, 2)*x_1 + (101.0L/1024.0L)*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) - 53.0L/720.0L*x_0*pow(x_1, 4) + (13.0L/45.0L)*x_0*pow(x_1, 3) - 779.0L/1440.0L*x_0*pow(x_1, 2) + (1117.0L/2304.0L)*x_0*x_1 - 40487.0L/230400.0L*x_0 + (7.0L/3600.0L)*pow(x_1, 6) - 11.0L/600.0L*pow(x_1, 5) + (13.0L/240.0L)*pow(x_1, 4) - 31.0L/1440.0L*pow(x_1, 3) - 31.0L/256.0L*pow(x_1, 2) + (2021.0L/38400.0L)*x_1 + 151043.0L/921600.0L;
       return __pp_r188___result;
    }
    static inline O __pp_r189__(const I &x_0, const I &x_1) {
       O __pp_r189___result;
       __pp_r189___result = (7.0L/14400.0L)*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 167.0L/14400.0L*pow(x_0, 5) + (1.0L/120.0L)*pow(x_0, 4)*pow(x_1, 2) - 19.0L/480.0L*pow(x_0, 4)*x_1 + (203.0L/3840.0L)*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) - 4.0L/45.0L*pow(x_0, 3)*pow(x_1, 2) + (269.0L/1440.0L)*pow(x_0, 3)*x_1 - 1979.0L/17280.0L*pow(x_0, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 4) - 13.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (53.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 373.0L/960.0L*pow(x_0, 2)*x_1 + (319.0L/3072.0L)*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) - 53.0L/720.0L*x_0*pow(x_1, 4) + (13.0L/45.0L)*x_0*pow(x_1, 3) - 779.0L/1440.0L*x_0*pow(x_1, 2) + (1117.0L/2304.0L)*x_0*x_1 - 40787.0L/230400.0L*x_0 + (7.0L/3600.0L)*pow(x_1, 6) - 11.0L/600.0L*pow(x_1, 5) + (13.0L/240.0L)*pow(x_1, 4) - 31.0L/1440.0L*pow(x_1, 3) - 31.0L/256.0L*pow(x_1, 2) + (2021.0L/38400.0L)*x_1 + 151163.0L/921600.0L;
       return __pp_r189___result;
    }
    static inline O __pp_r190__(const I &x_0, const I &x_1) {
       O __pp_r190___result;
       __pp_r190___result = -43.0L/43200.0L*pow(x_0, 6) + (3.0L/400.0L)*pow(x_0, 5)*x_1 - 263.0L/14400.0L*pow(x_0, 5) + (1.0L/360.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/480.0L*pow(x_0, 4)*x_1 + (31.0L/768.0L)*pow(x_0, 4) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 3) - 19.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (61.0L/288.0L)*pow(x_0, 3)*x_1 - 439.0L/3456.0L*pow(x_0, 3) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/10.0L*pow(x_0, 2)*pow(x_1, 3) + (5.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 71.0L/192.0L*pow(x_0, 2)*x_1 + (1487.0L/15360.0L)*pow(x_0, 2) + (3.0L/400.0L)*x_0*pow(x_1, 5) - 109.0L/1440.0L*x_0*pow(x_1, 4) + (85.0L/288.0L)*x_0*pow(x_1, 3) - 317.0L/576.0L*x_0*pow(x_1, 2) + (2833.0L/5760.0L)*x_0*x_1 - 41273.0L/230400.0L*x_0 + (83.0L/43200.0L)*pow(x_1, 6) - 29.0L/1600.0L*pow(x_1, 5) + (41.0L/768.0L)*pow(x_1, 4) - 23.0L/1152.0L*pow(x_1, 3) - 629.0L/5120.0L*pow(x_1, 2) + (4123.0L/76800.0L)*x_1 + 3773.0L/23040.0L;
       return __pp_r190___result;
    }
    static inline O __pp_r191__(const I &x_0, const I &x_1) {
       O __pp_r191___result;
       __pp_r191___result = -11.0L/14400.0L*pow(x_0, 6) + (11.0L/1800.0L)*pow(x_0, 5)*x_1 - 223.0L/14400.0L*pow(x_0, 5) + (1.0L/160.0L)*pow(x_0, 4)*pow(x_1, 2) - 53.0L/1440.0L*pow(x_0, 4)*x_1 + (125.0L/2304.0L)*pow(x_0, 4) + (13.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 - 311.0L/3456.0L*pow(x_0, 3) + (7.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 23.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (19.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 277.0L/576.0L*pow(x_0, 2)*x_1 + (7021.0L/46080.0L)*pow(x_0, 2) + (11.0L/1800.0L)*x_0*pow(x_1, 5) - 89.0L/1440.0L*x_0*pow(x_1, 4) + (23.0L/96.0L)*x_0*pow(x_1, 3) - 253.0L/576.0L*x_0*pow(x_1, 2) + (731.0L/1920.0L)*x_0*x_1 - 31033.0L/230400.0L*x_0 + (31.0L/14400.0L)*pow(x_1, 6) - 301.0L/14400.0L*pow(x_1, 5) + (155.0L/2304.0L)*pow(x_1, 4) - 197.0L/3456.0L*pow(x_1, 3) - 3101.0L/46080.0L*pow(x_1, 2) + (2129.0L/230400.0L)*x_1 + 12343.0L/69120.0L;
       return __pp_r191___result;
    }
    static inline O __pp_r192__(const I &x_0, const I &x_1) {
       O __pp_r192___result;
       __pp_r192___result = (31.0L/43200.0L)*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 - 7.0L/1600.0L*pow(x_0, 5) + (17.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 31.0L/480.0L*pow(x_0, 4)*x_1 + (205.0L/2304.0L)*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/20.0L*pow(x_0, 3)*pow(x_1, 2) + (25.0L/288.0L)*pow(x_0, 3)*x_1 - 37.0L/1152.0L*pow(x_0, 3) + (23.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 17.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (43.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 109.0L/192.0L*pow(x_0, 2)*x_1 + (9521.0L/46080.0L)*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) - 7.0L/120.0L*x_0*pow(x_1, 4) + (2.0L/9.0L)*x_0*pow(x_1, 3) - 19.0L/48.0L*x_0*pow(x_1, 2) + (3761.0L/11520.0L)*x_0*x_1 - 8261.0L/76800.0L*x_0 + (47.0L/21600.0L)*pow(x_1, 6) - 17.0L/800.0L*pow(x_1, 5) + (5.0L/72.0L)*pow(x_1, 4) - 37.0L/576.0L*pow(x_1, 3) - 619.0L/11520.0L*pow(x_1, 2) - 83.0L/19200.0L*x_1 + 101869.0L/552960.0L;
       return __pp_r192___result;
    }
    static inline O __pp_r193__(const I &x_0, const I &x_1) {
       O __pp_r193___result;
       __pp_r193___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/960.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 - 25.0L/128.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/12.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/64.0L*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) - 1.0L/12.0L*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) - 75.0L/32.0L*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 - 3375.0L/1024.0L*x_0 + (1.0L/675.0L)*pow(x_1, 6) - 1.0L/30.0L*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) - 25.0L/16.0L*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) - 3375.0L/512.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r193___result;
    }
    static inline O __pp_r194__(const I &x_0, const I &x_1) {
       O __pp_r194___result;
       __pp_r194___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/960.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 - 25.0L/128.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/12.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/64.0L*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) - 1.0L/12.0L*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) - 75.0L/32.0L*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 - 3375.0L/1024.0L*x_0 + (1.0L/675.0L)*pow(x_1, 6) - 1.0L/30.0L*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) - 25.0L/16.0L*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) - 3375.0L/512.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r194___result;
    }
    static inline O __pp_r195__(const I &x_0, const I &x_1) {
       O __pp_r195___result;
       __pp_r195___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/960.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 - 25.0L/128.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/12.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/64.0L*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) - 1.0L/12.0L*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) - 75.0L/32.0L*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 - 3375.0L/1024.0L*x_0 + (1.0L/675.0L)*pow(x_1, 6) - 1.0L/30.0L*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) - 25.0L/16.0L*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) - 3375.0L/512.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r195___result;
    }
    static inline O __pp_r196__(const I &x_0, const I &x_1) {
       O __pp_r196___result;
       __pp_r196___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/960.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 - 25.0L/128.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/12.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/64.0L*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) - 1.0L/12.0L*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) - 75.0L/32.0L*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 - 3375.0L/1024.0L*x_0 + (1.0L/675.0L)*pow(x_1, 6) - 1.0L/30.0L*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) - 25.0L/16.0L*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) - 3375.0L/512.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r196___result;
    }
    static inline O __pp_r197__(const I &x_0, const I &x_1) {
       O __pp_r197___result;
       __pp_r197___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/960.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 - 25.0L/128.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/12.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/64.0L*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) - 1.0L/12.0L*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) - 75.0L/32.0L*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 - 3375.0L/1024.0L*x_0 + (1.0L/675.0L)*pow(x_1, 6) - 1.0L/30.0L*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) - 25.0L/16.0L*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) - 3375.0L/512.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r197___result;
    }
    static inline O __pp_r198__(const I &x_0, const I &x_1) {
       O __pp_r198___result;
       __pp_r198___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (13.0L/2880.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (5.0L/288.0L)*pow(x_0, 4)*x_1 - 83.0L/2304.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) - 19.0L/288.0L*pow(x_0, 3)*x_1 + (349.0L/3456.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/36.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 163.0L/576.0L*pow(x_0, 2)*x_1 + (1933.0L/9216.0L)*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) - 1.0L/18.0L*x_0*pow(x_1, 4) + (29.0L/72.0L)*x_0*pow(x_1, 3) - 419.0L/288.0L*x_0*pow(x_1, 2) + (6029.0L/2304.0L)*x_0*x_1 - 86339.0L/46080.0L*x_0 + (1.0L/800.0L)*pow(x_1, 6) - 1.0L/36.0L*pow(x_1, 5) + (37.0L/144.0L)*pow(x_1, 4) - 547.0L/432.0L*pow(x_1, 3) + (8077.0L/2304.0L)*pow(x_1, 2) - 119107.0L/23040.0L*x_1 + 1753837.0L/552960.0L;
       return __pp_r198___result;
    }
    static inline O __pp_r199__(const I &x_0, const I &x_1) {
       O __pp_r199___result;
       __pp_r199___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (91.0L/14400.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (17.0L/480.0L)*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (31.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 + (6139.0L/17280.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (1193.0L/960.0L)*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) + (4.0L/45.0L)*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) + (2299.0L/1440.0L)*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 + (310891.0L/230400.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (3.0L/100.0L)*pow(x_1, 5) - 17.0L/80.0L*pow(x_1, 4) + (553.0L/720.0L)*pow(x_1, 3) - 5579.0L/3840.0L*pow(x_1, 2) + (16339.0L/12800.0L)*x_1 - 294811.0L/921600.0L;
       return __pp_r199___result;
    }
    static inline O __pp_r200__(const I &x_0, const I &x_1) {
       O __pp_r200___result;
       __pp_r200___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/960.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 - 25.0L/128.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/12.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/64.0L*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) - 1.0L/12.0L*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) - 75.0L/32.0L*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 - 3375.0L/1024.0L*x_0 + (1.0L/675.0L)*pow(x_1, 6) - 1.0L/30.0L*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) - 25.0L/16.0L*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) - 3375.0L/512.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r200___result;
    }
    static inline O __pp_r201__(const I &x_0, const I &x_1) {
       O __pp_r201___result;
       __pp_r201___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (13.0L/2880.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (5.0L/288.0L)*pow(x_0, 4)*x_1 - 83.0L/2304.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) - 19.0L/288.0L*pow(x_0, 3)*x_1 + (349.0L/3456.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/36.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 163.0L/576.0L*pow(x_0, 2)*x_1 + (1933.0L/9216.0L)*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) - 1.0L/18.0L*x_0*pow(x_1, 4) + (29.0L/72.0L)*x_0*pow(x_1, 3) - 419.0L/288.0L*x_0*pow(x_1, 2) + (6029.0L/2304.0L)*x_0*x_1 - 86339.0L/46080.0L*x_0 + (1.0L/800.0L)*pow(x_1, 6) - 1.0L/36.0L*pow(x_1, 5) + (37.0L/144.0L)*pow(x_1, 4) - 547.0L/432.0L*pow(x_1, 3) + (8077.0L/2304.0L)*pow(x_1, 2) - 119107.0L/23040.0L*x_1 + 1753837.0L/552960.0L;
       return __pp_r201___result;
    }
    static inline O __pp_r202__(const I &x_0, const I &x_1) {
       O __pp_r202___result;
       __pp_r202___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (91.0L/14400.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (17.0L/480.0L)*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (31.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 + (6139.0L/17280.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (1193.0L/960.0L)*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) + (4.0L/45.0L)*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) + (2299.0L/1440.0L)*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 + (310891.0L/230400.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (3.0L/100.0L)*pow(x_1, 5) - 17.0L/80.0L*pow(x_1, 4) + (553.0L/720.0L)*pow(x_1, 3) - 5579.0L/3840.0L*pow(x_1, 2) + (16339.0L/12800.0L)*x_1 - 294811.0L/921600.0L;
       return __pp_r202___result;
    }
    static inline O __pp_r203__(const I &x_0, const I &x_1) {
       O __pp_r203___result;
       __pp_r203___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/960.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 - 25.0L/128.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/12.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/64.0L*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) - 1.0L/12.0L*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) - 75.0L/32.0L*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 - 3375.0L/1024.0L*x_0 + (1.0L/675.0L)*pow(x_1, 6) - 1.0L/30.0L*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) - 25.0L/16.0L*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) - 3375.0L/512.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r203___result;
    }
    static inline O __pp_r204__(const I &x_0, const I &x_1) {
       O __pp_r204___result;
       __pp_r204___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (13.0L/2880.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (5.0L/288.0L)*pow(x_0, 4)*x_1 - 83.0L/2304.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) - 19.0L/288.0L*pow(x_0, 3)*x_1 + (349.0L/3456.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/36.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 163.0L/576.0L*pow(x_0, 2)*x_1 + (1933.0L/9216.0L)*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) - 1.0L/18.0L*x_0*pow(x_1, 4) + (29.0L/72.0L)*x_0*pow(x_1, 3) - 419.0L/288.0L*x_0*pow(x_1, 2) + (6029.0L/2304.0L)*x_0*x_1 - 86339.0L/46080.0L*x_0 + (1.0L/800.0L)*pow(x_1, 6) - 1.0L/36.0L*pow(x_1, 5) + (37.0L/144.0L)*pow(x_1, 4) - 547.0L/432.0L*pow(x_1, 3) + (8077.0L/2304.0L)*pow(x_1, 2) - 119107.0L/23040.0L*x_1 + 1753837.0L/552960.0L;
       return __pp_r204___result;
    }
    static inline O __pp_r205__(const I &x_0, const I &x_1) {
       O __pp_r205___result;
       __pp_r205___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/960.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 - 25.0L/128.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/12.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/64.0L*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) - 1.0L/12.0L*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) - 75.0L/32.0L*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 - 3375.0L/1024.0L*x_0 + (1.0L/675.0L)*pow(x_1, 6) - 1.0L/30.0L*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) - 25.0L/16.0L*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) - 3375.0L/512.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r205___result;
    }
    static inline O __pp_r206__(const I &x_0, const I &x_1) {
       O __pp_r206___result;
       __pp_r206___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (13.0L/2880.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (5.0L/288.0L)*pow(x_0, 4)*x_1 - 83.0L/2304.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) - 19.0L/288.0L*pow(x_0, 3)*x_1 + (349.0L/3456.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/36.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 163.0L/576.0L*pow(x_0, 2)*x_1 + (1933.0L/9216.0L)*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) - 1.0L/18.0L*x_0*pow(x_1, 4) + (29.0L/72.0L)*x_0*pow(x_1, 3) - 419.0L/288.0L*x_0*pow(x_1, 2) + (6029.0L/2304.0L)*x_0*x_1 - 86339.0L/46080.0L*x_0 + (1.0L/800.0L)*pow(x_1, 6) - 1.0L/36.0L*pow(x_1, 5) + (37.0L/144.0L)*pow(x_1, 4) - 547.0L/432.0L*pow(x_1, 3) + (8077.0L/2304.0L)*pow(x_1, 2) - 119107.0L/23040.0L*x_1 + 1753837.0L/552960.0L;
       return __pp_r206___result;
    }
    static inline O __pp_r207__(const I &x_0, const I &x_1) {
       O __pp_r207___result;
       __pp_r207___result = -1.0L/800.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/144.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/144.0L)*pow(x_0, 4)*x_1 - 1.0L/288.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 3)*x_1 + (29.0L/432.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/36.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 83.0L/288.0L*pow(x_0, 2)*x_1 + (523.0L/2304.0L)*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) - 1.0L/18.0L*x_0*pow(x_1, 4) + (29.0L/72.0L)*x_0*pow(x_1, 3) - 419.0L/288.0L*x_0*pow(x_1, 2) + (377.0L/144.0L)*x_0*x_1 - 4327.0L/2304.0L*x_0 + (1.0L/800.0L)*pow(x_1, 6) - 1.0L/36.0L*pow(x_1, 5) + (37.0L/144.0L)*pow(x_1, 4) - 547.0L/432.0L*pow(x_1, 3) + (8077.0L/2304.0L)*pow(x_1, 2) - 11911.0L/2304.0L*x_1 + 10963.0L/3456.0L;
       return __pp_r207___result;
    }
    static inline O __pp_r208__(const I &x_0, const I &x_1) {
       O __pp_r208___result;
       __pp_r208___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (91.0L/14400.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (17.0L/480.0L)*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (31.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 + (6139.0L/17280.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (1193.0L/960.0L)*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) + (4.0L/45.0L)*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) + (2299.0L/1440.0L)*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 + (310891.0L/230400.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (3.0L/100.0L)*pow(x_1, 5) - 17.0L/80.0L*pow(x_1, 4) + (553.0L/720.0L)*pow(x_1, 3) - 5579.0L/3840.0L*pow(x_1, 2) + (16339.0L/12800.0L)*x_1 - 294811.0L/921600.0L;
       return __pp_r208___result;
    }
    static inline O __pp_r209__(const I &x_0, const I &x_1) {
       O __pp_r209___result;
       __pp_r209___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (91.0L/14400.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (17.0L/480.0L)*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (31.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 + (6139.0L/17280.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (1193.0L/960.0L)*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) + (4.0L/45.0L)*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) + (2299.0L/1440.0L)*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 + (310891.0L/230400.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (3.0L/100.0L)*pow(x_1, 5) - 17.0L/80.0L*pow(x_1, 4) + (553.0L/720.0L)*pow(x_1, 3) - 5579.0L/3840.0L*pow(x_1, 2) + (16339.0L/12800.0L)*x_1 - 294811.0L/921600.0L;
       return __pp_r209___result;
    }
    static inline O __pp_r210__(const I &x_0, const I &x_1) {
       O __pp_r210___result;
       __pp_r210___result = -7.0L/5400.0L*pow(x_0, 6) + (1.0L/400.0L)*pow(x_0, 5)*x_1 - 37.0L/7200.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/40.0L)*pow(x_0, 4)*x_1 - 21.0L/640.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (31.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 209.0L/720.0L*pow(x_0, 3)*x_1 + (2777.0L/8640.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (99.0L/80.0L)*pow(x_0, 2)*x_1 - 7777.0L/7680.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) + (4.0L/45.0L)*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) + (2299.0L/1440.0L)*x_0*pow(x_1, 2) - 13481.0L/5760.0L*x_0*x_1 + (154943.0L/115200.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) + (3.0L/100.0L)*pow(x_1, 5) - 17.0L/80.0L*pow(x_1, 4) + (553.0L/720.0L)*pow(x_1, 3) - 5579.0L/3840.0L*pow(x_1, 2) + (12253.0L/9600.0L)*x_1 - 147203.0L/460800.0L;
       return __pp_r210___result;
    }
    static inline O __pp_r211__(const I &x_0, const I &x_1) {
       O __pp_r211___result;
       __pp_r211___result = -1.0L/1200.0L*pow(x_0, 6) + (19.0L/3600.0L)*pow(x_0, 5)*x_1 - 97.0L/7200.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/60.0L*pow(x_0, 4)*x_1 + (19.0L/640.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 29.0L/720.0L*pow(x_0, 3)*x_1 + (617.0L/8640.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) - 31.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (39.0L/80.0L)*pow(x_0, 2)*x_1 - 3457.0L/7680.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (17.0L/360.0L)*x_0*pow(x_1, 4) - 103.0L/360.0L*x_0*pow(x_1, 3) + (1219.0L/1440.0L)*x_0*pow(x_1, 2) - 7001.0L/5760.0L*x_0*x_1 + (77183.0L/115200.0L)*x_0 - 1.0L/800.0L*pow(x_1, 6) + (13.0L/600.0L)*pow(x_1, 5) - 3.0L/20.0L*pow(x_1, 4) + (373.0L/720.0L)*pow(x_1, 3) - 3419.0L/3840.0L*pow(x_1, 2) + (5773.0L/9600.0L)*x_1 + 8317.0L/460800.0L;
       return __pp_r211___result;
    }
    static inline O __pp_r212__(const I &x_0, const I &x_1) {
       O __pp_r212___result;
       __pp_r212___result = -7.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 - 41.0L/2880.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 7.0L/288.0L*pow(x_0, 4)*x_1 + (463.0L/11520.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/36.0L*pow(x_0, 3)*pow(x_1, 2) + (7.0L/160.0L)*pow(x_0, 3)*x_1 - 97.0L/17280.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/36.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) + (73.0L/2880.0L)*pow(x_0, 2)*x_1 - 6101.0L/46080.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 1.0L/72.0L*x_0*pow(x_1, 4) + (1.0L/20.0L)*x_0*pow(x_1, 3) - 7.0L/90.0L*x_0*pow(x_1, 2) + (71.0L/1280.0L)*x_0*x_1 - 1337.0L/46080.0L*x_0 + (1.0L/4320.0L)*pow(x_1, 6) - 1.0L/360.0L*pow(x_1, 5) + (13.0L/720.0L)*pow(x_1, 4) - 53.0L/540.0L*pow(x_1, 3) + (137.0L/360.0L)*pow(x_1, 2) - 3671.0L/4608.0L*x_1 + 1821463.0L/2764800.0L;
       return __pp_r212___result;
    }
    static inline O __pp_r213__(const I &x_0, const I &x_1) {
       O __pp_r213___result;
       __pp_r213___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (1.0L/240.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/48.0L*pow(x_0, 4)*x_1 + (1.0L/32.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) - 1.0L/8.0L*pow(x_0, 3)*x_1 + (1.0L/8.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (3.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 3.0L/8.0L*pow(x_0, 2)*x_1 + (9.0L/32.0L)*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) - 5.0L/96.0L*x_0*pow(x_1, 4) + (37.0L/96.0L)*x_0*pow(x_1, 3) - 271.0L/192.0L*x_0*pow(x_1, 2) + (1969.0L/768.0L)*x_0*x_1 - 2843.0L/1536.0L*x_0 + (11.0L/8640.0L)*pow(x_1, 6) - 9.0L/320.0L*pow(x_1, 5) + (199.0L/768.0L)*pow(x_1, 4) - 163.0L/128.0L*pow(x_1, 3) + (10811.0L/3072.0L)*pow(x_1, 2) - 15923.0L/3072.0L*x_1 + 39049.0L/12288.0L;
       return __pp_r213___result;
    }
    static inline O __pp_r214__(const I &x_0, const I &x_1) {
       O __pp_r214___result;
       __pp_r214___result = (1.0L/5400.0L)*pow(x_0, 6) - 7.0L/3600.0L*pow(x_0, 5)*x_1 + (43.0L/7200.0L)*pow(x_0, 5) + (1.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/360.0L*pow(x_0, 4)*x_1 + (11.0L/5760.0L)*pow(x_0, 4) - 13.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (41.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 259.0L/720.0L*pow(x_0, 3)*x_1 + (3277.0L/8640.0L)*pow(x_0, 3) - 11.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (37.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 31.0L/60.0L*pow(x_0, 2)*pow(x_1, 2) + (1657.0L/1440.0L)*pow(x_0, 2)*x_1 - 22081.0L/23040.0L*pow(x_0, 2) - 11.0L/1800.0L*x_0*pow(x_1, 5) + (133.0L/1440.0L)*x_0*pow(x_1, 4) - 797.0L/1440.0L*x_0*pow(x_1, 3) + (4723.0L/2880.0L)*x_0*pow(x_1, 2) - 27587.0L/11520.0L*x_0*x_1 + (39517.0L/28800.0L)*x_0 - 73.0L/43200.0L*pow(x_1, 6) + (427.0L/14400.0L)*pow(x_1, 5) - 2423.0L/11520.0L*pow(x_1, 4) + (13147.0L/17280.0L)*pow(x_1, 3) - 66323.0L/46080.0L*pow(x_1, 2) + (290947.0L/230400.0L)*x_1 - 867593.0L/2764800.0L;
       return __pp_r214___result;
    }
    static inline O __pp_r215__(const I &x_0, const I &x_1) {
       O __pp_r215___result;
       __pp_r215___result = (1.0L/240.0L)*x_0*pow(x_1, 5) - 7.0L/96.0L*x_0*pow(x_1, 4) + (49.0L/96.0L)*x_0*pow(x_1, 3) - 343.0L/192.0L*x_0*pow(x_1, 2) + (2401.0L/768.0L)*x_0*x_1 - 16807.0L/7680.0L*x_0 + (1.0L/960.0L)*pow(x_1, 6) - 23.0L/960.0L*pow(x_1, 5) + (175.0L/768.0L)*pow(x_1, 4) - 147.0L/128.0L*pow(x_1, 3) + (9947.0L/3072.0L)*pow(x_1, 2) - 74431.0L/15360.0L*x_1 + 184877.0L/61440.0L;
       return __pp_r215___result;
    }
    static inline O __pp_r216__(const I &x_0, const I &x_1) {
       O __pp_r216___result;
       __pp_r216___result = -1.0L/21600.0L*pow(x_0, 6) - 1.0L/1800.0L*pow(x_0, 5)*x_1 + (13.0L/7200.0L)*pow(x_0, 5) - 1.0L/360.0L*pow(x_0, 4)*pow(x_1, 2) + (13.0L/720.0L)*pow(x_0, 4)*x_1 - 169.0L/5760.0L*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) + (13.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) - 169.0L/720.0L*pow(x_0, 3)*x_1 + (2197.0L/8640.0L)*pow(x_0, 3) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 4) + (13.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) - 169.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) + (2197.0L/1440.0L)*pow(x_0, 2)*x_1 - 28561.0L/23040.0L*pow(x_0, 2) - 17.0L/3600.0L*x_0*pow(x_1, 5) + (103.0L/1440.0L)*x_0*pow(x_1, 4) - 617.0L/1440.0L*x_0*pow(x_1, 3) + (3643.0L/2880.0L)*x_0*pow(x_1, 2) - 21107.0L/11520.0L*x_0*x_1 + (29797.0L/28800.0L)*x_0 - 83.0L/43200.0L*pow(x_1, 6) + (487.0L/14400.0L)*pow(x_1, 5) - 2783.0L/11520.0L*pow(x_1, 4) + (15307.0L/17280.0L)*pow(x_1, 3) - 79283.0L/46080.0L*pow(x_1, 2) + (368707.0L/230400.0L)*x_1 - 1334153.0L/2764800.0L;
       return __pp_r216___result;
    }
    static inline O __pp_r217__(const I &x_0, const I &x_1) {
       O __pp_r217___result;
       __pp_r217___result = (7.0L/10800.0L)*pow(x_0, 6) + (1.0L/1200.0L)*pow(x_0, 5)*x_1 - 17.0L/7200.0L*pow(x_0, 5) + (11.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 2.0L/45.0L*pow(x_0, 4)*x_1 + (371.0L/5760.0L)*pow(x_0, 4) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 3) + (11.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 79.0L/720.0L*pow(x_0, 3)*x_1 + (1117.0L/8640.0L)*pow(x_0, 3) - 1.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 17.0L/120.0L*pow(x_0, 2)*pow(x_1, 2) + (577.0L/1440.0L)*pow(x_0, 2)*x_1 - 9121.0L/23040.0L*pow(x_0, 2) - 1.0L/300.0L*x_0*pow(x_1, 5) + (73.0L/1440.0L)*x_0*pow(x_1, 4) - 437.0L/1440.0L*x_0*pow(x_1, 3) + (2563.0L/2880.0L)*x_0*pow(x_1, 2) - 14627.0L/11520.0L*x_0*x_1 + (20077.0L/28800.0L)*x_0 - 53.0L/43200.0L*pow(x_1, 6) + (307.0L/14400.0L)*pow(x_1, 5) - 1703.0L/11520.0L*pow(x_1, 4) + (8827.0L/17280.0L)*pow(x_1, 3) - 40403.0L/46080.0L*pow(x_1, 2) + (135427.0L/230400.0L)*x_1 + 65527.0L/2764800.0L;
       return __pp_r217___result;
    }
    static inline O __pp_r218__(const I &x_0, const I &x_1) {
       O __pp_r218___result;
       __pp_r218___result = (29.0L/43200.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 - 1.0L/320.0L*pow(x_0, 5) + (13.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/96.0L*pow(x_0, 4)*x_1 + (863.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 37.0L/1440.0L*pow(x_0, 3)*x_1 + (301.0L/5760.0L)*pow(x_0, 3) + (7.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (53.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 59.0L/960.0L*pow(x_0, 2)*x_1 - 3601.0L/46080.0L*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 1.0L/96.0L*x_0*pow(x_1, 4) + (47.0L/1440.0L)*x_0*pow(x_1, 3) - 11.0L/320.0L*x_0*pow(x_1, 2) + (7.0L/5760.0L)*x_0*x_1 - 29.0L/15360.0L*x_0 + (11.0L/43200.0L)*pow(x_1, 6) - 1.0L/320.0L*pow(x_1, 5) + (233.0L/11520.0L)*pow(x_1, 4) - 607.0L/5760.0L*pow(x_1, 3) + (18161.0L/46080.0L)*pow(x_1, 2) - 2489.0L/3072.0L*x_1 + 57409.0L/86400.0L;
       return __pp_r218___result;
    }
    static inline O __pp_r219__(const I &x_0, const I &x_1) {
       O __pp_r219___result;
       __pp_r219___result = (1.0L/2400.0L)*pow(x_0, 6) + (1.0L/450.0L)*pow(x_0, 5)*x_1 - 47.0L/7200.0L*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) - 17.0L/720.0L*pow(x_0, 4)*x_1 + (191.0L/5760.0L)*pow(x_0, 4) + (1.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (11.0L/720.0L)*pow(x_0, 3)*x_1 + (37.0L/8640.0L)*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) + (11.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 79.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) + (1117.0L/1440.0L)*pow(x_0, 2)*x_1 - 15601.0L/23040.0L*pow(x_0, 2) - 7.0L/3600.0L*x_0*pow(x_1, 5) + (43.0L/1440.0L)*x_0*pow(x_1, 4) - 257.0L/1440.0L*x_0*pow(x_1, 3) + (1483.0L/2880.0L)*x_0*pow(x_1, 2) - 8147.0L/11520.0L*x_0*x_1 + (10357.0L/28800.0L)*x_0 - 7.0L/4800.0L*pow(x_1, 6) + (367.0L/14400.0L)*pow(x_1, 5) - 2063.0L/11520.0L*pow(x_1, 4) + (10987.0L/17280.0L)*pow(x_1, 3) - 53363.0L/46080.0L*pow(x_1, 2) + (213187.0L/230400.0L)*x_1 - 401033.0L/2764800.0L;
       return __pp_r219___result;
    }
    static inline O __pp_r220__(const I &x_0, const I &x_1) {
       O __pp_r220___result;
       __pp_r220___result = (19.0L/43200.0L)*pow(x_0, 6) + (1.0L/400.0L)*pow(x_0, 5)*x_1 - 7.0L/960.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/32.0L*pow(x_0, 4)*x_1 + (503.0L/11520.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (143.0L/1440.0L)*pow(x_0, 3)*x_1 - 419.0L/5760.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 37.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (301.0L/960.0L)*pow(x_0, 2)*x_1 - 16561.0L/46080.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) - 1.0L/32.0L*x_0*pow(x_1, 4) + (227.0L/1440.0L)*x_0*pow(x_1, 3) - 131.0L/320.0L*x_0*pow(x_1, 2) + (3247.0L/5760.0L)*x_0*x_1 - 5213.0L/15360.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (1.0L/960.0L)*pow(x_1, 5) - 127.0L/11520.0L*pow(x_1, 4) + (113.0L/5760.0L)*pow(x_1, 3) + (5201.0L/46080.0L)*pow(x_1, 2) - 7261.0L/15360.0L*x_1 + 42829.0L/86400.0L;
       return __pp_r220___result;
    }
    static inline O __pp_r221__(const I &x_0, const I &x_1) {
       O __pp_r221___result;
       __pp_r221___result = (1.0L/240.0L)*x_0*pow(x_1, 5) - 7.0L/96.0L*x_0*pow(x_1, 4) + (49.0L/96.0L)*x_0*pow(x_1, 3) - 343.0L/192.0L*x_0*pow(x_1, 2) + (2401.0L/768.0L)*x_0*x_1 - 16807.0L/7680.0L*x_0 + (1.0L/960.0L)*pow(x_1, 6) - 23.0L/960.0L*pow(x_1, 5) + (175.0L/768.0L)*pow(x_1, 4) - 147.0L/128.0L*pow(x_1, 3) + (9947.0L/3072.0L)*pow(x_1, 2) - 74431.0L/15360.0L*x_1 + 184877.0L/61440.0L;
       return __pp_r221___result;
    }
    static inline O __pp_r222__(const I &x_0, const I &x_1) {
       O __pp_r222___result;
       __pp_r222___result = -1.0L/21600.0L*pow(x_0, 6) - 1.0L/1800.0L*pow(x_0, 5)*x_1 + (13.0L/7200.0L)*pow(x_0, 5) - 1.0L/360.0L*pow(x_0, 4)*pow(x_1, 2) + (13.0L/720.0L)*pow(x_0, 4)*x_1 - 169.0L/5760.0L*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) + (13.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) - 169.0L/720.0L*pow(x_0, 3)*x_1 + (2197.0L/8640.0L)*pow(x_0, 3) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 4) + (13.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) - 169.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) + (2197.0L/1440.0L)*pow(x_0, 2)*x_1 - 28561.0L/23040.0L*pow(x_0, 2) - 17.0L/3600.0L*x_0*pow(x_1, 5) + (103.0L/1440.0L)*x_0*pow(x_1, 4) - 617.0L/1440.0L*x_0*pow(x_1, 3) + (3643.0L/2880.0L)*x_0*pow(x_1, 2) - 21107.0L/11520.0L*x_0*x_1 + (29797.0L/28800.0L)*x_0 - 83.0L/43200.0L*pow(x_1, 6) + (487.0L/14400.0L)*pow(x_1, 5) - 2783.0L/11520.0L*pow(x_1, 4) + (15307.0L/17280.0L)*pow(x_1, 3) - 79283.0L/46080.0L*pow(x_1, 2) + (368707.0L/230400.0L)*x_1 - 1334153.0L/2764800.0L;
       return __pp_r222___result;
    }
    static inline O __pp_r223__(const I &x_0, const I &x_1) {
       O __pp_r223___result;
       __pp_r223___result = (1.0L/2400.0L)*pow(x_0, 6) + (1.0L/450.0L)*pow(x_0, 5)*x_1 - 47.0L/7200.0L*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) - 17.0L/720.0L*pow(x_0, 4)*x_1 + (191.0L/5760.0L)*pow(x_0, 4) + (1.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (11.0L/720.0L)*pow(x_0, 3)*x_1 + (37.0L/8640.0L)*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) + (11.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 79.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) + (1117.0L/1440.0L)*pow(x_0, 2)*x_1 - 15601.0L/23040.0L*pow(x_0, 2) - 7.0L/3600.0L*x_0*pow(x_1, 5) + (43.0L/1440.0L)*x_0*pow(x_1, 4) - 257.0L/1440.0L*x_0*pow(x_1, 3) + (1483.0L/2880.0L)*x_0*pow(x_1, 2) - 8147.0L/11520.0L*x_0*x_1 + (10357.0L/28800.0L)*x_0 - 7.0L/4800.0L*pow(x_1, 6) + (367.0L/14400.0L)*pow(x_1, 5) - 2063.0L/11520.0L*pow(x_1, 4) + (10987.0L/17280.0L)*pow(x_1, 3) - 53363.0L/46080.0L*pow(x_1, 2) + (213187.0L/230400.0L)*x_1 - 401033.0L/2764800.0L;
       return __pp_r223___result;
    }
    static inline O __pp_r224__(const I &x_0, const I &x_1) {
       O __pp_r224___result;
       __pp_r224___result = (19.0L/43200.0L)*pow(x_0, 6) + (1.0L/400.0L)*pow(x_0, 5)*x_1 - 7.0L/960.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/32.0L*pow(x_0, 4)*x_1 + (503.0L/11520.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (143.0L/1440.0L)*pow(x_0, 3)*x_1 - 419.0L/5760.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 37.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (301.0L/960.0L)*pow(x_0, 2)*x_1 - 16561.0L/46080.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) - 1.0L/32.0L*x_0*pow(x_1, 4) + (227.0L/1440.0L)*x_0*pow(x_1, 3) - 131.0L/320.0L*x_0*pow(x_1, 2) + (3247.0L/5760.0L)*x_0*x_1 - 5213.0L/15360.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (1.0L/960.0L)*pow(x_1, 5) - 127.0L/11520.0L*pow(x_1, 4) + (113.0L/5760.0L)*pow(x_1, 3) + (5201.0L/46080.0L)*pow(x_1, 2) - 7261.0L/15360.0L*x_1 + 42829.0L/86400.0L;
       return __pp_r224___result;
    }
    static inline O __pp_r225__(const I &x_0, const I &x_1) {
       O __pp_r225___result;
       __pp_r225___result = (1.0L/400.0L)*pow(x_0, 6) - 11.0L/1800.0L*pow(x_0, 5)*x_1 + (13.0L/450.0L)*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) - 2.0L/45.0L*pow(x_0, 4)*x_1 + (41.0L/360.0L)*pow(x_0, 4) + (1.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) - 1.0L/180.0L*pow(x_0, 3)*x_1 + (89.0L/1080.0L)*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) + (11.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 79.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) + (551.0L/720.0L)*pow(x_0, 2)*x_1 - 3679.0L/5760.0L*pow(x_0, 2) - 7.0L/3600.0L*x_0*pow(x_1, 5) + (43.0L/1440.0L)*x_0*pow(x_1, 4) - 257.0L/1440.0L*x_0*pow(x_1, 3) + (1483.0L/2880.0L)*x_0*pow(x_1, 2) - 8177.0L/11520.0L*x_0*x_1 + (42523.0L/115200.0L)*x_0 - 7.0L/4800.0L*pow(x_1, 6) + (367.0L/14400.0L)*pow(x_1, 5) - 2063.0L/11520.0L*pow(x_1, 4) + (10987.0L/17280.0L)*pow(x_1, 3) - 53363.0L/46080.0L*pow(x_1, 2) + (213127.0L/230400.0L)*x_1 - 398423.0L/2764800.0L;
       return __pp_r225___result;
    }
    static inline O __pp_r226__(const I &x_0, const I &x_1) {
       O __pp_r226___result;
       __pp_r226___result = (109.0L/43200.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 + (9.0L/320.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/96.0L*pow(x_0, 4)*x_1 + (1433.0L/11520.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) + (113.0L/1440.0L)*pow(x_0, 3)*x_1 + (31.0L/5760.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 37.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (97.0L/320.0L)*pow(x_0, 2)*x_1 - 14791.0L/46080.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) - 1.0L/32.0L*x_0*pow(x_1, 4) + (227.0L/1440.0L)*x_0*pow(x_1, 3) - 131.0L/320.0L*x_0*pow(x_1, 2) + (101.0L/180.0L)*x_0*x_1 - 1689.0L/5120.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (1.0L/960.0L)*pow(x_1, 5) - 127.0L/11520.0L*pow(x_1, 4) + (113.0L/5760.0L)*pow(x_1, 3) + (5201.0L/46080.0L)*pow(x_1, 2) - 1453.0L/3072.0L*x_1 + 686569.0L/1382400.0L;
       return __pp_r226___result;
    }
    static inline O __pp_r227__(const I &x_0, const I &x_1) {
       O __pp_r227___result;
       __pp_r227___result = -1.0L/2160.0L*pow(x_0, 6) + (1.0L/360.0L)*pow(x_0, 5)*x_1 - 1.0L/90.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/18.0L)*pow(x_0, 4)*x_1 - 1.0L/9.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/9.0L*pow(x_0, 3)*pow(x_1, 2) + (4.0L/9.0L)*pow(x_0, 3)*x_1 - 16.0L/27.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/9.0L)*pow(x_0, 2)*pow(x_1, 3) - 2.0L/3.0L*pow(x_0, 2)*pow(x_1, 2) + (16.0L/9.0L)*pow(x_0, 2)*x_1 - 16.0L/9.0L*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) + (5.0L/288.0L)*x_0*pow(x_1, 4) - 19.0L/288.0L*x_0*pow(x_1, 3) + (5.0L/576.0L)*x_0*pow(x_1, 2) + (989.0L/2304.0L)*x_0*x_1 - 3023.0L/4608.0L*x_0 - 13.0L/8640.0L*pow(x_1, 6) + (77.0L/2880.0L)*pow(x_1, 5) - 445.0L/2304.0L*pow(x_1, 4) + (2489.0L/3456.0L)*pow(x_1, 3) - 13297.0L/9216.0L*pow(x_1, 2) + (13249.0L/9216.0L)*x_1 - 292261.0L/552960.0L;
       return __pp_r227___result;
    }
    static inline O __pp_r228__(const I &x_0, const I &x_1) {
       O __pp_r228___result;
       __pp_r228___result = -1.0L/240.0L*x_0*pow(x_1, 5) + (7.0L/96.0L)*x_0*pow(x_1, 4) - 49.0L/96.0L*x_0*pow(x_1, 3) + (343.0L/192.0L)*x_0*pow(x_1, 2) - 2401.0L/768.0L*x_0*x_1 + (16807.0L/7680.0L)*x_0 - 1.0L/960.0L*pow(x_1, 6) + (1.0L/64.0L)*pow(x_1, 5) - 21.0L/256.0L*pow(x_1, 4) + (49.0L/384.0L)*pow(x_1, 3) + (343.0L/1024.0L)*pow(x_1, 2) - 7203.0L/5120.0L*x_1 + 16807.0L/12288.0L;
       return __pp_r228___result;
    }
    static inline O __pp_r229__(const I &x_0, const I &x_1) {
       O __pp_r229___result;
       __pp_r229___result = -1.0L/240.0L*x_0*pow(x_1, 5) + (7.0L/96.0L)*x_0*pow(x_1, 4) - 49.0L/96.0L*x_0*pow(x_1, 3) + (343.0L/192.0L)*x_0*pow(x_1, 2) - 2401.0L/768.0L*x_0*x_1 + (16807.0L/7680.0L)*x_0 - 1.0L/960.0L*pow(x_1, 6) + (1.0L/64.0L)*pow(x_1, 5) - 21.0L/256.0L*pow(x_1, 4) + (49.0L/384.0L)*pow(x_1, 3) + (343.0L/1024.0L)*pow(x_1, 2) - 7203.0L/5120.0L*x_1 + 16807.0L/12288.0L;
       return __pp_r229___result;
    }
    static inline O __pp_r230__(const I &x_0, const I &x_1) {
       O __pp_r230___result;
       __pp_r230___result = -19.0L/43200.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 19.0L/1600.0L*pow(x_0, 5) - 1.0L/180.0L*pow(x_0, 4)*pow(x_1, 2) + (23.0L/480.0L)*pow(x_0, 4)*x_1 - 1159.0L/11520.0L*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) - 17.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (761.0L/1440.0L)*pow(x_0, 3)*x_1 - 3857.0L/5760.0L*pow(x_0, 3) - 1.0L/720.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/20.0L)*pow(x_0, 2)*pow(x_1, 3) - 199.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (421.0L/320.0L)*pow(x_0, 2)*x_1 - 67279.0L/46080.0L*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) - 7.0L/160.0L*x_0*pow(x_1, 4) + (389.0L/1440.0L)*x_0*pow(x_1, 3) - 293.0L/320.0L*x_0*pow(x_1, 2) + (9793.0L/5760.0L)*x_0*x_1 - 34689.0L/25600.0L*x_0 - 1.0L/43200.0L*pow(x_1, 6) + (11.0L/4800.0L)*pow(x_1, 5) - 289.0L/11520.0L*pow(x_1, 4) + (599.0L/5760.0L)*pow(x_1, 3) - 7921.0L/46080.0L*pow(x_1, 2) + (3041.0L/76800.0L)*x_1 + 19391.0L/172800.0L;
       return __pp_r230___result;
    }
    static inline O __pp_r231__(const I &x_0, const I &x_1) {
       O __pp_r231___result;
       __pp_r231___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 11.0L/14400.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/1440.0L*pow(x_0, 4)*x_1 + (121.0L/11520.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 11.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (121.0L/1440.0L)*pow(x_0, 3)*x_1 - 1331.0L/17280.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 11.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (121.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 1331.0L/2880.0L*pow(x_0, 2)*x_1 + (14641.0L/46080.0L)*pow(x_0, 2) + (1.0L/3600.0L)*x_0*pow(x_1, 5) + (17.0L/1440.0L)*x_0*pow(x_1, 4) - 251.0L/1440.0L*x_0*pow(x_1, 3) + (2483.0L/2880.0L)*x_0*pow(x_1, 2) - 10687.0L/5760.0L*x_0*x_1 + (343159.0L/230400.0L)*x_0 + (19.0L/43200.0L)*pow(x_1, 6) - 127.0L/14400.0L*pow(x_1, 5) + (991.0L/11520.0L)*pow(x_1, 4) - 8443.0L/17280.0L*pow(x_1, 3) + (73999.0L/46080.0L)*pow(x_1, 2) - 646237.0L/230400.0L*x_1 + 347071.0L/172800.0L;
       return __pp_r231___result;
    }
    static inline O __pp_r232__(const I &x_0, const I &x_1) {
       O __pp_r232___result;
       __pp_r232___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 11.0L/14400.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/1440.0L*pow(x_0, 4)*x_1 + (121.0L/11520.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 11.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (121.0L/1440.0L)*pow(x_0, 3)*x_1 - 1331.0L/17280.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 11.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (121.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 1331.0L/2880.0L*pow(x_0, 2)*x_1 + (14641.0L/46080.0L)*pow(x_0, 2) + (1.0L/3600.0L)*x_0*pow(x_1, 5) + (17.0L/1440.0L)*x_0*pow(x_1, 4) - 251.0L/1440.0L*x_0*pow(x_1, 3) + (2483.0L/2880.0L)*x_0*pow(x_1, 2) - 10687.0L/5760.0L*x_0*x_1 + (343159.0L/230400.0L)*x_0 + (19.0L/43200.0L)*pow(x_1, 6) - 127.0L/14400.0L*pow(x_1, 5) + (991.0L/11520.0L)*pow(x_1, 4) - 8443.0L/17280.0L*pow(x_1, 3) + (73999.0L/46080.0L)*pow(x_1, 2) - 646237.0L/230400.0L*x_1 + 347071.0L/172800.0L;
       return __pp_r232___result;
    }
    static inline O __pp_r233__(const I &x_0, const I &x_1) {
       O __pp_r233___result;
       __pp_r233___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 - 1.0L/14400.0L*pow(x_0, 5) - 41.0L/1440.0L*pow(x_0, 4)*x_1 + (17.0L/768.0L)*pow(x_0, 4) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/96.0L)*pow(x_0, 3)*x_1 - 77.0L/3456.0L*pow(x_0, 3) - 7.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 31.0L/576.0L*pow(x_0, 2)*x_1 - 887.0L/15360.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 4) + (1.0L/24.0L)*x_0*pow(x_1, 3) - 23.0L/288.0L*x_0*pow(x_1, 2) + (71.0L/1280.0L)*x_0*x_1 - 3601.0L/230400.0L*x_0 - 71.0L/7200.0L*pow(x_1, 5) + (5.0L/96.0L)*pow(x_1, 4) - 121.0L/1728.0L*pow(x_1, 3) - 3.0L/160.0L*pow(x_1, 2) - 839.0L/28800.0L*x_1 + 34433.0L/184320.0L;
       return __pp_r233___result;
    }
    static inline O __pp_r234__(const I &x_0, const I &x_1) {
       O __pp_r234___result;
       __pp_r234___result = (19.0L/43200.0L)*pow(x_0, 6) - 1.0L/3600.0L*pow(x_0, 5)*x_1 - 11.0L/4800.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/480.0L*pow(x_0, 4)*x_1 + (271.0L/11520.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (71.0L/1440.0L)*pow(x_0, 3)*x_1 - 131.0L/5760.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) + (31.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 17.0L/320.0L*pow(x_0, 2)*x_1 - 2657.0L/46080.0L*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/160.0L*x_0*pow(x_1, 4) + (59.0L/1440.0L)*x_0*pow(x_1, 3) - 77.0L/960.0L*x_0*pow(x_1, 2) + (319.0L/5760.0L)*x_0*x_1 - 1201.0L/76800.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 47.0L/4800.0L*pow(x_1, 5) + (601.0L/11520.0L)*pow(x_1, 4) - 403.0L/5760.0L*pow(x_1, 3) - 863.0L/46080.0L*pow(x_1, 2) - 2237.0L/76800.0L*x_1 + 32281.0L/172800.0L;
       return __pp_r234___result;
    }
    static inline O __pp_r235__(const I &x_0, const I &x_1) {
       O __pp_r235___result;
       __pp_r235___result = (109.0L/43200.0L)*pow(x_0, 6) - 31.0L/3600.0L*pow(x_0, 5)*x_1 - 1.0L/4800.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/480.0L*pow(x_0, 4)*x_1 + (121.0L/11520.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 - 41.0L/5760.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) + (31.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 41.0L/960.0L*pow(x_0, 2)*x_1 - 3047.0L/46080.0L*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/160.0L*x_0*pow(x_1, 4) + (59.0L/1440.0L)*x_0*pow(x_1, 3) - 77.0L/960.0L*x_0*pow(x_1, 2) + (19.0L/360.0L)*x_0*x_1 - 1031.0L/76800.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 47.0L/4800.0L*pow(x_1, 5) + (601.0L/11520.0L)*pow(x_1, 4) - 403.0L/5760.0L*pow(x_1, 3) - 863.0L/46080.0L*pow(x_1, 2) - 739.0L/25600.0L*x_1 + 257933.0L/1382400.0L;
       return __pp_r235___result;
    }
    static inline O __pp_r236__(const I &x_0, const I &x_1) {
       O __pp_r236___result;
       __pp_r236___result = -23.0L/21600.0L*pow(x_0, 6) + (7.0L/1800.0L)*pow(x_0, 5)*x_1 + (1.0L/3600.0L)*pow(x_0, 5) - 1.0L/720.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/40.0L*pow(x_0, 4)*x_1 + (23.0L/1152.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 3)*x_1 - 13.0L/864.0L*pow(x_0, 3) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 2)*x_1 - 1643.0L/23040.0L*pow(x_0, 2) - 1.0L/225.0L*x_0*pow(x_1, 5) + (1.0L/45.0L)*x_0*pow(x_1, 4) - 1.0L/36.0L*x_0*pow(x_1, 3) + (1.0L/144.0L)*x_0*pow(x_1, 2) + (7.0L/5760.0L)*x_0*x_1 - 119.0L/57600.0L*x_0 - 1.0L/675.0L*pow(x_1, 6) + (1.0L/800.0L)*pow(x_1, 5) + (5.0L/288.0L)*pow(x_1, 4) - 7.0L/576.0L*pow(x_1, 3) - 841.0L/11520.0L*pow(x_1, 2) - 77.0L/38400.0L*x_1 + 50087.0L/276480.0L;
       return __pp_r236___result;
    }
    static inline O __pp_r237__(const I &x_0, const I &x_1) {
       O __pp_r237___result;
       __pp_r237___result = (1.0L/2400.0L)*pow(x_0, 6) - 1.0L/1800.0L*pow(x_0, 5)*x_1 - 7.0L/3600.0L*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) - 7.0L/360.0L*pow(x_0, 4)*x_1 + (41.0L/1920.0L)*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (23.0L/720.0L)*pow(x_0, 3)*x_1 - 67.0L/4320.0L*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) - 7.0L/720.0L*pow(x_0, 2)*x_1 - 547.0L/7680.0L*pow(x_0, 2) - 17.0L/3600.0L*x_0*pow(x_1, 5) + (31.0L/1440.0L)*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) + (19.0L/2880.0L)*x_0*pow(x_1, 2) + (13.0L/11520.0L)*x_0*x_1 - 239.0L/115200.0L*x_0 - 7.0L/4800.0L*pow(x_1, 6) + (19.0L/14400.0L)*pow(x_1, 5) + (67.0L/3840.0L)*pow(x_1, 4) - 209.0L/17280.0L*pow(x_1, 3) - 1121.0L/15360.0L*pow(x_1, 2) - 461.0L/230400.0L*x_1 + 166957.0L/921600.0L;
       return __pp_r237___result;
    }
    static inline O __pp_r238__(const I &x_0, const I &x_1) {
       O __pp_r238___result;
       __pp_r238___result = (1.0L/400.0L)*pow(x_0, 6) - 2.0L/225.0L*pow(x_0, 5)*x_1 + (1.0L/7200.0L)*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/720.0L)*pow(x_0, 4)*x_1 + (1.0L/120.0L)*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/90.0L)*pow(x_0, 3)*x_1 + (1.0L/8640.0L)*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/1440.0L)*pow(x_0, 2)*x_1 - 51.0L/640.0L*pow(x_0, 2) - 17.0L/3600.0L*x_0*pow(x_1, 5) + (31.0L/1440.0L)*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) + (19.0L/2880.0L)*x_0*pow(x_1, 2) - 17.0L/11520.0L*x_0*x_1 + (1.0L/7200.0L)*x_0 - 7.0L/4800.0L*pow(x_1, 6) + (19.0L/14400.0L)*pow(x_1, 5) + (67.0L/3840.0L)*pow(x_1, 4) - 209.0L/17280.0L*pow(x_1, 3) - 1121.0L/15360.0L*pow(x_1, 2) - 401.0L/230400.0L*x_1 + 166747.0L/921600.0L;
       return __pp_r238___result;
    }
    static inline O __pp_r239__(const I &x_0, const I &x_1) {
       O __pp_r239___result;
       __pp_r239___result = -19.0L/43200.0L*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 - 67.0L/14400.0L*pow(x_0, 5) - 1.0L/180.0L*pow(x_0, 4)*pow(x_1, 2) + (13.0L/1440.0L)*pow(x_0, 4)*x_1 + (89.0L/11520.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (49.0L/1440.0L)*pow(x_0, 3)*x_1 - 139.0L/17280.0L*pow(x_0, 3) - 1.0L/720.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (29.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 119.0L/2880.0L*pow(x_0, 2)*x_1 - 611.0L/9216.0L*pow(x_0, 2) + (1.0L/3600.0L)*x_0*pow(x_1, 5) - 11.0L/1440.0L*x_0*pow(x_1, 4) + (61.0L/1440.0L)*x_0*pow(x_1, 3) - 233.0L/2880.0L*x_0*pow(x_1, 2) + (61.0L/1152.0L)*x_0*x_1 - 3097.0L/230400.0L*x_0 - 1.0L/43200.0L*pow(x_1, 6) - 139.0L/14400.0L*pow(x_1, 5) + (599.0L/11520.0L)*pow(x_1, 4) - 1207.0L/17280.0L*pow(x_1, 3) - 173.0L/9216.0L*pow(x_1, 2) - 6649.0L/230400.0L*x_1 + 64483.0L/345600.0L;
       return __pp_r239___result;
    }
    static inline O __pp_r240__(const I &x_0, const I &x_1) {
       O __pp_r240___result;
       __pp_r240___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 47.0L/14400.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/480.0L)*pow(x_0, 4)*x_1 + (43.0L/3840.0L)*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (29.0L/1440.0L)*pow(x_0, 3)*x_1 - 59.0L/17280.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/40.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 53.0L/960.0L*pow(x_0, 2)*x_1 - 193.0L/3072.0L*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) - 1.0L/1440.0L*x_0*pow(x_1, 4) + (41.0L/1440.0L)*x_0*pow(x_1, 3) - 193.0L/2880.0L*x_0*pow(x_1, 2) + (53.0L/1152.0L)*x_0*x_1 - 2777.0L/230400.0L*x_0 + (1.0L/4800.0L)*pow(x_1, 6) - 53.0L/4800.0L*pow(x_1, 5) + (71.0L/1280.0L)*pow(x_1, 4) - 143.0L/1920.0L*pow(x_1, 3) - 47.0L/3072.0L*pow(x_1, 2) - 2323.0L/76800.0L*x_1 + 21521.0L/115200.0L;
       return __pp_r240___result;
    }
    static inline O __pp_r241__(const I &x_0, const I &x_1) {
       O __pp_r241___result;
       __pp_r241___result = -73.0L/43200.0L*pow(x_0, 6) + (1.0L/300.0L)*pow(x_0, 5)*x_1 - 143.0L/14400.0L*pow(x_0, 5) - 11.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (3.0L/160.0L)*pow(x_0, 4)*x_1 - 1.0L/768.0L*pow(x_0, 4) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/45.0L*pow(x_0, 3)*pow(x_1, 2) + (13.0L/288.0L)*pow(x_0, 3)*x_1 - 55.0L/3456.0L*pow(x_0, 3) + (1.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 7.0L/192.0L*pow(x_0, 2)*x_1 - 1073.0L/15360.0L*pow(x_0, 2) - 1.0L/1200.0L*x_0*pow(x_1, 5) - 1.0L/360.0L*x_0*pow(x_1, 4) + (5.0L/144.0L)*x_0*pow(x_1, 3) - 11.0L/144.0L*x_0*pow(x_1, 2) + (611.0L/11520.0L)*x_0*x_1 - 3263.0L/230400.0L*x_0 + (1.0L/5400.0L)*pow(x_1, 6) - 13.0L/1200.0L*pow(x_1, 5) + (7.0L/128.0L)*pow(x_1, 4) - 7.0L/96.0L*pow(x_1, 3) - 131.0L/7680.0L*pow(x_1, 2) - 1121.0L/38400.0L*x_1 + 6877.0L/36864.0L;
       return __pp_r241___result;
    }
    static inline O __pp_r242__(const I &x_0, const I &x_1) {
       O __pp_r242___result;
       __pp_r242___result = -1.0L/2160.0L*pow(x_0, 6) - 31.0L/7200.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/80.0L)*pow(x_0, 4)*x_1 + (1.0L/180.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 3)*x_1 - 7.0L/8640.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/480.0L)*pow(x_0, 2)*x_1 - 23.0L/288.0L*pow(x_0, 2) - 1.0L/240.0L*x_0*pow(x_1, 5) + (29.0L/1440.0L)*x_0*pow(x_1, 4) - 13.0L/480.0L*x_0*pow(x_1, 3) + (17.0L/2880.0L)*x_0*pow(x_1, 2) - 1.0L/768.0L*x_0*x_1 + (7.0L/57600.0L)*x_0 - 13.0L/8640.0L*pow(x_1, 6) + (7.0L/4800.0L)*pow(x_1, 5) + (199.0L/11520.0L)*pow(x_1, 4) - 23.0L/1920.0L*pow(x_1, 3) - 673.0L/9216.0L*pow(x_1, 2) - 133.0L/76800.0L*x_1 + 500239.0L/2764800.0L;
       return __pp_r242___result;
    }
    static inline O __pp_r243__(const I &x_0, const I &x_1) {
       O __pp_r243___result;
       __pp_r243___result = -1.0L/4320.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 7.0L/2400.0L*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/180.0L)*pow(x_0, 4)*x_1 + (13.0L/1440.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/360.0L)*pow(x_0, 3)*x_1 + (11.0L/2880.0L)*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/240.0L)*pow(x_0, 2)*pow(x_1, 2) - 17.0L/1440.0L*pow(x_0, 2)*x_1 - 11.0L/144.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) + (13.0L/480.0L)*x_0*pow(x_1, 4) - 59.0L/1440.0L*x_0*pow(x_1, 3) + (19.0L/960.0L)*x_0*pow(x_1, 2) - 19.0L/2304.0L*x_0*x_1 + (29.0L/19200.0L)*x_0 - 11.0L/8640.0L*pow(x_1, 6) + (1.0L/14400.0L)*pow(x_1, 5) + (239.0L/11520.0L)*pow(x_1, 4) - 287.0L/17280.0L*pow(x_1, 3) - 641.0L/9216.0L*pow(x_1, 2) - 719.0L/230400.0L*x_1 + 500879.0L/2764800.0L;
       return __pp_r243___result;
    }
    static inline O __pp_r244__(const I &x_0, const I &x_1) {
       O __pp_r244___result;
       __pp_r244___result = -37.0L/21600.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 23.0L/2400.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/45.0L)*pow(x_0, 4)*x_1 - 1.0L/288.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/36.0L)*pow(x_0, 3)*x_1 - 5.0L/576.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/144.0L)*pow(x_0, 2)*x_1 - 961.0L/11520.0L*pow(x_0, 2) - 19.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/40.0L)*x_0*pow(x_1, 4) - 5.0L/144.0L*x_0*pow(x_1, 3) + (1.0L/96.0L)*x_0*pow(x_1, 2) - 7.0L/5760.0L*x_0*x_1 - 23.0L/38400.0L*x_0 - 7.0L/5400.0L*pow(x_1, 6) + (1.0L/3600.0L)*pow(x_1, 5) + (23.0L/1152.0L)*pow(x_1, 4) - 13.0L/864.0L*pow(x_1, 3) - 1643.0L/23040.0L*pow(x_1, 2) - 119.0L/57600.0L*x_1 + 10003.0L/55296.0L;
       return __pp_r244___result;
    }
    static inline O __pp_r245__(const I &x_0, const I &x_1) {
       O __pp_r245___result;
       __pp_r245___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 47.0L/14400.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/480.0L)*pow(x_0, 4)*x_1 + (43.0L/3840.0L)*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (29.0L/1440.0L)*pow(x_0, 3)*x_1 - 59.0L/17280.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/40.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 53.0L/960.0L*pow(x_0, 2)*x_1 - 193.0L/3072.0L*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) - 23.0L/720.0L*x_0*pow(x_1, 4) + (11.0L/90.0L)*x_0*pow(x_1, 3) - 299.0L/1440.0L*x_0*pow(x_1, 2) + (349.0L/2304.0L)*x_0*x_1 - 10067.0L/230400.0L*x_0 + (1.0L/800.0L)*pow(x_1, 6) - 1.0L/100.0L*pow(x_1, 5) + (1.0L/80.0L)*pow(x_1, 4) + (43.0L/480.0L)*pow(x_1, 3) - 221.0L/768.0L*pow(x_1, 2) + (7141.0L/38400.0L)*x_1 + 110203.0L/921600.0L;
       return __pp_r245___result;
    }
    static inline O __pp_r246__(const I &x_0, const I &x_1) {
       O __pp_r246___result;
       __pp_r246___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 47.0L/14400.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/480.0L)*pow(x_0, 4)*x_1 + (43.0L/3840.0L)*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (29.0L/1440.0L)*pow(x_0, 3)*x_1 - 59.0L/17280.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/40.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 53.0L/960.0L*pow(x_0, 2)*x_1 - 193.0L/3072.0L*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) - 23.0L/720.0L*x_0*pow(x_1, 4) + (11.0L/90.0L)*x_0*pow(x_1, 3) - 299.0L/1440.0L*x_0*pow(x_1, 2) + (349.0L/2304.0L)*x_0*x_1 - 10067.0L/230400.0L*x_0 + (1.0L/800.0L)*pow(x_1, 6) - 1.0L/100.0L*pow(x_1, 5) + (1.0L/80.0L)*pow(x_1, 4) + (43.0L/480.0L)*pow(x_1, 3) - 221.0L/768.0L*pow(x_1, 2) + (7141.0L/38400.0L)*x_1 + 110203.0L/921600.0L;
       return __pp_r246___result;
    }
    static inline O __pp_r247__(const I &x_0, const I &x_1) {
       O __pp_r247___result;
       __pp_r247___result = -73.0L/43200.0L*pow(x_0, 6) + (1.0L/300.0L)*pow(x_0, 5)*x_1 - 143.0L/14400.0L*pow(x_0, 5) - 11.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (3.0L/160.0L)*pow(x_0, 4)*x_1 - 1.0L/768.0L*pow(x_0, 4) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/45.0L*pow(x_0, 3)*pow(x_1, 2) + (13.0L/288.0L)*pow(x_0, 3)*x_1 - 55.0L/3456.0L*pow(x_0, 3) + (1.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 7.0L/192.0L*pow(x_0, 2)*x_1 - 1073.0L/15360.0L*pow(x_0, 2) + (1.0L/300.0L)*x_0*pow(x_1, 5) - 49.0L/1440.0L*x_0*pow(x_1, 4) + (37.0L/288.0L)*x_0*pow(x_1, 3) - 125.0L/576.0L*x_0*pow(x_1, 2) + (913.0L/5760.0L)*x_0*x_1 - 10553.0L/230400.0L*x_0 + (53.0L/43200.0L)*pow(x_1, 6) - 47.0L/4800.0L*pow(x_1, 5) + (3.0L/256.0L)*pow(x_1, 4) + (35.0L/384.0L)*pow(x_1, 3) - 4447.0L/15360.0L*pow(x_1, 2) + (14363.0L/76800.0L)*x_1 + 2749.0L/23040.0L;
       return __pp_r247___result;
    }
    static inline O __pp_r248__(const I &x_0, const I &x_1) {
       O __pp_r248___result;
       __pp_r248___result = -73.0L/43200.0L*pow(x_0, 6) + (1.0L/300.0L)*pow(x_0, 5)*x_1 - 143.0L/14400.0L*pow(x_0, 5) - 11.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (3.0L/160.0L)*pow(x_0, 4)*x_1 - 1.0L/768.0L*pow(x_0, 4) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/45.0L*pow(x_0, 3)*pow(x_1, 2) + (13.0L/288.0L)*pow(x_0, 3)*x_1 - 55.0L/3456.0L*pow(x_0, 3) + (1.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 7.0L/192.0L*pow(x_0, 2)*x_1 - 1073.0L/15360.0L*pow(x_0, 2) + (1.0L/300.0L)*x_0*pow(x_1, 5) - 49.0L/1440.0L*x_0*pow(x_1, 4) + (37.0L/288.0L)*x_0*pow(x_1, 3) - 125.0L/576.0L*x_0*pow(x_1, 2) + (913.0L/5760.0L)*x_0*x_1 - 10553.0L/230400.0L*x_0 + (53.0L/43200.0L)*pow(x_1, 6) - 47.0L/4800.0L*pow(x_1, 5) + (3.0L/256.0L)*pow(x_1, 4) + (35.0L/384.0L)*pow(x_1, 3) - 4447.0L/15360.0L*pow(x_1, 2) + (14363.0L/76800.0L)*x_1 + 2749.0L/23040.0L;
       return __pp_r248___result;
    }
    static inline O __pp_r249__(const I &x_0, const I &x_1) {
       O __pp_r249___result;
       __pp_r249___result = -7.0L/4800.0L*pow(x_0, 6) + (7.0L/3600.0L)*pow(x_0, 5)*x_1 - 103.0L/14400.0L*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) + (7.0L/1440.0L)*pow(x_0, 4)*x_1 + (29.0L/2304.0L)*pow(x_0, 4) - 1.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 3)*x_1 + (73.0L/3456.0L)*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) - 2.0L/45.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 85.0L/576.0L*pow(x_0, 2)*x_1 - 659.0L/46080.0L*pow(x_0, 2) + (7.0L/3600.0L)*x_0*pow(x_1, 5) - 29.0L/1440.0L*x_0*pow(x_1, 4) + (7.0L/96.0L)*x_0*pow(x_1, 3) - 61.0L/576.0L*x_0*pow(x_1, 2) + (91.0L/1920.0L)*x_0*x_1 - 313.0L/230400.0L*x_0 + (7.0L/4800.0L)*pow(x_1, 6) - 181.0L/14400.0L*pow(x_1, 5) + (59.0L/2304.0L)*pow(x_1, 4) + (187.0L/3456.0L)*pow(x_1, 3) - 10781.0L/46080.0L*pow(x_1, 2) + (32849.0L/230400.0L)*x_1 + 9271.0L/69120.0L;
       return __pp_r249___result;
    }
    static inline O __pp_r250__(const I &x_0, const I &x_1) {
       O __pp_r250___result;
       __pp_r250___result = -7.0L/4800.0L*pow(x_0, 6) + (7.0L/3600.0L)*pow(x_0, 5)*x_1 - 103.0L/14400.0L*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) + (7.0L/1440.0L)*pow(x_0, 4)*x_1 + (29.0L/2304.0L)*pow(x_0, 4) - 1.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 3)*x_1 + (73.0L/3456.0L)*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) - 2.0L/45.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 85.0L/576.0L*pow(x_0, 2)*x_1 - 659.0L/46080.0L*pow(x_0, 2) + (7.0L/3600.0L)*x_0*pow(x_1, 5) - 29.0L/1440.0L*x_0*pow(x_1, 4) + (7.0L/96.0L)*x_0*pow(x_1, 3) - 61.0L/576.0L*x_0*pow(x_1, 2) + (91.0L/1920.0L)*x_0*x_1 - 313.0L/230400.0L*x_0 + (7.0L/4800.0L)*pow(x_1, 6) - 181.0L/14400.0L*pow(x_1, 5) + (59.0L/2304.0L)*pow(x_1, 4) + (187.0L/3456.0L)*pow(x_1, 3) - 10781.0L/46080.0L*pow(x_1, 2) + (32849.0L/230400.0L)*x_1 + 9271.0L/69120.0L;
       return __pp_r250___result;
    }
    static inline O __pp_r251__(const I &x_0, const I &x_1) {
       O __pp_r251___result;
       __pp_r251___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 + (19.0L/4800.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/480.0L*pow(x_0, 4)*x_1 + (109.0L/2304.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) - 23.0L/288.0L*pow(x_0, 3)*x_1 + (91.0L/1152.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 7.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 15.0L/64.0L*pow(x_0, 2)*x_1 + (1841.0L/46080.0L)*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 1.0L/60.0L*x_0*pow(x_1, 4) + (1.0L/18.0L)*x_0*pow(x_1, 3) - 1.0L/16.0L*x_0*pow(x_1, 2) - 79.0L/11520.0L*x_0*x_1 + (1979.0L/76800.0L)*x_0 + (1.0L/675.0L)*pow(x_1, 6) - 31.0L/2400.0L*pow(x_1, 5) + (1.0L/36.0L)*pow(x_1, 4) + (3.0L/64.0L)*pow(x_1, 3) - 2539.0L/11520.0L*pow(x_1, 2) + (2477.0L/19200.0L)*x_1 + 77293.0L/552960.0L;
       return __pp_r251___result;
    }
    static inline O __pp_r252__(const I &x_0, const I &x_1) {
       O __pp_r252___result;
       __pp_r252___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 + (19.0L/4800.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/480.0L*pow(x_0, 4)*x_1 + (109.0L/2304.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) - 23.0L/288.0L*pow(x_0, 3)*x_1 + (91.0L/1152.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 7.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 15.0L/64.0L*pow(x_0, 2)*x_1 + (1841.0L/46080.0L)*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 1.0L/60.0L*x_0*pow(x_1, 4) + (1.0L/18.0L)*x_0*pow(x_1, 3) - 1.0L/16.0L*x_0*pow(x_1, 2) - 79.0L/11520.0L*x_0*x_1 + (1979.0L/76800.0L)*x_0 + (1.0L/675.0L)*pow(x_1, 6) - 31.0L/2400.0L*pow(x_1, 5) + (1.0L/36.0L)*pow(x_1, 4) + (3.0L/64.0L)*pow(x_1, 3) - 2539.0L/11520.0L*pow(x_1, 2) + (2477.0L/19200.0L)*x_1 + 77293.0L/552960.0L;
       return __pp_r252___result;
    }
    static inline O __pp_r253__(const I &x_0, const I &x_1) {
       O __pp_r253___result;
       __pp_r253___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 + (59.0L/4800.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/480.0L*pow(x_0, 4)*x_1 + (157.0L/2304.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) - 23.0L/288.0L*pow(x_0, 3)*x_1 + (115.0L/1152.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 7.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 15.0L/64.0L*pow(x_0, 2)*x_1 + (2321.0L/46080.0L)*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 1.0L/60.0L*x_0*pow(x_1, 4) + (1.0L/18.0L)*x_0*pow(x_1, 3) - 1.0L/16.0L*x_0*pow(x_1, 2) - 79.0L/11520.0L*x_0*x_1 + (2179.0L/76800.0L)*x_0 + (1.0L/675.0L)*pow(x_1, 6) - 31.0L/2400.0L*pow(x_1, 5) + (1.0L/36.0L)*pow(x_1, 4) + (3.0L/64.0L)*pow(x_1, 3) - 2539.0L/11520.0L*pow(x_1, 2) + (2477.0L/19200.0L)*x_1 + 77437.0L/552960.0L;
       return __pp_r253___result;
    }
    static inline O __pp_r254__(const I &x_0, const I &x_1) {
       O __pp_r254___result;
       __pp_r254___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 + (59.0L/4800.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/480.0L*pow(x_0, 4)*x_1 + (157.0L/2304.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) - 23.0L/288.0L*pow(x_0, 3)*x_1 + (115.0L/1152.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 7.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 15.0L/64.0L*pow(x_0, 2)*x_1 + (2321.0L/46080.0L)*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 1.0L/60.0L*x_0*pow(x_1, 4) + (1.0L/18.0L)*x_0*pow(x_1, 3) - 1.0L/16.0L*x_0*pow(x_1, 2) - 79.0L/11520.0L*x_0*x_1 + (2179.0L/76800.0L)*x_0 + (1.0L/675.0L)*pow(x_1, 6) - 31.0L/2400.0L*pow(x_1, 5) + (1.0L/36.0L)*pow(x_1, 4) + (3.0L/64.0L)*pow(x_1, 3) - 2539.0L/11520.0L*pow(x_1, 2) + (2477.0L/19200.0L)*x_1 + 77437.0L/552960.0L;
       return __pp_r254___result;
    }
    static inline O __pp_r255__(const I &x_0, const I &x_1) {
       O __pp_r255___result;
       __pp_r255___result = (43.0L/14400.0L)*pow(x_0, 6) - 41.0L/3600.0L*pow(x_0, 5)*x_1 + (25.0L/576.0L)*pow(x_0, 5) + (1.0L/80.0L)*pow(x_0, 4)*pow(x_1, 2) - 29.0L/288.0L*pow(x_0, 4)*x_1 + (2353.0L/11520.0L)*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/9.0L)*pow(x_0, 3)*pow(x_1, 2) - 169.0L/480.0L*pow(x_0, 3)*x_1 + (7213.0L/17280.0L)*pow(x_0, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 4) - 7.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (193.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 2047.0L/2880.0L*pow(x_0, 2)*x_1 + (21529.0L/46080.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 1.0L/144.0L*x_0*pow(x_1, 4) - 1.0L/80.0L*x_0*pow(x_1, 3) + (253.0L/1440.0L)*x_0*pow(x_1, 2) - 1627.0L/3840.0L*x_0*x_1 + (14753.0L/46080.0L)*x_0 + (11.0L/7200.0L)*pow(x_1, 6) - 1.0L/72.0L*pow(x_1, 5) + (209.0L/5760.0L)*pow(x_1, 4) + (31.0L/4320.0L)*pow(x_1, 3) - 2677.0L/23040.0L*pow(x_1, 2) - 389.0L/23040.0L*x_1 + 622483.0L/2764800.0L;
       return __pp_r255___result;
    }
    static inline O __pp_r256__(const I &x_0, const I &x_1) {
       O __pp_r256___result;
       __pp_r256___result = (43.0L/14400.0L)*pow(x_0, 6) - 41.0L/3600.0L*pow(x_0, 5)*x_1 + (25.0L/576.0L)*pow(x_0, 5) + (1.0L/80.0L)*pow(x_0, 4)*pow(x_1, 2) - 29.0L/288.0L*pow(x_0, 4)*x_1 + (2353.0L/11520.0L)*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/9.0L)*pow(x_0, 3)*pow(x_1, 2) - 169.0L/480.0L*pow(x_0, 3)*x_1 + (7213.0L/17280.0L)*pow(x_0, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 4) - 7.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (193.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 2047.0L/2880.0L*pow(x_0, 2)*x_1 + (21529.0L/46080.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 1.0L/144.0L*x_0*pow(x_1, 4) - 1.0L/80.0L*x_0*pow(x_1, 3) + (253.0L/1440.0L)*x_0*pow(x_1, 2) - 1627.0L/3840.0L*x_0*x_1 + (14753.0L/46080.0L)*x_0 + (11.0L/7200.0L)*pow(x_1, 6) - 1.0L/72.0L*pow(x_1, 5) + (209.0L/5760.0L)*pow(x_1, 4) + (31.0L/4320.0L)*pow(x_1, 3) - 2677.0L/23040.0L*pow(x_1, 2) - 389.0L/23040.0L*x_1 + 622483.0L/2764800.0L;
       return __pp_r256___result;
    }
    static inline O __pp_r257__(const I &x_0, const I &x_1) {
       O __pp_r257___result;
       __pp_r257___result = (11.0L/4800.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 + (89.0L/2880.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/288.0L*pow(x_0, 4)*x_1 + (1273.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 3)*x_1 + (733.0L/17280.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/36.0L)*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (1193.0L/2880.0L)*pow(x_0, 2)*x_1 - 17351.0L/46080.0L*pow(x_0, 2) + (19.0L/3600.0L)*x_0*pow(x_1, 5) - 5.0L/72.0L*x_0*pow(x_1, 4) + (29.0L/80.0L)*x_0*pow(x_1, 3) - 1367.0L/1440.0L*x_0*pow(x_1, 2) + (4853.0L/3840.0L)*x_0*x_1 - 31903.0L/46080.0L*x_0 + (1.0L/1200.0L)*pow(x_1, 6) - 1.0L/720.0L*pow(x_1, 5) - 331.0L/5760.0L*pow(x_1, 4) + (1651.0L/4320.0L)*pow(x_1, 3) - 22117.0L/23040.0L*pow(x_1, 2) + (22939.0L/23040.0L)*x_1 - 777197.0L/2764800.0L;
       return __pp_r257___result;
    }
    static inline O __pp_r258__(const I &x_0, const I &x_1) {
       O __pp_r258___result;
       __pp_r258___result = (11.0L/4800.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 + (89.0L/2880.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/288.0L*pow(x_0, 4)*x_1 + (1273.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 3)*x_1 + (733.0L/17280.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/36.0L)*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (1193.0L/2880.0L)*pow(x_0, 2)*x_1 - 17351.0L/46080.0L*pow(x_0, 2) + (19.0L/3600.0L)*x_0*pow(x_1, 5) - 5.0L/72.0L*x_0*pow(x_1, 4) + (29.0L/80.0L)*x_0*pow(x_1, 3) - 1367.0L/1440.0L*x_0*pow(x_1, 2) + (4853.0L/3840.0L)*x_0*x_1 - 31903.0L/46080.0L*x_0 + (1.0L/1200.0L)*pow(x_1, 6) - 1.0L/720.0L*pow(x_1, 5) - 331.0L/5760.0L*pow(x_1, 4) + (1651.0L/4320.0L)*pow(x_1, 3) - 22117.0L/23040.0L*pow(x_1, 2) + (22939.0L/23040.0L)*x_1 - 777197.0L/2764800.0L;
       return __pp_r258___result;
    }
    static inline O __pp_r259__(const I &x_0, const I &x_1) {
       O __pp_r259___result;
       __pp_r259___result = -29.0L/43200.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 - 131.0L/14400.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (89.0L/1440.0L)*pow(x_0, 4)*x_1 - 1319.0L/11520.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 41.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (227.0L/480.0L)*pow(x_0, 3)*x_1 - 10931.0L/17280.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (4109.0L/2880.0L)*pow(x_0, 2)*x_1 - 69839.0L/46080.0L*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) - 59.0L/720.0L*x_0*pow(x_1, 4) + (19.0L/40.0L)*x_0*pow(x_1, 3) - 131.0L/90.0L*x_0*pow(x_1, 2) + (9227.0L/3840.0L)*x_0*x_1 - 395711.0L/230400.0L*x_0 + (17.0L/21600.0L)*pow(x_1, 6) - 1.0L/7200.0L*pow(x_1, 5) - 103.0L/1440.0L*pow(x_1, 4) + (4031.0L/8640.0L)*pow(x_1, 3) - 14339.0L/11520.0L*pow(x_1, 2) + (10859.0L/7200.0L)*x_1 - 1840079.0L/2764800.0L;
       return __pp_r259___result;
    }
    static inline O __pp_r260__(const I &x_0, const I &x_1) {
       O __pp_r260___result;
       __pp_r260___result = -29.0L/43200.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 - 131.0L/14400.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (89.0L/1440.0L)*pow(x_0, 4)*x_1 - 1319.0L/11520.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 41.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (227.0L/480.0L)*pow(x_0, 3)*x_1 - 10931.0L/17280.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (4109.0L/2880.0L)*pow(x_0, 2)*x_1 - 69839.0L/46080.0L*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) - 59.0L/720.0L*x_0*pow(x_1, 4) + (19.0L/40.0L)*x_0*pow(x_1, 3) - 131.0L/90.0L*x_0*pow(x_1, 2) + (9227.0L/3840.0L)*x_0*x_1 - 395711.0L/230400.0L*x_0 + (17.0L/21600.0L)*pow(x_1, 6) - 1.0L/7200.0L*pow(x_1, 5) - 103.0L/1440.0L*pow(x_1, 4) + (4031.0L/8640.0L)*pow(x_1, 3) - 14339.0L/11520.0L*pow(x_1, 2) + (10859.0L/7200.0L)*x_1 - 1840079.0L/2764800.0L;
       return __pp_r260___result;
    }
    static inline O __pp_r261__(const I &x_0, const I &x_1) {
       O __pp_r261___result;
       __pp_r261___result = (11.0L/4800.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 + (89.0L/2880.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/288.0L*pow(x_0, 4)*x_1 + (1273.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 3)*x_1 + (733.0L/17280.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/36.0L)*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (1193.0L/2880.0L)*pow(x_0, 2)*x_1 - 17351.0L/46080.0L*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 5.0L/288.0L*x_0*pow(x_1, 4) + (49.0L/480.0L)*x_0*pow(x_1, 3) - 859.0L/2880.0L*x_0*pow(x_1, 2) + (9.0L/20.0L)*x_0*x_1 - 13153.0L/46080.0L*x_0 - 1.0L/4800.0L*pow(x_1, 6) + (11.0L/2880.0L)*pow(x_1, 5) - 287.0L/11520.0L*pow(x_1, 4) + (979.0L/17280.0L)*pow(x_1, 3) + (2641.0L/46080.0L)*pow(x_1, 2) - 19747.0L/46080.0L*x_1 + 666089.0L/1382400.0L;
       return __pp_r261___result;
    }
    static inline O __pp_r262__(const I &x_0, const I &x_1) {
       O __pp_r262___result;
       __pp_r262___result = -29.0L/43200.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 - 131.0L/14400.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (89.0L/1440.0L)*pow(x_0, 4)*x_1 - 1319.0L/11520.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 41.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (227.0L/480.0L)*pow(x_0, 3)*x_1 - 10931.0L/17280.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (4109.0L/2880.0L)*pow(x_0, 2)*x_1 - 69839.0L/46080.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 43.0L/1440.0L*x_0*pow(x_1, 4) + (103.0L/480.0L)*x_0*pow(x_1, 3) - 2317.0L/2880.0L*x_0*pow(x_1, 2) + (1017.0L/640.0L)*x_0*x_1 - 301961.0L/230400.0L*x_0 - 11.0L/43200.0L*pow(x_1, 6) + (73.0L/14400.0L)*pow(x_1, 5) - 449.0L/11520.0L*pow(x_1, 4) + (2437.0L/17280.0L)*pow(x_1, 3) - 10481.0L/46080.0L*pow(x_1, 2) + (19363.0L/230400.0L)*x_1 + 16831.0L/172800.0L;
       return __pp_r262___result;
    }
    static inline O __pp_r263__(const I &x_0, const I &x_1) {
       O __pp_r263___result;
       __pp_r263___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (29.0L/14400.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/160.0L)*pow(x_0, 4)*x_1 - 13.0L/3840.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 - 691.0L/17280.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 337.0L/960.0L*pow(x_0, 2)*x_1 + (4027.0L/15360.0L)*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) + (37.0L/1440.0L)*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) + (2803.0L/2880.0L)*x_0*pow(x_1, 2) - 11327.0L/5760.0L*x_0*x_1 + (353399.0L/230400.0L)*x_0 + (1.0L/4800.0L)*pow(x_1, 6) - 29.0L/4800.0L*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) - 289.0L/640.0L*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) - 211999.0L/76800.0L*x_1 + 38279.0L/19200.0L;
       return __pp_r263___result;
    }
    static inline O __pp_r264__(const I &x_0, const I &x_1) {
       O __pp_r264___result;
       __pp_r264___result = (11.0L/4800.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 + (89.0L/2880.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/288.0L*pow(x_0, 4)*x_1 + (1273.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 3)*x_1 + (733.0L/17280.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/36.0L)*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (1193.0L/2880.0L)*pow(x_0, 2)*x_1 - 17351.0L/46080.0L*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 5.0L/288.0L*x_0*pow(x_1, 4) + (49.0L/480.0L)*x_0*pow(x_1, 3) - 859.0L/2880.0L*x_0*pow(x_1, 2) + (9.0L/20.0L)*x_0*x_1 - 13153.0L/46080.0L*x_0 - 1.0L/4800.0L*pow(x_1, 6) + (11.0L/2880.0L)*pow(x_1, 5) - 287.0L/11520.0L*pow(x_1, 4) + (979.0L/17280.0L)*pow(x_1, 3) + (2641.0L/46080.0L)*pow(x_1, 2) - 19747.0L/46080.0L*x_1 + 666089.0L/1382400.0L;
       return __pp_r264___result;
    }
    static inline O __pp_r265__(const I &x_0, const I &x_1) {
       O __pp_r265___result;
       __pp_r265___result = -29.0L/43200.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 - 131.0L/14400.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (89.0L/1440.0L)*pow(x_0, 4)*x_1 - 1319.0L/11520.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 41.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (227.0L/480.0L)*pow(x_0, 3)*x_1 - 10931.0L/17280.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (4109.0L/2880.0L)*pow(x_0, 2)*x_1 - 69839.0L/46080.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 43.0L/1440.0L*x_0*pow(x_1, 4) + (103.0L/480.0L)*x_0*pow(x_1, 3) - 2317.0L/2880.0L*x_0*pow(x_1, 2) + (1017.0L/640.0L)*x_0*x_1 - 301961.0L/230400.0L*x_0 - 11.0L/43200.0L*pow(x_1, 6) + (73.0L/14400.0L)*pow(x_1, 5) - 449.0L/11520.0L*pow(x_1, 4) + (2437.0L/17280.0L)*pow(x_1, 3) - 10481.0L/46080.0L*pow(x_1, 2) + (19363.0L/230400.0L)*x_1 + 16831.0L/172800.0L;
       return __pp_r265___result;
    }
    static inline O __pp_r266__(const I &x_0, const I &x_1) {
       O __pp_r266___result;
       __pp_r266___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (29.0L/14400.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/160.0L)*pow(x_0, 4)*x_1 - 13.0L/3840.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 - 691.0L/17280.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 337.0L/960.0L*pow(x_0, 2)*x_1 + (4027.0L/15360.0L)*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) + (37.0L/1440.0L)*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) + (2803.0L/2880.0L)*x_0*pow(x_1, 2) - 11327.0L/5760.0L*x_0*x_1 + (353399.0L/230400.0L)*x_0 + (1.0L/4800.0L)*pow(x_1, 6) - 29.0L/4800.0L*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) - 289.0L/640.0L*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) - 211999.0L/76800.0L*x_1 + 38279.0L/19200.0L;
       return __pp_r266___result;
    }
    static inline O __pp_r267__(const I &x_0, const I &x_1) {
       O __pp_r267___result;
       __pp_r267___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (29.0L/14400.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/160.0L)*pow(x_0, 4)*x_1 - 13.0L/3840.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 - 691.0L/17280.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 337.0L/960.0L*pow(x_0, 2)*x_1 + (4027.0L/15360.0L)*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) + (37.0L/1440.0L)*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) + (2803.0L/2880.0L)*x_0*pow(x_1, 2) - 11327.0L/5760.0L*x_0*x_1 + (353399.0L/230400.0L)*x_0 + (1.0L/4800.0L)*pow(x_1, 6) - 29.0L/4800.0L*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) - 289.0L/640.0L*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) - 211999.0L/76800.0L*x_1 + 38279.0L/19200.0L;
       return __pp_r267___result;
    }
    static inline O __pp_r268__(const I &x_0, const I &x_1) {
       O __pp_r268___result;
       __pp_r268___result = -1.0L/800.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 79.0L/3600.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (3.0L/80.0L)*pow(x_0, 4)*x_1 - 71.0L/480.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (11.0L/90.0L)*pow(x_0, 3)*x_1 - 947.0L/2160.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 101.0L/480.0L*pow(x_0, 2)*x_1 - 1187.0L/3840.0L*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) + (37.0L/1440.0L)*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) + (2803.0L/2880.0L)*x_0*pow(x_1, 2) - 21439.0L/11520.0L*x_0*x_1 + (128707.0L/115200.0L)*x_0 + (1.0L/4800.0L)*pow(x_1, 6) - 29.0L/4800.0L*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) - 289.0L/640.0L*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) - 209569.0L/76800.0L*x_1 + 574799.0L/307200.0L;
       return __pp_r268___result;
    }
    static inline O __pp_r269__(const I &x_0, const I &x_1) {
       O __pp_r269___result;
       __pp_r269___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (1.0L/144.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/144.0L*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) - 25.0L/72.0L*pow(x_0, 3)*x_1 + (125.0L/216.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 125.0L/72.0L*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) + (5.0L/144.0L)*x_0*pow(x_1, 4) - 25.0L/72.0L*x_0*pow(x_1, 3) + (125.0L/72.0L)*x_0*pow(x_1, 2) - 625.0L/144.0L*x_0*x_1 + (625.0L/144.0L)*x_0 + (1.0L/4320.0L)*pow(x_1, 6) - 1.0L/144.0L*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) - 125.0L/216.0L*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) - 625.0L/144.0L*x_1 + 3125.0L/864.0L;
       return __pp_r269___result;
    }
    static inline O __pp_r270__(const I &x_0, const I &x_1) {
       O __pp_r270___result;
       __pp_r270___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (29.0L/14400.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/160.0L)*pow(x_0, 4)*x_1 - 13.0L/3840.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 - 691.0L/17280.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 337.0L/960.0L*pow(x_0, 2)*x_1 + (4027.0L/15360.0L)*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) + (37.0L/1440.0L)*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) + (2803.0L/2880.0L)*x_0*pow(x_1, 2) - 11327.0L/5760.0L*x_0*x_1 + (353399.0L/230400.0L)*x_0 + (1.0L/4800.0L)*pow(x_1, 6) - 29.0L/4800.0L*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) - 289.0L/640.0L*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) - 211999.0L/76800.0L*x_1 + 38279.0L/19200.0L;
       return __pp_r270___result;
    }
    static inline O __pp_r271__(const I &x_0, const I &x_1) {
       O __pp_r271___result;
       __pp_r271___result = -1.0L/800.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 79.0L/3600.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (3.0L/80.0L)*pow(x_0, 4)*x_1 - 71.0L/480.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (11.0L/90.0L)*pow(x_0, 3)*x_1 - 947.0L/2160.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 101.0L/480.0L*pow(x_0, 2)*x_1 - 1187.0L/3840.0L*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) + (37.0L/1440.0L)*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) + (2803.0L/2880.0L)*x_0*pow(x_1, 2) - 21439.0L/11520.0L*x_0*x_1 + (128707.0L/115200.0L)*x_0 + (1.0L/4800.0L)*pow(x_1, 6) - 29.0L/4800.0L*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) - 289.0L/640.0L*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) - 209569.0L/76800.0L*x_1 + 574799.0L/307200.0L;
       return __pp_r271___result;
    }
    static inline O __pp_r272__(const I &x_0, const I &x_1) {
       O __pp_r272___result;
       __pp_r272___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (1.0L/144.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/144.0L*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) - 25.0L/72.0L*pow(x_0, 3)*x_1 + (125.0L/216.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 125.0L/72.0L*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) + (5.0L/144.0L)*x_0*pow(x_1, 4) - 25.0L/72.0L*x_0*pow(x_1, 3) + (125.0L/72.0L)*x_0*pow(x_1, 2) - 625.0L/144.0L*x_0*x_1 + (625.0L/144.0L)*x_0 + (1.0L/4320.0L)*pow(x_1, 6) - 1.0L/144.0L*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) - 125.0L/216.0L*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) - 625.0L/144.0L*x_1 + 3125.0L/864.0L;
       return __pp_r272___result;
    }
    static inline O __pp_r273__(const I &x_0, const I &x_1) {
       O __pp_r273___result;
       __pp_r273___result = (1.0L/400.0L)*pow(x_0, 6) - 2.0L/225.0L*pow(x_0, 5)*x_1 + (1.0L/7200.0L)*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/720.0L)*pow(x_0, 4)*x_1 + (1.0L/120.0L)*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/90.0L)*pow(x_0, 3)*x_1 + (1.0L/8640.0L)*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/1440.0L)*pow(x_0, 2)*x_1 - 51.0L/640.0L*pow(x_0, 2) - 1.0L/1800.0L*x_0*pow(x_1, 5) + (1.0L/90.0L)*x_0*pow(x_1, 4) - 13.0L/720.0L*x_0*pow(x_1, 3) + (1.0L/720.0L)*x_0*pow(x_1, 2) - 1.0L/5760.0L*x_0*x_1 + (1.0L/115200.0L)*x_0 - 1.0L/2400.0L*pow(x_1, 6) + (1.0L/225.0L)*pow(x_1, 5) + (11.0L/1920.0L)*pow(x_1, 4) + (1.0L/1080.0L)*pow(x_1, 3) - 613.0L/7680.0L*pow(x_1, 2) + (1.0L/57600.0L)*x_1 + 83291.0L/460800.0L;
       return __pp_r273___result;
    }
    static inline O __pp_r274__(const I &x_0, const I &x_1) {
       O __pp_r274___result;
       __pp_r274___result = (11.0L/4320.0L)*pow(x_0, 6) - 1.0L/120.0L*pow(x_0, 5)*x_1 + (1.0L/144.0L)*pow(x_0, 4)*pow(x_1, 2) + (49.0L/5760.0L)*pow(x_0, 4) + (1.0L/80.0L)*pow(x_0, 3)*x_1 + (1.0L/144.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) - 367.0L/4608.0L*pow(x_0, 2) + (1.0L/120.0L)*x_0*pow(x_1, 5) - 1.0L/80.0L*x_0*pow(x_1, 3) + (11.0L/4320.0L)*pow(x_1, 6) + (49.0L/5760.0L)*pow(x_1, 4) - 367.0L/4608.0L*pow(x_1, 2) + 124937.0L/691200.0L;
       return __pp_r274___result;
    }
    static inline O __pp_r275__(const I &x_0, const I &x_1) {
       O __pp_r275___result;
       __pp_r275___result = -1.0L/2160.0L*pow(x_0, 6) - 31.0L/7200.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/80.0L)*pow(x_0, 4)*x_1 + (1.0L/180.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 3)*x_1 - 7.0L/8640.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/480.0L)*pow(x_0, 2)*x_1 - 23.0L/288.0L*pow(x_0, 2) + (7.0L/720.0L)*x_0*pow(x_1, 4) - 1.0L/60.0L*x_0*pow(x_1, 3) + (1.0L/1440.0L)*x_0*pow(x_1, 2) - 1.0L/115200.0L*x_0 - 1.0L/2160.0L*pow(x_1, 6) + (11.0L/2400.0L)*pow(x_1, 5) + (1.0L/180.0L)*pow(x_1, 4) + (1.0L/960.0L)*pow(x_1, 3) - 23.0L/288.0L*pow(x_1, 2) + (1.0L/38400.0L)*x_1 + 15617.0L/86400.0L;
       return __pp_r275___result;
    }
    static inline O __pp_r276__(const I &x_0, const I &x_1) {
       O __pp_r276___result;
       __pp_r276___result = -1.0L/2400.0L*pow(x_0, 6) + (1.0L/1800.0L)*pow(x_0, 5)*x_1 - 1.0L/225.0L*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/90.0L)*pow(x_0, 4)*x_1 + (11.0L/1920.0L)*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (13.0L/720.0L)*pow(x_0, 3)*x_1 - 1.0L/1080.0L*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/720.0L)*pow(x_0, 2)*x_1 - 613.0L/7680.0L*pow(x_0, 2) + (2.0L/225.0L)*x_0*pow(x_1, 5) - 1.0L/720.0L*x_0*pow(x_1, 4) - 1.0L/90.0L*x_0*pow(x_1, 3) - 1.0L/1440.0L*x_0*pow(x_1, 2) + (1.0L/5760.0L)*x_0*x_1 - 1.0L/57600.0L*x_0 + (1.0L/400.0L)*pow(x_1, 6) + (1.0L/7200.0L)*pow(x_1, 5) + (1.0L/120.0L)*pow(x_1, 4) + (1.0L/8640.0L)*pow(x_1, 3) - 51.0L/640.0L*pow(x_1, 2) + (1.0L/115200.0L)*x_1 + 83291.0L/460800.0L;
       return __pp_r276___result;
    }
    static inline O __pp_r277__(const I &x_0, const I &x_1) {
       O __pp_r277___result;
       __pp_r277___result = (1.0L/400.0L)*pow(x_0, 6) - 2.0L/225.0L*pow(x_0, 5)*x_1 + (1.0L/7200.0L)*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/720.0L)*pow(x_0, 4)*x_1 + (1.0L/120.0L)*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/90.0L)*pow(x_0, 3)*x_1 + (1.0L/8640.0L)*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/1440.0L)*pow(x_0, 2)*x_1 - 51.0L/640.0L*pow(x_0, 2) - 17.0L/3600.0L*x_0*pow(x_1, 5) + (31.0L/1440.0L)*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) + (19.0L/2880.0L)*x_0*pow(x_1, 2) - 17.0L/11520.0L*x_0*x_1 + (1.0L/7200.0L)*x_0 - 7.0L/4800.0L*pow(x_1, 6) + (19.0L/14400.0L)*pow(x_1, 5) + (67.0L/3840.0L)*pow(x_1, 4) - 209.0L/17280.0L*pow(x_1, 3) - 1121.0L/15360.0L*pow(x_1, 2) - 401.0L/230400.0L*x_1 + 166747.0L/921600.0L;
       return __pp_r277___result;
    }
    static inline O __pp_r278__(const I &x_0, const I &x_1) {
       O __pp_r278___result;
       __pp_r278___result = -1.0L/2160.0L*pow(x_0, 6) - 31.0L/7200.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/80.0L)*pow(x_0, 4)*x_1 + (1.0L/180.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 3)*x_1 - 7.0L/8640.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/480.0L)*pow(x_0, 2)*x_1 - 23.0L/288.0L*pow(x_0, 2) - 1.0L/240.0L*x_0*pow(x_1, 5) + (29.0L/1440.0L)*x_0*pow(x_1, 4) - 13.0L/480.0L*x_0*pow(x_1, 3) + (17.0L/2880.0L)*x_0*pow(x_1, 2) - 1.0L/768.0L*x_0*x_1 + (7.0L/57600.0L)*x_0 - 13.0L/8640.0L*pow(x_1, 6) + (7.0L/4800.0L)*pow(x_1, 5) + (199.0L/11520.0L)*pow(x_1, 4) - 23.0L/1920.0L*pow(x_1, 3) - 673.0L/9216.0L*pow(x_1, 2) - 133.0L/76800.0L*x_1 + 500239.0L/2764800.0L;
       return __pp_r278___result;
    }
    static inline O __pp_r279__(const I &x_0, const I &x_1) {
       O __pp_r279___result;
       __pp_r279___result = -1.0L/4320.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 7.0L/2400.0L*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/180.0L)*pow(x_0, 4)*x_1 + (13.0L/1440.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/360.0L)*pow(x_0, 3)*x_1 + (11.0L/2880.0L)*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/240.0L)*pow(x_0, 2)*pow(x_1, 2) - 17.0L/1440.0L*pow(x_0, 2)*x_1 - 11.0L/144.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) + (13.0L/480.0L)*x_0*pow(x_1, 4) - 59.0L/1440.0L*x_0*pow(x_1, 3) + (19.0L/960.0L)*x_0*pow(x_1, 2) - 19.0L/2304.0L*x_0*x_1 + (29.0L/19200.0L)*x_0 - 11.0L/8640.0L*pow(x_1, 6) + (1.0L/14400.0L)*pow(x_1, 5) + (239.0L/11520.0L)*pow(x_1, 4) - 287.0L/17280.0L*pow(x_1, 3) - 641.0L/9216.0L*pow(x_1, 2) - 719.0L/230400.0L*x_1 + 500879.0L/2764800.0L;
       return __pp_r279___result;
    }
    static inline O __pp_r280__(const I &x_0, const I &x_1) {
       O __pp_r280___result;
       __pp_r280___result = (1.0L/400.0L)*pow(x_0, 6) - 2.0L/225.0L*pow(x_0, 5)*x_1 + (1.0L/7200.0L)*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/720.0L)*pow(x_0, 4)*x_1 + (1.0L/120.0L)*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/90.0L)*pow(x_0, 3)*x_1 + (1.0L/8640.0L)*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/1440.0L)*pow(x_0, 2)*x_1 - 51.0L/640.0L*pow(x_0, 2) - 17.0L/3600.0L*x_0*pow(x_1, 5) + (31.0L/1440.0L)*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) + (19.0L/2880.0L)*x_0*pow(x_1, 2) - 17.0L/11520.0L*x_0*x_1 + (1.0L/7200.0L)*x_0 - 7.0L/4800.0L*pow(x_1, 6) + (19.0L/14400.0L)*pow(x_1, 5) + (67.0L/3840.0L)*pow(x_1, 4) - 209.0L/17280.0L*pow(x_1, 3) - 1121.0L/15360.0L*pow(x_1, 2) - 401.0L/230400.0L*x_1 + 166747.0L/921600.0L;
       return __pp_r280___result;
    }
    static inline O __pp_r281__(const I &x_0, const I &x_1) {
       O __pp_r281___result;
       __pp_r281___result = -1.0L/2160.0L*pow(x_0, 6) - 31.0L/7200.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/80.0L)*pow(x_0, 4)*x_1 + (1.0L/180.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 3)*x_1 - 7.0L/8640.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/480.0L)*pow(x_0, 2)*x_1 - 23.0L/288.0L*pow(x_0, 2) - 1.0L/240.0L*x_0*pow(x_1, 5) + (29.0L/1440.0L)*x_0*pow(x_1, 4) - 13.0L/480.0L*x_0*pow(x_1, 3) + (17.0L/2880.0L)*x_0*pow(x_1, 2) - 1.0L/768.0L*x_0*x_1 + (7.0L/57600.0L)*x_0 - 13.0L/8640.0L*pow(x_1, 6) + (7.0L/4800.0L)*pow(x_1, 5) + (199.0L/11520.0L)*pow(x_1, 4) - 23.0L/1920.0L*pow(x_1, 3) - 673.0L/9216.0L*pow(x_1, 2) - 133.0L/76800.0L*x_1 + 500239.0L/2764800.0L;
       return __pp_r281___result;
    }
    static inline O __pp_r282__(const I &x_0, const I &x_1) {
       O __pp_r282___result;
       __pp_r282___result = -1.0L/4320.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 7.0L/2400.0L*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/180.0L)*pow(x_0, 4)*x_1 + (13.0L/1440.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/360.0L)*pow(x_0, 3)*x_1 + (11.0L/2880.0L)*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/240.0L)*pow(x_0, 2)*pow(x_1, 2) - 17.0L/1440.0L*pow(x_0, 2)*x_1 - 11.0L/144.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) + (13.0L/480.0L)*x_0*pow(x_1, 4) - 59.0L/1440.0L*x_0*pow(x_1, 3) + (19.0L/960.0L)*x_0*pow(x_1, 2) - 19.0L/2304.0L*x_0*x_1 + (29.0L/19200.0L)*x_0 - 11.0L/8640.0L*pow(x_1, 6) + (1.0L/14400.0L)*pow(x_1, 5) + (239.0L/11520.0L)*pow(x_1, 4) - 287.0L/17280.0L*pow(x_1, 3) - 641.0L/9216.0L*pow(x_1, 2) - 719.0L/230400.0L*x_1 + 500879.0L/2764800.0L;
       return __pp_r282___result;
    }
    static inline O __pp_r283__(const I &x_0, const I &x_1) {
       O __pp_r283___result;
       __pp_r283___result = -37.0L/21600.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 23.0L/2400.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/45.0L)*pow(x_0, 4)*x_1 - 1.0L/288.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/36.0L)*pow(x_0, 3)*x_1 - 5.0L/576.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/144.0L)*pow(x_0, 2)*x_1 - 961.0L/11520.0L*pow(x_0, 2) - 19.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/40.0L)*x_0*pow(x_1, 4) - 5.0L/144.0L*x_0*pow(x_1, 3) + (1.0L/96.0L)*x_0*pow(x_1, 2) - 7.0L/5760.0L*x_0*x_1 - 23.0L/38400.0L*x_0 - 7.0L/5400.0L*pow(x_1, 6) + (1.0L/3600.0L)*pow(x_1, 5) + (23.0L/1152.0L)*pow(x_1, 4) - 13.0L/864.0L*pow(x_1, 3) - 1643.0L/23040.0L*pow(x_1, 2) - 119.0L/57600.0L*x_1 + 10003.0L/55296.0L;
       return __pp_r283___result;
    }
    static inline O __pp_r284__(const I &x_0, const I &x_1) {
       O __pp_r284___result;
       __pp_r284___result = -37.0L/21600.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/800.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/45.0L)*pow(x_0, 4)*x_1 + (5.0L/288.0L)*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/36.0L)*pow(x_0, 3)*x_1 + (7.0L/576.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/144.0L)*pow(x_0, 2)*x_1 - 841.0L/11520.0L*pow(x_0, 2) - 19.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/40.0L)*x_0*pow(x_1, 4) - 5.0L/144.0L*x_0*pow(x_1, 3) + (1.0L/96.0L)*x_0*pow(x_1, 2) - 7.0L/5760.0L*x_0*x_1 + (77.0L/38400.0L)*x_0 - 7.0L/5400.0L*pow(x_1, 6) + (1.0L/3600.0L)*pow(x_1, 5) + (23.0L/1152.0L)*pow(x_1, 4) - 13.0L/864.0L*pow(x_1, 3) - 1643.0L/23040.0L*pow(x_1, 2) - 119.0L/57600.0L*x_1 + 50087.0L/276480.0L;
       return __pp_r284___result;
    }
    static inline O __pp_r285__(const I &x_0, const I &x_1) {
       O __pp_r285___result;
       __pp_r285___result = -1.0L/4320.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (71.0L/7200.0L)*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/180.0L*pow(x_0, 4)*x_1 + (5.0L/96.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (7.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 1.0L/24.0L*pow(x_0, 3)*x_1 + (121.0L/1728.0L)*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 23.0L/288.0L*pow(x_0, 2)*x_1 - 3.0L/160.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) + (41.0L/1440.0L)*x_0*pow(x_1, 4) - 5.0L/96.0L*x_0*pow(x_1, 3) + (31.0L/576.0L)*x_0*pow(x_1, 2) - 71.0L/1280.0L*x_0*x_1 + (839.0L/28800.0L)*x_0 - 11.0L/8640.0L*pow(x_1, 6) - 1.0L/14400.0L*pow(x_1, 5) + (17.0L/768.0L)*pow(x_1, 4) - 77.0L/3456.0L*pow(x_1, 3) - 887.0L/15360.0L*pow(x_1, 2) - 3601.0L/230400.0L*x_1 + 34433.0L/184320.0L;
       return __pp_r285___result;
    }
    static inline O __pp_r286__(const I &x_0, const I &x_1) {
       O __pp_r286___result;
       __pp_r286___result = -1.0L/360.0L*pow(x_0, 5)*x_1 + (91.0L/7200.0L)*pow(x_0, 5) - 7.0L/360.0L*pow(x_0, 4)*x_1 + (19.0L/288.0L)*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) + (17.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 7.0L/72.0L*pow(x_0, 3)*x_1 + (185.0L/1728.0L)*pow(x_0, 3) - 11.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 55.0L/288.0L*pow(x_0, 2)*x_1 + (53.0L/1440.0L)*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) + (61.0L/1440.0L)*x_0*pow(x_1, 4) - 31.0L/288.0L*x_0*pow(x_1, 3) + (95.0L/576.0L)*x_0*pow(x_1, 2) - 1919.0L/11520.0L*x_0*x_1 + (2119.0L/28800.0L)*x_0 - 1.0L/960.0L*pow(x_1, 6) - 41.0L/14400.0L*pow(x_1, 5) + (83.0L/2304.0L)*pow(x_1, 4) - 205.0L/3456.0L*pow(x_1, 3) - 101.0L/46080.0L*pow(x_1, 2) - 13841.0L/230400.0L*x_1 + 111491.0L/552960.0L;
       return __pp_r286___result;
    }
    static inline O __pp_r287__(const I &x_0, const I &x_1) {
       O __pp_r287___result;
       __pp_r287___result = (2.0L/675.0L)*pow(x_0, 6) - 7.0L/600.0L*pow(x_0, 5)*x_1 + (7.0L/160.0L)*pow(x_0, 5) + (1.0L/90.0L)*pow(x_0, 4)*pow(x_1, 2) - 7.0L/72.0L*pow(x_0, 4)*x_1 + (97.0L/480.0L)*pow(x_0, 4) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/8.0L)*pow(x_0, 3)*pow(x_1, 2) - 133.0L/360.0L*pow(x_0, 3)*x_1 + (1223.0L/2880.0L)*pow(x_0, 3) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) - 961.0L/1440.0L*pow(x_0, 2)*x_1 + (871.0L/1920.0L)*pow(x_0, 2) - 3.0L/400.0L*x_0*pow(x_1, 5) + (5.0L/96.0L)*x_0*pow(x_1, 4) - 253.0L/1440.0L*x_0*pow(x_1, 3) + (129.0L/320.0L)*x_0*pow(x_1, 2) - 6721.0L/11520.0L*x_0*x_1 + (1403.0L/3840.0L)*x_0 - 43.0L/43200.0L*pow(x_1, 6) - 11.0L/2880.0L*pow(x_1, 5) + (57.0L/1280.0L)*pow(x_1, 4) - 1711.0L/17280.0L*pow(x_1, 3) + (1567.0L/15360.0L)*pow(x_1, 2) - 9491.0L/46080.0L*x_1 + 264251.0L/921600.0L;
       return __pp_r287___result;
    }
    static inline O __pp_r288__(const I &x_0, const I &x_1) {
       O __pp_r288___result;
       __pp_r288___result = -37.0L/21600.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 23.0L/2400.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/45.0L)*pow(x_0, 4)*x_1 - 1.0L/288.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/36.0L)*pow(x_0, 3)*x_1 - 5.0L/576.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/144.0L)*pow(x_0, 2)*x_1 - 961.0L/11520.0L*pow(x_0, 2) - 19.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/40.0L)*x_0*pow(x_1, 4) - 5.0L/144.0L*x_0*pow(x_1, 3) + (1.0L/96.0L)*x_0*pow(x_1, 2) - 7.0L/5760.0L*x_0*x_1 - 23.0L/38400.0L*x_0 - 7.0L/5400.0L*pow(x_1, 6) + (1.0L/3600.0L)*pow(x_1, 5) + (23.0L/1152.0L)*pow(x_1, 4) - 13.0L/864.0L*pow(x_1, 3) - 1643.0L/23040.0L*pow(x_1, 2) - 119.0L/57600.0L*x_1 + 10003.0L/55296.0L;
       return __pp_r288___result;
    }
    static inline O __pp_r289__(const I &x_0, const I &x_1) {
       O __pp_r289___result;
       __pp_r289___result = -37.0L/21600.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/800.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/45.0L)*pow(x_0, 4)*x_1 + (5.0L/288.0L)*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/36.0L)*pow(x_0, 3)*x_1 + (7.0L/576.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/144.0L)*pow(x_0, 2)*x_1 - 841.0L/11520.0L*pow(x_0, 2) - 19.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/40.0L)*x_0*pow(x_1, 4) - 5.0L/144.0L*x_0*pow(x_1, 3) + (1.0L/96.0L)*x_0*pow(x_1, 2) - 7.0L/5760.0L*x_0*x_1 + (77.0L/38400.0L)*x_0 - 7.0L/5400.0L*pow(x_1, 6) + (1.0L/3600.0L)*pow(x_1, 5) + (23.0L/1152.0L)*pow(x_1, 4) - 13.0L/864.0L*pow(x_1, 3) - 1643.0L/23040.0L*pow(x_1, 2) - 119.0L/57600.0L*x_1 + 50087.0L/276480.0L;
       return __pp_r289___result;
    }
    static inline O __pp_r290__(const I &x_0, const I &x_1) {
       O __pp_r290___result;
       __pp_r290___result = -1.0L/4320.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (71.0L/7200.0L)*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/180.0L*pow(x_0, 4)*x_1 + (5.0L/96.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (7.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 1.0L/24.0L*pow(x_0, 3)*x_1 + (121.0L/1728.0L)*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 23.0L/288.0L*pow(x_0, 2)*x_1 - 3.0L/160.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) + (41.0L/1440.0L)*x_0*pow(x_1, 4) - 5.0L/96.0L*x_0*pow(x_1, 3) + (31.0L/576.0L)*x_0*pow(x_1, 2) - 71.0L/1280.0L*x_0*x_1 + (839.0L/28800.0L)*x_0 - 11.0L/8640.0L*pow(x_1, 6) - 1.0L/14400.0L*pow(x_1, 5) + (17.0L/768.0L)*pow(x_1, 4) - 77.0L/3456.0L*pow(x_1, 3) - 887.0L/15360.0L*pow(x_1, 2) - 3601.0L/230400.0L*x_1 + 34433.0L/184320.0L;
       return __pp_r290___result;
    }
    static inline O __pp_r291__(const I &x_0, const I &x_1) {
       O __pp_r291___result;
       __pp_r291___result = -1.0L/360.0L*pow(x_0, 5)*x_1 + (91.0L/7200.0L)*pow(x_0, 5) - 7.0L/360.0L*pow(x_0, 4)*x_1 + (19.0L/288.0L)*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) + (17.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 7.0L/72.0L*pow(x_0, 3)*x_1 + (185.0L/1728.0L)*pow(x_0, 3) - 11.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 55.0L/288.0L*pow(x_0, 2)*x_1 + (53.0L/1440.0L)*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) + (61.0L/1440.0L)*x_0*pow(x_1, 4) - 31.0L/288.0L*x_0*pow(x_1, 3) + (95.0L/576.0L)*x_0*pow(x_1, 2) - 1919.0L/11520.0L*x_0*x_1 + (2119.0L/28800.0L)*x_0 - 1.0L/960.0L*pow(x_1, 6) - 41.0L/14400.0L*pow(x_1, 5) + (83.0L/2304.0L)*pow(x_1, 4) - 205.0L/3456.0L*pow(x_1, 3) - 101.0L/46080.0L*pow(x_1, 2) - 13841.0L/230400.0L*x_1 + 111491.0L/552960.0L;
       return __pp_r291___result;
    }
    static inline O __pp_r292__(const I &x_0, const I &x_1) {
       O __pp_r292___result;
       __pp_r292___result = (2.0L/675.0L)*pow(x_0, 6) - 7.0L/600.0L*pow(x_0, 5)*x_1 + (7.0L/160.0L)*pow(x_0, 5) + (1.0L/90.0L)*pow(x_0, 4)*pow(x_1, 2) - 7.0L/72.0L*pow(x_0, 4)*x_1 + (97.0L/480.0L)*pow(x_0, 4) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/8.0L)*pow(x_0, 3)*pow(x_1, 2) - 133.0L/360.0L*pow(x_0, 3)*x_1 + (1223.0L/2880.0L)*pow(x_0, 3) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) - 961.0L/1440.0L*pow(x_0, 2)*x_1 + (871.0L/1920.0L)*pow(x_0, 2) - 3.0L/400.0L*x_0*pow(x_1, 5) + (5.0L/96.0L)*x_0*pow(x_1, 4) - 253.0L/1440.0L*x_0*pow(x_1, 3) + (129.0L/320.0L)*x_0*pow(x_1, 2) - 6721.0L/11520.0L*x_0*x_1 + (1403.0L/3840.0L)*x_0 - 43.0L/43200.0L*pow(x_1, 6) - 11.0L/2880.0L*pow(x_1, 5) + (57.0L/1280.0L)*pow(x_1, 4) - 1711.0L/17280.0L*pow(x_1, 3) + (1567.0L/15360.0L)*pow(x_1, 2) - 9491.0L/46080.0L*x_1 + 264251.0L/921600.0L;
       return __pp_r292___result;
    }
    static inline O __pp_r293__(const I &x_0, const I &x_1) {
       O __pp_r293___result;
       __pp_r293___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 + (59.0L/4800.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/480.0L*pow(x_0, 4)*x_1 + (157.0L/2304.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) - 23.0L/288.0L*pow(x_0, 3)*x_1 + (115.0L/1152.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 7.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 15.0L/64.0L*pow(x_0, 2)*x_1 + (2321.0L/46080.0L)*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 1.0L/60.0L*x_0*pow(x_1, 4) + (1.0L/18.0L)*x_0*pow(x_1, 3) - 1.0L/16.0L*x_0*pow(x_1, 2) - 79.0L/11520.0L*x_0*x_1 + (2179.0L/76800.0L)*x_0 + (1.0L/675.0L)*pow(x_1, 6) - 31.0L/2400.0L*pow(x_1, 5) + (1.0L/36.0L)*pow(x_1, 4) + (3.0L/64.0L)*pow(x_1, 3) - 2539.0L/11520.0L*pow(x_1, 2) + (2477.0L/19200.0L)*x_1 + 77437.0L/552960.0L;
       return __pp_r293___result;
    }
    static inline O __pp_r294__(const I &x_0, const I &x_1) {
       O __pp_r294___result;
       __pp_r294___result = (43.0L/14400.0L)*pow(x_0, 6) - 41.0L/3600.0L*pow(x_0, 5)*x_1 + (25.0L/576.0L)*pow(x_0, 5) + (1.0L/80.0L)*pow(x_0, 4)*pow(x_1, 2) - 29.0L/288.0L*pow(x_0, 4)*x_1 + (2353.0L/11520.0L)*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/9.0L)*pow(x_0, 3)*pow(x_1, 2) - 169.0L/480.0L*pow(x_0, 3)*x_1 + (7213.0L/17280.0L)*pow(x_0, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 4) - 7.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (193.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 2047.0L/2880.0L*pow(x_0, 2)*x_1 + (21529.0L/46080.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 1.0L/144.0L*x_0*pow(x_1, 4) - 1.0L/80.0L*x_0*pow(x_1, 3) + (253.0L/1440.0L)*x_0*pow(x_1, 2) - 1627.0L/3840.0L*x_0*x_1 + (14753.0L/46080.0L)*x_0 + (11.0L/7200.0L)*pow(x_1, 6) - 1.0L/72.0L*pow(x_1, 5) + (209.0L/5760.0L)*pow(x_1, 4) + (31.0L/4320.0L)*pow(x_1, 3) - 2677.0L/23040.0L*pow(x_1, 2) - 389.0L/23040.0L*x_1 + 622483.0L/2764800.0L;
       return __pp_r294___result;
    }
    static inline O __pp_r295__(const I &x_0, const I &x_1) {
       O __pp_r295___result;
       __pp_r295___result = (11.0L/4800.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 + (89.0L/2880.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 11.0L/288.0L*pow(x_0, 4)*x_1 + (1273.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 3)*x_1 + (733.0L/17280.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/36.0L)*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (1193.0L/2880.0L)*pow(x_0, 2)*x_1 - 17351.0L/46080.0L*pow(x_0, 2) + (19.0L/3600.0L)*x_0*pow(x_1, 5) - 5.0L/72.0L*x_0*pow(x_1, 4) + (29.0L/80.0L)*x_0*pow(x_1, 3) - 1367.0L/1440.0L*x_0*pow(x_1, 2) + (4853.0L/3840.0L)*x_0*x_1 - 31903.0L/46080.0L*x_0 + (1.0L/1200.0L)*pow(x_1, 6) - 1.0L/720.0L*pow(x_1, 5) - 331.0L/5760.0L*pow(x_1, 4) + (1651.0L/4320.0L)*pow(x_1, 3) - 22117.0L/23040.0L*pow(x_1, 2) + (22939.0L/23040.0L)*x_1 - 777197.0L/2764800.0L;
       return __pp_r295___result;
    }
    static inline O __pp_r296__(const I &x_0, const I &x_1) {
       O __pp_r296___result;
       __pp_r296___result = -29.0L/43200.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 - 131.0L/14400.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (89.0L/1440.0L)*pow(x_0, 4)*x_1 - 1319.0L/11520.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 41.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (227.0L/480.0L)*pow(x_0, 3)*x_1 - 10931.0L/17280.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (4109.0L/2880.0L)*pow(x_0, 2)*x_1 - 69839.0L/46080.0L*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) - 59.0L/720.0L*x_0*pow(x_1, 4) + (19.0L/40.0L)*x_0*pow(x_1, 3) - 131.0L/90.0L*x_0*pow(x_1, 2) + (9227.0L/3840.0L)*x_0*x_1 - 395711.0L/230400.0L*x_0 + (17.0L/21600.0L)*pow(x_1, 6) - 1.0L/7200.0L*pow(x_1, 5) - 103.0L/1440.0L*pow(x_1, 4) + (4031.0L/8640.0L)*pow(x_1, 3) - 14339.0L/11520.0L*pow(x_1, 2) + (10859.0L/7200.0L)*x_1 - 1840079.0L/2764800.0L;
       return __pp_r296___result;
    }
    static inline O __pp_r297__(const I &x_0, const I &x_1) {
       O __pp_r297___result;
       __pp_r297___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 - 119.0L/3600.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (67.0L/720.0L)*pow(x_0, 4)*x_1 - 373.0L/1440.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 41.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (17.0L/30.0L)*pow(x_0, 3)*x_1 - 2227.0L/2160.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (2257.0L/1440.0L)*pow(x_0, 2)*x_1 - 24041.0L/11520.0L*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) - 59.0L/720.0L*x_0*pow(x_1, 4) + (19.0L/40.0L)*x_0*pow(x_1, 3) - 131.0L/90.0L*x_0*pow(x_1, 2) + (301.0L/120.0L)*x_0*x_1 - 30731.0L/14400.0L*x_0 + (17.0L/21600.0L)*pow(x_1, 6) - 1.0L/7200.0L*pow(x_1, 5) - 103.0L/1440.0L*pow(x_1, 4) + (4031.0L/8640.0L)*pow(x_1, 3) - 14339.0L/11520.0L*pow(x_1, 2) + (177389.0L/115200.0L)*x_1 - 272383.0L/345600.0L;
       return __pp_r297___result;
    }
    static inline O __pp_r298__(const I &x_0, const I &x_1) {
       O __pp_r298___result;
       __pp_r298___result = -1.0L/360.0L*pow(x_0, 5)*x_1 + (91.0L/7200.0L)*pow(x_0, 5) - 7.0L/360.0L*pow(x_0, 4)*x_1 + (19.0L/288.0L)*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) + (17.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 7.0L/72.0L*pow(x_0, 3)*x_1 + (185.0L/1728.0L)*pow(x_0, 3) - 11.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 55.0L/288.0L*pow(x_0, 2)*x_1 + (53.0L/1440.0L)*pow(x_0, 2) - 1.0L/360.0L*x_0*pow(x_1, 5) + (1.0L/90.0L)*x_0*pow(x_1, 4) - 1.0L/72.0L*x_0*pow(x_1, 3) + (7.0L/288.0L)*x_0*pow(x_1, 2) - 11.0L/180.0L*x_0*x_1 + (4831.0L/115200.0L)*x_0 - 13.0L/7200.0L*pow(x_1, 5) - 1.0L/144.0L*pow(x_1, 4) + (181.0L/1728.0L)*pow(x_1, 3) - 791.0L/2880.0L*pow(x_1, 2) + (17987.0L/115200.0L)*x_1 + 9289.0L/69120.0L;
       return __pp_r298___result;
    }
    static inline O __pp_r299__(const I &x_0, const I &x_1) {
       O __pp_r299___result;
       __pp_r299___result = (2.0L/675.0L)*pow(x_0, 6) - 7.0L/600.0L*pow(x_0, 5)*x_1 + (7.0L/160.0L)*pow(x_0, 5) + (1.0L/90.0L)*pow(x_0, 4)*pow(x_1, 2) - 7.0L/72.0L*pow(x_0, 4)*x_1 + (97.0L/480.0L)*pow(x_0, 4) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/8.0L)*pow(x_0, 3)*pow(x_1, 2) - 133.0L/360.0L*pow(x_0, 3)*x_1 + (1223.0L/2880.0L)*pow(x_0, 3) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) - 961.0L/1440.0L*pow(x_0, 2)*x_1 + (871.0L/1920.0L)*pow(x_0, 2) - 1.0L/300.0L*x_0*pow(x_1, 5) + (1.0L/48.0L)*x_0*pow(x_1, 4) - 59.0L/720.0L*x_0*pow(x_1, 3) + (21.0L/80.0L)*x_0*pow(x_1, 2) - 2753.0L/5760.0L*x_0*x_1 + (2563.0L/7680.0L)*x_0 + (1.0L/21600.0L)*pow(x_1, 6) - 1.0L/360.0L*pow(x_1, 5) + (1.0L/640.0L)*pow(x_1, 4) + (281.0L/4320.0L)*pow(x_1, 3) - 1309.0L/7680.0L*pow(x_1, 2) + (59.0L/5760.0L)*x_1 + 101143.0L/460800.0L;
       return __pp_r299___result;
    }
    static inline O __pp_r300__(const I &x_0, const I &x_1) {
       O __pp_r300___result;
       __pp_r300___result = (49.0L/21600.0L)*pow(x_0, 6) - 3.0L/400.0L*pow(x_0, 5)*x_1 + (1.0L/32.0L)*pow(x_0, 5) + (1.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/144.0L*pow(x_0, 4)*x_1 + (13.0L/120.0L)*pow(x_0, 4) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/180.0L)*pow(x_0, 3)*x_1 + (143.0L/2880.0L)*pow(x_0, 3) - 11.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/18.0L)*pow(x_0, 2)*pow(x_1, 3) - 17.0L/80.0L*pow(x_0, 2)*pow(x_1, 2) + (659.0L/1440.0L)*pow(x_0, 2)*x_1 - 749.0L/1920.0L*pow(x_0, 2) + (1.0L/1200.0L)*x_0*pow(x_1, 5) - 1.0L/24.0L*x_0*pow(x_1, 4) + (211.0L/720.0L)*x_0*pow(x_1, 3) - 69.0L/80.0L*x_0*pow(x_1, 2) + (6967.0L/5760.0L)*x_0*x_1 - 5213.0L/7680.0L*x_0 - 7.0L/10800.0L*pow(x_1, 6) + (7.0L/720.0L)*pow(x_1, 5) - 59.0L/640.0L*pow(x_1, 4) + (1901.0L/4320.0L)*pow(x_1, 3) - 7789.0L/7680.0L*pow(x_1, 2) + (5891.0L/5760.0L)*x_1 - 132137.0L/460800.0L;
       return __pp_r300___result;
    }
    static inline O __pp_r301__(const I &x_0, const I &x_1) {
       O __pp_r301___result;
       __pp_r301___result = -1.0L/1440.0L*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 - 7.0L/800.0L*pow(x_0, 5) - 1.0L/96.0L*pow(x_0, 4)*pow(x_1, 2) + (47.0L/720.0L)*pow(x_0, 4)*x_1 - 7.0L/60.0L*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/10.0L*pow(x_0, 3)*pow(x_1, 2) + (41.0L/90.0L)*pow(x_0, 3)*x_1 - 1801.0L/2880.0L*pow(x_0, 3) - 1.0L/96.0L*pow(x_0, 2)*pow(x_1, 4) + (19.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (2117.0L/1440.0L)*pow(x_0, 2)*x_1 - 367.0L/240.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 13.0L/240.0L*x_0*pow(x_1, 4) + (73.0L/180.0L)*x_0*pow(x_1, 3) - 219.0L/160.0L*x_0*pow(x_1, 2) + (1691.0L/720.0L)*x_0*x_1 - 65431.0L/38400.0L*x_0 - 1.0L/1440.0L*pow(x_1, 6) + (79.0L/7200.0L)*pow(x_1, 5) - 17.0L/160.0L*pow(x_1, 4) + (4531.0L/8640.0L)*pow(x_1, 3) - 1247.0L/960.0L*pow(x_1, 2) + (176869.0L/115200.0L)*x_1 - 77321.0L/115200.0L;
       return __pp_r301___result;
    }
    static inline O __pp_r302__(const I &x_0, const I &x_1) {
       O __pp_r302___result;
       __pp_r302___result = -1.0L/576.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 - 157.0L/4800.0L*pow(x_0, 5) - 1.0L/96.0L*pow(x_0, 4)*pow(x_1, 2) + (139.0L/1440.0L)*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/10.0L*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 - 5897.0L/5760.0L*pow(x_0, 3) - 1.0L/96.0L*pow(x_0, 2)*pow(x_1, 4) + (19.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (4639.0L/2880.0L)*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 13.0L/240.0L*x_0*pow(x_1, 4) + (73.0L/180.0L)*x_0*pow(x_1, 3) - 219.0L/160.0L*x_0*pow(x_1, 2) + (28271.0L/11520.0L)*x_0*x_1 - 162857.0L/76800.0L*x_0 - 1.0L/1440.0L*pow(x_1, 6) + (79.0L/7200.0L)*pow(x_1, 5) - 17.0L/160.0L*pow(x_1, 4) + (4531.0L/8640.0L)*pow(x_1, 3) - 1247.0L/960.0L*pow(x_1, 2) + (90257.0L/57600.0L)*x_1 - 731563.0L/921600.0L;
       return __pp_r302___result;
    }
    static inline O __pp_r303__(const I &x_0, const I &x_1) {
       O __pp_r303___result;
       __pp_r303___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 - 119.0L/3600.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (67.0L/720.0L)*pow(x_0, 4)*x_1 - 373.0L/1440.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 41.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (17.0L/30.0L)*pow(x_0, 3)*x_1 - 2227.0L/2160.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (2257.0L/1440.0L)*pow(x_0, 2)*x_1 - 24041.0L/11520.0L*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) - 59.0L/720.0L*x_0*pow(x_1, 4) + (19.0L/40.0L)*x_0*pow(x_1, 3) - 131.0L/90.0L*x_0*pow(x_1, 2) + (301.0L/120.0L)*x_0*x_1 - 30731.0L/14400.0L*x_0 + (17.0L/21600.0L)*pow(x_1, 6) - 1.0L/7200.0L*pow(x_1, 5) - 103.0L/1440.0L*pow(x_1, 4) + (4031.0L/8640.0L)*pow(x_1, 3) - 14339.0L/11520.0L*pow(x_1, 2) + (177389.0L/115200.0L)*x_1 - 272383.0L/345600.0L;
       return __pp_r303___result;
    }
    static inline O __pp_r304__(const I &x_0, const I &x_1) {
       O __pp_r304___result;
       __pp_r304___result = -1.0L/800.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 79.0L/3600.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (3.0L/80.0L)*pow(x_0, 4)*x_1 - 71.0L/480.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (11.0L/90.0L)*pow(x_0, 3)*x_1 - 947.0L/2160.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 101.0L/480.0L*pow(x_0, 2)*x_1 - 1187.0L/3840.0L*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) - 19.0L/720.0L*x_0*pow(x_1, 4) + (11.0L/360.0L)*x_0*pow(x_1, 3) + (29.0L/90.0L)*x_0*pow(x_1, 2) - 377.0L/360.0L*x_0*x_1 + (10229.0L/14400.0L)*x_0 + (1.0L/800.0L)*pow(x_1, 6) - 9.0L/800.0L*pow(x_1, 5) + (19.0L/480.0L)*pow(x_1, 4) - 121.0L/960.0L*pow(x_1, 3) + (2047.0L/3840.0L)*pow(x_1, 2) - 16699.0L/12800.0L*x_1 + 42553.0L/38400.0L;
       return __pp_r304___result;
    }
    static inline O __pp_r305__(const I &x_0, const I &x_1) {
       O __pp_r305___result;
       __pp_r305___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (1.0L/144.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/144.0L*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) - 25.0L/72.0L*pow(x_0, 3)*x_1 + (125.0L/216.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 125.0L/72.0L*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) - 5.0L/288.0L*x_0*pow(x_1, 4) - 25.0L/288.0L*x_0*pow(x_1, 3) + (625.0L/576.0L)*x_0*pow(x_1, 2) - 8125.0L/2304.0L*x_0*x_1 + (18125.0L/4608.0L)*x_0 + (11.0L/8640.0L)*pow(x_1, 6) - 7.0L/576.0L*pow(x_1, 5) + (125.0L/2304.0L)*pow(x_1, 4) - 875.0L/3456.0L*pow(x_1, 3) + (10625.0L/9216.0L)*pow(x_1, 2) - 26875.0L/9216.0L*x_1 + 315625.0L/110592.0L;
       return __pp_r305___result;
    }
    static inline O __pp_r306__(const I &x_0, const I &x_1) {
       O __pp_r306___result;
       __pp_r306___result = -1.0L/576.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 - 157.0L/4800.0L*pow(x_0, 5) - 1.0L/96.0L*pow(x_0, 4)*pow(x_1, 2) + (139.0L/1440.0L)*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/10.0L*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 - 5897.0L/5760.0L*pow(x_0, 3) - 1.0L/96.0L*pow(x_0, 2)*pow(x_1, 4) + (19.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (4639.0L/2880.0L)*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 13.0L/240.0L*x_0*pow(x_1, 4) + (73.0L/180.0L)*x_0*pow(x_1, 3) - 219.0L/160.0L*x_0*pow(x_1, 2) + (28271.0L/11520.0L)*x_0*x_1 - 162857.0L/76800.0L*x_0 - 1.0L/1440.0L*pow(x_1, 6) + (79.0L/7200.0L)*pow(x_1, 5) - 17.0L/160.0L*pow(x_1, 4) + (4531.0L/8640.0L)*pow(x_1, 3) - 1247.0L/960.0L*pow(x_1, 2) + (90257.0L/57600.0L)*x_1 - 731563.0L/921600.0L;
       return __pp_r306___result;
    }
    static inline O __pp_r307__(const I &x_0, const I &x_1) {
       O __pp_r307___result;
       __pp_r307___result = -11.0L/8640.0L*pow(x_0, 6) + (1.0L/360.0L)*pow(x_0, 5)*x_1 - 311.0L/14400.0L*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) + (59.0L/1440.0L)*pow(x_0, 4)*x_1 - 1729.0L/11520.0L*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (151.0L/1440.0L)*pow(x_0, 3)*x_1 - 7451.0L/17280.0L*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) - 481.0L/2880.0L*pow(x_0, 2)*x_1 - 14869.0L/46080.0L*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) + (1.0L/720.0L)*x_0*pow(x_1, 4) - 7.0L/180.0L*x_0*pow(x_1, 3) + (589.0L/1440.0L)*x_0*pow(x_1, 2) - 12689.0L/11520.0L*x_0*x_1 + (166789.0L/230400.0L)*x_0 - 1.0L/4320.0L*pow(x_1, 6) - 1.0L/7200.0L*pow(x_1, 5) + (7.0L/1440.0L)*pow(x_1, 4) - 589.0L/8640.0L*pow(x_1, 3) + (1379.0L/2880.0L)*pow(x_1, 2) - 73583.0L/57600.0L*x_1 + 3048191.0L/2764800.0L;
       return __pp_r307___result;
    }
    static inline O __pp_r308__(const I &x_0, const I &x_1) {
       O __pp_r308___result;
       __pp_r308___result = (1.0L/4800.0L)*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (7.0L/960.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/32.0L*pow(x_0, 4)*x_1 + (65.0L/768.0L)*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/12.0L)*pow(x_0, 3)*pow(x_1, 2) - 35.0L/96.0L*pow(x_0, 3)*x_1 + (75.0L/128.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 325.0L/192.0L*pow(x_0, 2)*x_1 + (6625.0L/3072.0L)*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) + (1.0L/96.0L)*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) + (75.0L/64.0L)*x_0*pow(x_1, 2) - 1375.0L/384.0L*x_0*x_1 + (12125.0L/3072.0L)*x_0 - 1.0L/4800.0L*pow(x_1, 6) - 1.0L/960.0L*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) - 25.0L/128.0L*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) - 8875.0L/3072.0L*x_1 + 4375.0L/1536.0L;
       return __pp_r308___result;
    }
    static inline O __pp_r309__(const I &x_0, const I &x_1) {
       O __pp_r309___result;
       __pp_r309___result = -1.0L/800.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 79.0L/3600.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (3.0L/80.0L)*pow(x_0, 4)*x_1 - 71.0L/480.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (11.0L/90.0L)*pow(x_0, 3)*x_1 - 947.0L/2160.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 101.0L/480.0L*pow(x_0, 2)*x_1 - 1187.0L/3840.0L*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) + (37.0L/1440.0L)*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) + (2803.0L/2880.0L)*x_0*pow(x_1, 2) - 21439.0L/11520.0L*x_0*x_1 + (128707.0L/115200.0L)*x_0 + (1.0L/4800.0L)*pow(x_1, 6) - 29.0L/4800.0L*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) - 289.0L/640.0L*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) - 209569.0L/76800.0L*x_1 + 574799.0L/307200.0L;
       return __pp_r309___result;
    }
    static inline O __pp_r310__(const I &x_0, const I &x_1) {
       O __pp_r310___result;
       __pp_r310___result = -1.0L/800.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 79.0L/3600.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (3.0L/80.0L)*pow(x_0, 4)*x_1 - 71.0L/480.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (11.0L/90.0L)*pow(x_0, 3)*x_1 - 947.0L/2160.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 101.0L/480.0L*pow(x_0, 2)*x_1 - 1187.0L/3840.0L*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) + (37.0L/1440.0L)*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) + (2803.0L/2880.0L)*x_0*pow(x_1, 2) - 21439.0L/11520.0L*x_0*x_1 + (128707.0L/115200.0L)*x_0 + (1.0L/4800.0L)*pow(x_1, 6) - 29.0L/4800.0L*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) - 289.0L/640.0L*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) - 209569.0L/76800.0L*x_1 + 574799.0L/307200.0L;
       return __pp_r310___result;
    }
    static inline O __pp_r311__(const I &x_0, const I &x_1) {
       O __pp_r311___result;
       __pp_r311___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (1.0L/144.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/144.0L*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) - 25.0L/72.0L*pow(x_0, 3)*x_1 + (125.0L/216.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 125.0L/72.0L*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) + (5.0L/144.0L)*x_0*pow(x_1, 4) - 25.0L/72.0L*x_0*pow(x_1, 3) + (125.0L/72.0L)*x_0*pow(x_1, 2) - 625.0L/144.0L*x_0*x_1 + (625.0L/144.0L)*x_0 + (1.0L/4320.0L)*pow(x_1, 6) - 1.0L/144.0L*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) - 125.0L/216.0L*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) - 625.0L/144.0L*x_1 + 3125.0L/864.0L;
       return __pp_r311___result;
    }
    static inline O __pp_r312__(const I &x_0, const I &x_1) {
       O __pp_r312___result;
       __pp_r312___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (1.0L/144.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/144.0L*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) - 25.0L/72.0L*pow(x_0, 3)*x_1 + (125.0L/216.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 125.0L/72.0L*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) + (5.0L/144.0L)*x_0*pow(x_1, 4) - 25.0L/72.0L*x_0*pow(x_1, 3) + (125.0L/72.0L)*x_0*pow(x_1, 2) - 625.0L/144.0L*x_0*x_1 + (625.0L/144.0L)*x_0 + (1.0L/4320.0L)*pow(x_1, 6) - 1.0L/144.0L*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) - 125.0L/216.0L*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) - 625.0L/144.0L*x_1 + 3125.0L/864.0L;
       return __pp_r312___result;
    }
    static inline O __pp_r313__(const I &x_0, const I &x_1) {
       O __pp_r313___result;
       __pp_r313___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 1.0L/144.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) - 25.0L/72.0L*pow(x_0, 3)*x_1 - 125.0L/216.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (125.0L/72.0L)*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) - 5.0L/144.0L*x_0*pow(x_1, 4) - 25.0L/72.0L*x_0*pow(x_1, 3) - 125.0L/72.0L*x_0*pow(x_1, 2) - 625.0L/144.0L*x_0*x_1 - 625.0L/144.0L*x_0 + (1.0L/4320.0L)*pow(x_1, 6) + (1.0L/144.0L)*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) + (125.0L/216.0L)*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) + (625.0L/144.0L)*x_1 + 3125.0L/864.0L;
       return __pp_r313___result;
    }
    static inline O __pp_r314__(const I &x_0, const I &x_1) {
       O __pp_r314___result;
       __pp_r314___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 1.0L/144.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) - 25.0L/72.0L*pow(x_0, 3)*x_1 - 125.0L/216.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (125.0L/72.0L)*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) - 5.0L/144.0L*x_0*pow(x_1, 4) - 25.0L/72.0L*x_0*pow(x_1, 3) - 125.0L/72.0L*x_0*pow(x_1, 2) - 625.0L/144.0L*x_0*x_1 - 625.0L/144.0L*x_0 + (1.0L/4320.0L)*pow(x_1, 6) + (1.0L/144.0L)*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) + (125.0L/216.0L)*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) + (625.0L/144.0L)*x_1 + 3125.0L/864.0L;
       return __pp_r314___result;
    }
    static inline O __pp_r315__(const I &x_0, const I &x_1) {
       O __pp_r315___result;
       __pp_r315___result = -1.0L/800.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (79.0L/3600.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 3.0L/80.0L*pow(x_0, 4)*x_1 - 71.0L/480.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (11.0L/90.0L)*pow(x_0, 3)*x_1 + (947.0L/2160.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (101.0L/480.0L)*pow(x_0, 2)*x_1 - 1187.0L/3840.0L*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) - 37.0L/1440.0L*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) - 2803.0L/2880.0L*x_0*pow(x_1, 2) - 21439.0L/11520.0L*x_0*x_1 - 128707.0L/115200.0L*x_0 + (1.0L/4800.0L)*pow(x_1, 6) + (29.0L/4800.0L)*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) + (289.0L/640.0L)*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) + (209569.0L/76800.0L)*x_1 + 574799.0L/307200.0L;
       return __pp_r315___result;
    }
    static inline O __pp_r316__(const I &x_0, const I &x_1) {
       O __pp_r316___result;
       __pp_r316___result = -1.0L/800.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (79.0L/3600.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 3.0L/80.0L*pow(x_0, 4)*x_1 - 71.0L/480.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (11.0L/90.0L)*pow(x_0, 3)*x_1 + (947.0L/2160.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (101.0L/480.0L)*pow(x_0, 2)*x_1 - 1187.0L/3840.0L*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) - 37.0L/1440.0L*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) - 2803.0L/2880.0L*x_0*pow(x_1, 2) - 21439.0L/11520.0L*x_0*x_1 - 128707.0L/115200.0L*x_0 + (1.0L/4800.0L)*pow(x_1, 6) + (29.0L/4800.0L)*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) + (289.0L/640.0L)*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) + (209569.0L/76800.0L)*x_1 + 574799.0L/307200.0L;
       return __pp_r316___result;
    }
    static inline O __pp_r317__(const I &x_0, const I &x_1) {
       O __pp_r317___result;
       __pp_r317___result = (1.0L/4800.0L)*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 7.0L/960.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/32.0L)*pow(x_0, 4)*x_1 + (65.0L/768.0L)*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/12.0L*pow(x_0, 3)*pow(x_1, 2) - 35.0L/96.0L*pow(x_0, 3)*x_1 - 75.0L/128.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (325.0L/192.0L)*pow(x_0, 2)*x_1 + (6625.0L/3072.0L)*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) - 1.0L/96.0L*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) - 75.0L/64.0L*x_0*pow(x_1, 2) - 1375.0L/384.0L*x_0*x_1 - 12125.0L/3072.0L*x_0 - 1.0L/4800.0L*pow(x_1, 6) + (1.0L/960.0L)*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) + (25.0L/128.0L)*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) + (8875.0L/3072.0L)*x_1 + 4375.0L/1536.0L;
       return __pp_r317___result;
    }
    static inline O __pp_r318__(const I &x_0, const I &x_1) {
       O __pp_r318___result;
       __pp_r318___result = -11.0L/8640.0L*pow(x_0, 6) + (1.0L/360.0L)*pow(x_0, 5)*x_1 + (311.0L/14400.0L)*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) - 59.0L/1440.0L*pow(x_0, 4)*x_1 - 1729.0L/11520.0L*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (151.0L/1440.0L)*pow(x_0, 3)*x_1 + (7451.0L/17280.0L)*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) + (481.0L/2880.0L)*pow(x_0, 2)*x_1 - 14869.0L/46080.0L*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) - 1.0L/720.0L*x_0*pow(x_1, 4) - 7.0L/180.0L*x_0*pow(x_1, 3) - 589.0L/1440.0L*x_0*pow(x_1, 2) - 12689.0L/11520.0L*x_0*x_1 - 166789.0L/230400.0L*x_0 - 1.0L/4320.0L*pow(x_1, 6) + (1.0L/7200.0L)*pow(x_1, 5) + (7.0L/1440.0L)*pow(x_1, 4) + (589.0L/8640.0L)*pow(x_1, 3) + (1379.0L/2880.0L)*pow(x_1, 2) + (73583.0L/57600.0L)*x_1 + 3048191.0L/2764800.0L;
       return __pp_r318___result;
    }
    static inline O __pp_r319__(const I &x_0, const I &x_1) {
       O __pp_r319___result;
       __pp_r319___result = -1.0L/576.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 + (157.0L/4800.0L)*pow(x_0, 5) - 1.0L/96.0L*pow(x_0, 4)*pow(x_1, 2) - 139.0L/1440.0L*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/10.0L)*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 + (5897.0L/5760.0L)*pow(x_0, 3) - 1.0L/96.0L*pow(x_0, 2)*pow(x_1, 4) - 19.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 4639.0L/2880.0L*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (13.0L/240.0L)*x_0*pow(x_1, 4) + (73.0L/180.0L)*x_0*pow(x_1, 3) + (219.0L/160.0L)*x_0*pow(x_1, 2) + (28271.0L/11520.0L)*x_0*x_1 + (162857.0L/76800.0L)*x_0 - 1.0L/1440.0L*pow(x_1, 6) - 79.0L/7200.0L*pow(x_1, 5) - 17.0L/160.0L*pow(x_1, 4) - 4531.0L/8640.0L*pow(x_1, 3) - 1247.0L/960.0L*pow(x_1, 2) - 90257.0L/57600.0L*x_1 - 731563.0L/921600.0L;
       return __pp_r319___result;
    }
    static inline O __pp_r320__(const I &x_0, const I &x_1) {
       O __pp_r320___result;
       __pp_r320___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 1.0L/144.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) - 25.0L/72.0L*pow(x_0, 3)*x_1 - 125.0L/216.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (125.0L/72.0L)*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) + (5.0L/288.0L)*x_0*pow(x_1, 4) - 25.0L/288.0L*x_0*pow(x_1, 3) - 625.0L/576.0L*x_0*pow(x_1, 2) - 8125.0L/2304.0L*x_0*x_1 - 18125.0L/4608.0L*x_0 + (11.0L/8640.0L)*pow(x_1, 6) + (7.0L/576.0L)*pow(x_1, 5) + (125.0L/2304.0L)*pow(x_1, 4) + (875.0L/3456.0L)*pow(x_1, 3) + (10625.0L/9216.0L)*pow(x_1, 2) + (26875.0L/9216.0L)*x_1 + 315625.0L/110592.0L;
       return __pp_r320___result;
    }
    static inline O __pp_r321__(const I &x_0, const I &x_1) {
       O __pp_r321___result;
       __pp_r321___result = -1.0L/800.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (79.0L/3600.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 3.0L/80.0L*pow(x_0, 4)*x_1 - 71.0L/480.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (11.0L/90.0L)*pow(x_0, 3)*x_1 + (947.0L/2160.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (101.0L/480.0L)*pow(x_0, 2)*x_1 - 1187.0L/3840.0L*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) + (19.0L/720.0L)*x_0*pow(x_1, 4) + (11.0L/360.0L)*x_0*pow(x_1, 3) - 29.0L/90.0L*x_0*pow(x_1, 2) - 377.0L/360.0L*x_0*x_1 - 10229.0L/14400.0L*x_0 + (1.0L/800.0L)*pow(x_1, 6) + (9.0L/800.0L)*pow(x_1, 5) + (19.0L/480.0L)*pow(x_1, 4) + (121.0L/960.0L)*pow(x_1, 3) + (2047.0L/3840.0L)*pow(x_1, 2) + (16699.0L/12800.0L)*x_1 + 42553.0L/38400.0L;
       return __pp_r321___result;
    }
    static inline O __pp_r322__(const I &x_0, const I &x_1) {
       O __pp_r322___result;
       __pp_r322___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 + (119.0L/3600.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 67.0L/720.0L*pow(x_0, 4)*x_1 - 373.0L/1440.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (41.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (17.0L/30.0L)*pow(x_0, 3)*x_1 + (2227.0L/2160.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 2257.0L/1440.0L*pow(x_0, 2)*x_1 - 24041.0L/11520.0L*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) + (59.0L/720.0L)*x_0*pow(x_1, 4) + (19.0L/40.0L)*x_0*pow(x_1, 3) + (131.0L/90.0L)*x_0*pow(x_1, 2) + (301.0L/120.0L)*x_0*x_1 + (30731.0L/14400.0L)*x_0 + (17.0L/21600.0L)*pow(x_1, 6) + (1.0L/7200.0L)*pow(x_1, 5) - 103.0L/1440.0L*pow(x_1, 4) - 4031.0L/8640.0L*pow(x_1, 3) - 14339.0L/11520.0L*pow(x_1, 2) - 177389.0L/115200.0L*x_1 - 272383.0L/345600.0L;
       return __pp_r322___result;
    }
    static inline O __pp_r323__(const I &x_0, const I &x_1) {
       O __pp_r323___result;
       __pp_r323___result = -1.0L/576.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 + (157.0L/4800.0L)*pow(x_0, 5) - 1.0L/96.0L*pow(x_0, 4)*pow(x_1, 2) - 139.0L/1440.0L*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/10.0L)*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 + (5897.0L/5760.0L)*pow(x_0, 3) - 1.0L/96.0L*pow(x_0, 2)*pow(x_1, 4) - 19.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 4639.0L/2880.0L*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (13.0L/240.0L)*x_0*pow(x_1, 4) + (73.0L/180.0L)*x_0*pow(x_1, 3) + (219.0L/160.0L)*x_0*pow(x_1, 2) + (28271.0L/11520.0L)*x_0*x_1 + (162857.0L/76800.0L)*x_0 - 1.0L/1440.0L*pow(x_1, 6) - 79.0L/7200.0L*pow(x_1, 5) - 17.0L/160.0L*pow(x_1, 4) - 4531.0L/8640.0L*pow(x_1, 3) - 1247.0L/960.0L*pow(x_1, 2) - 90257.0L/57600.0L*x_1 - 731563.0L/921600.0L;
       return __pp_r323___result;
    }
    static inline O __pp_r324__(const I &x_0, const I &x_1) {
       O __pp_r324___result;
       __pp_r324___result = -1.0L/1440.0L*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 + (7.0L/800.0L)*pow(x_0, 5) - 1.0L/96.0L*pow(x_0, 4)*pow(x_1, 2) - 47.0L/720.0L*pow(x_0, 4)*x_1 - 7.0L/60.0L*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/10.0L)*pow(x_0, 3)*pow(x_1, 2) + (41.0L/90.0L)*pow(x_0, 3)*x_1 + (1801.0L/2880.0L)*pow(x_0, 3) - 1.0L/96.0L*pow(x_0, 2)*pow(x_1, 4) - 19.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 2117.0L/1440.0L*pow(x_0, 2)*x_1 - 367.0L/240.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (13.0L/240.0L)*x_0*pow(x_1, 4) + (73.0L/180.0L)*x_0*pow(x_1, 3) + (219.0L/160.0L)*x_0*pow(x_1, 2) + (1691.0L/720.0L)*x_0*x_1 + (65431.0L/38400.0L)*x_0 - 1.0L/1440.0L*pow(x_1, 6) - 79.0L/7200.0L*pow(x_1, 5) - 17.0L/160.0L*pow(x_1, 4) - 4531.0L/8640.0L*pow(x_1, 3) - 1247.0L/960.0L*pow(x_1, 2) - 176869.0L/115200.0L*x_1 - 77321.0L/115200.0L;
       return __pp_r324___result;
    }
    static inline O __pp_r325__(const I &x_0, const I &x_1) {
       O __pp_r325___result;
       __pp_r325___result = (49.0L/21600.0L)*pow(x_0, 6) - 3.0L/400.0L*pow(x_0, 5)*x_1 - 1.0L/32.0L*pow(x_0, 5) + (1.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 4)*x_1 + (13.0L/120.0L)*pow(x_0, 4) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/180.0L)*pow(x_0, 3)*x_1 - 143.0L/2880.0L*pow(x_0, 3) - 11.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/18.0L*pow(x_0, 2)*pow(x_1, 3) - 17.0L/80.0L*pow(x_0, 2)*pow(x_1, 2) - 659.0L/1440.0L*pow(x_0, 2)*x_1 - 749.0L/1920.0L*pow(x_0, 2) + (1.0L/1200.0L)*x_0*pow(x_1, 5) + (1.0L/24.0L)*x_0*pow(x_1, 4) + (211.0L/720.0L)*x_0*pow(x_1, 3) + (69.0L/80.0L)*x_0*pow(x_1, 2) + (6967.0L/5760.0L)*x_0*x_1 + (5213.0L/7680.0L)*x_0 - 7.0L/10800.0L*pow(x_1, 6) - 7.0L/720.0L*pow(x_1, 5) - 59.0L/640.0L*pow(x_1, 4) - 1901.0L/4320.0L*pow(x_1, 3) - 7789.0L/7680.0L*pow(x_1, 2) - 5891.0L/5760.0L*x_1 - 132137.0L/460800.0L;
       return __pp_r325___result;
    }
    static inline O __pp_r326__(const I &x_0, const I &x_1) {
       O __pp_r326___result;
       __pp_r326___result = (2.0L/675.0L)*pow(x_0, 6) - 7.0L/600.0L*pow(x_0, 5)*x_1 - 7.0L/160.0L*pow(x_0, 5) + (1.0L/90.0L)*pow(x_0, 4)*pow(x_1, 2) + (7.0L/72.0L)*pow(x_0, 4)*x_1 + (97.0L/480.0L)*pow(x_0, 4) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/8.0L*pow(x_0, 3)*pow(x_1, 2) - 133.0L/360.0L*pow(x_0, 3)*x_1 - 1223.0L/2880.0L*pow(x_0, 3) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) + (961.0L/1440.0L)*pow(x_0, 2)*x_1 + (871.0L/1920.0L)*pow(x_0, 2) - 1.0L/300.0L*x_0*pow(x_1, 5) - 1.0L/48.0L*x_0*pow(x_1, 4) - 59.0L/720.0L*x_0*pow(x_1, 3) - 21.0L/80.0L*x_0*pow(x_1, 2) - 2753.0L/5760.0L*x_0*x_1 - 2563.0L/7680.0L*x_0 + (1.0L/21600.0L)*pow(x_1, 6) + (1.0L/360.0L)*pow(x_1, 5) + (1.0L/640.0L)*pow(x_1, 4) - 281.0L/4320.0L*pow(x_1, 3) - 1309.0L/7680.0L*pow(x_1, 2) - 59.0L/5760.0L*x_1 + 101143.0L/460800.0L;
       return __pp_r326___result;
    }
    static inline O __pp_r327__(const I &x_0, const I &x_1) {
       O __pp_r327___result;
       __pp_r327___result = -1.0L/360.0L*pow(x_0, 5)*x_1 - 91.0L/7200.0L*pow(x_0, 5) + (7.0L/360.0L)*pow(x_0, 4)*x_1 + (19.0L/288.0L)*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) - 17.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 7.0L/72.0L*pow(x_0, 3)*x_1 - 185.0L/1728.0L*pow(x_0, 3) + (11.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (55.0L/288.0L)*pow(x_0, 2)*x_1 + (53.0L/1440.0L)*pow(x_0, 2) - 1.0L/360.0L*x_0*pow(x_1, 5) - 1.0L/90.0L*x_0*pow(x_1, 4) - 1.0L/72.0L*x_0*pow(x_1, 3) - 7.0L/288.0L*x_0*pow(x_1, 2) - 11.0L/180.0L*x_0*x_1 - 4831.0L/115200.0L*x_0 + (13.0L/7200.0L)*pow(x_1, 5) - 1.0L/144.0L*pow(x_1, 4) - 181.0L/1728.0L*pow(x_1, 3) - 791.0L/2880.0L*pow(x_1, 2) - 17987.0L/115200.0L*x_1 + 9289.0L/69120.0L;
       return __pp_r327___result;
    }
    static inline O __pp_r328__(const I &x_0, const I &x_1) {
       O __pp_r328___result;
       __pp_r328___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 + (119.0L/3600.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 67.0L/720.0L*pow(x_0, 4)*x_1 - 373.0L/1440.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (41.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (17.0L/30.0L)*pow(x_0, 3)*x_1 + (2227.0L/2160.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 2257.0L/1440.0L*pow(x_0, 2)*x_1 - 24041.0L/11520.0L*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) + (59.0L/720.0L)*x_0*pow(x_1, 4) + (19.0L/40.0L)*x_0*pow(x_1, 3) + (131.0L/90.0L)*x_0*pow(x_1, 2) + (301.0L/120.0L)*x_0*x_1 + (30731.0L/14400.0L)*x_0 + (17.0L/21600.0L)*pow(x_1, 6) + (1.0L/7200.0L)*pow(x_1, 5) - 103.0L/1440.0L*pow(x_1, 4) - 4031.0L/8640.0L*pow(x_1, 3) - 14339.0L/11520.0L*pow(x_1, 2) - 177389.0L/115200.0L*x_1 - 272383.0L/345600.0L;
       return __pp_r328___result;
    }
    static inline O __pp_r329__(const I &x_0, const I &x_1) {
       O __pp_r329___result;
       __pp_r329___result = -29.0L/43200.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 + (131.0L/14400.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 89.0L/1440.0L*pow(x_0, 4)*x_1 - 1319.0L/11520.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (41.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (227.0L/480.0L)*pow(x_0, 3)*x_1 + (10931.0L/17280.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 4109.0L/2880.0L*pow(x_0, 2)*x_1 - 69839.0L/46080.0L*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) + (59.0L/720.0L)*x_0*pow(x_1, 4) + (19.0L/40.0L)*x_0*pow(x_1, 3) + (131.0L/90.0L)*x_0*pow(x_1, 2) + (9227.0L/3840.0L)*x_0*x_1 + (395711.0L/230400.0L)*x_0 + (17.0L/21600.0L)*pow(x_1, 6) + (1.0L/7200.0L)*pow(x_1, 5) - 103.0L/1440.0L*pow(x_1, 4) - 4031.0L/8640.0L*pow(x_1, 3) - 14339.0L/11520.0L*pow(x_1, 2) - 10859.0L/7200.0L*x_1 - 1840079.0L/2764800.0L;
       return __pp_r329___result;
    }
    static inline O __pp_r330__(const I &x_0, const I &x_1) {
       O __pp_r330___result;
       __pp_r330___result = (11.0L/4800.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 - 89.0L/2880.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/288.0L)*pow(x_0, 4)*x_1 + (1273.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 3)*x_1 - 733.0L/17280.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/36.0L*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 1193.0L/2880.0L*pow(x_0, 2)*x_1 - 17351.0L/46080.0L*pow(x_0, 2) + (19.0L/3600.0L)*x_0*pow(x_1, 5) + (5.0L/72.0L)*x_0*pow(x_1, 4) + (29.0L/80.0L)*x_0*pow(x_1, 3) + (1367.0L/1440.0L)*x_0*pow(x_1, 2) + (4853.0L/3840.0L)*x_0*x_1 + (31903.0L/46080.0L)*x_0 + (1.0L/1200.0L)*pow(x_1, 6) + (1.0L/720.0L)*pow(x_1, 5) - 331.0L/5760.0L*pow(x_1, 4) - 1651.0L/4320.0L*pow(x_1, 3) - 22117.0L/23040.0L*pow(x_1, 2) - 22939.0L/23040.0L*x_1 - 777197.0L/2764800.0L;
       return __pp_r330___result;
    }
    static inline O __pp_r331__(const I &x_0, const I &x_1) {
       O __pp_r331___result;
       __pp_r331___result = (43.0L/14400.0L)*pow(x_0, 6) - 41.0L/3600.0L*pow(x_0, 5)*x_1 - 25.0L/576.0L*pow(x_0, 5) + (1.0L/80.0L)*pow(x_0, 4)*pow(x_1, 2) + (29.0L/288.0L)*pow(x_0, 4)*x_1 + (2353.0L/11520.0L)*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/9.0L*pow(x_0, 3)*pow(x_1, 2) - 169.0L/480.0L*pow(x_0, 3)*x_1 - 7213.0L/17280.0L*pow(x_0, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 4) + (7.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (193.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (2047.0L/2880.0L)*pow(x_0, 2)*x_1 + (21529.0L/46080.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (1.0L/144.0L)*x_0*pow(x_1, 4) - 1.0L/80.0L*x_0*pow(x_1, 3) - 253.0L/1440.0L*x_0*pow(x_1, 2) - 1627.0L/3840.0L*x_0*x_1 - 14753.0L/46080.0L*x_0 + (11.0L/7200.0L)*pow(x_1, 6) + (1.0L/72.0L)*pow(x_1, 5) + (209.0L/5760.0L)*pow(x_1, 4) - 31.0L/4320.0L*pow(x_1, 3) - 2677.0L/23040.0L*pow(x_1, 2) + (389.0L/23040.0L)*x_1 + 622483.0L/2764800.0L;
       return __pp_r331___result;
    }
    static inline O __pp_r332__(const I &x_0, const I &x_1) {
       O __pp_r332___result;
       __pp_r332___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 - 59.0L/4800.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 4)*x_1 + (157.0L/2304.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) - 23.0L/288.0L*pow(x_0, 3)*x_1 - 115.0L/1152.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (7.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (15.0L/64.0L)*pow(x_0, 2)*x_1 + (2321.0L/46080.0L)*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (1.0L/60.0L)*x_0*pow(x_1, 4) + (1.0L/18.0L)*x_0*pow(x_1, 3) + (1.0L/16.0L)*x_0*pow(x_1, 2) - 79.0L/11520.0L*x_0*x_1 - 2179.0L/76800.0L*x_0 + (1.0L/675.0L)*pow(x_1, 6) + (31.0L/2400.0L)*pow(x_1, 5) + (1.0L/36.0L)*pow(x_1, 4) - 3.0L/64.0L*pow(x_1, 3) - 2539.0L/11520.0L*pow(x_1, 2) - 2477.0L/19200.0L*x_1 + 77437.0L/552960.0L;
       return __pp_r332___result;
    }
    static inline O __pp_r333__(const I &x_0, const I &x_1) {
       O __pp_r333___result;
       __pp_r333___result = (2.0L/675.0L)*pow(x_0, 6) - 7.0L/600.0L*pow(x_0, 5)*x_1 - 7.0L/160.0L*pow(x_0, 5) + (1.0L/90.0L)*pow(x_0, 4)*pow(x_1, 2) + (7.0L/72.0L)*pow(x_0, 4)*x_1 + (97.0L/480.0L)*pow(x_0, 4) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/8.0L*pow(x_0, 3)*pow(x_1, 2) - 133.0L/360.0L*pow(x_0, 3)*x_1 - 1223.0L/2880.0L*pow(x_0, 3) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) + (961.0L/1440.0L)*pow(x_0, 2)*x_1 + (871.0L/1920.0L)*pow(x_0, 2) - 3.0L/400.0L*x_0*pow(x_1, 5) - 5.0L/96.0L*x_0*pow(x_1, 4) - 253.0L/1440.0L*x_0*pow(x_1, 3) - 129.0L/320.0L*x_0*pow(x_1, 2) - 6721.0L/11520.0L*x_0*x_1 - 1403.0L/3840.0L*x_0 - 43.0L/43200.0L*pow(x_1, 6) + (11.0L/2880.0L)*pow(x_1, 5) + (57.0L/1280.0L)*pow(x_1, 4) + (1711.0L/17280.0L)*pow(x_1, 3) + (1567.0L/15360.0L)*pow(x_1, 2) + (9491.0L/46080.0L)*x_1 + 264251.0L/921600.0L;
       return __pp_r333___result;
    }
    static inline O __pp_r334__(const I &x_0, const I &x_1) {
       O __pp_r334___result;
       __pp_r334___result = -1.0L/360.0L*pow(x_0, 5)*x_1 - 91.0L/7200.0L*pow(x_0, 5) + (7.0L/360.0L)*pow(x_0, 4)*x_1 + (19.0L/288.0L)*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) - 17.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 7.0L/72.0L*pow(x_0, 3)*x_1 - 185.0L/1728.0L*pow(x_0, 3) + (11.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (55.0L/288.0L)*pow(x_0, 2)*x_1 + (53.0L/1440.0L)*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) - 61.0L/1440.0L*x_0*pow(x_1, 4) - 31.0L/288.0L*x_0*pow(x_1, 3) - 95.0L/576.0L*x_0*pow(x_1, 2) - 1919.0L/11520.0L*x_0*x_1 - 2119.0L/28800.0L*x_0 - 1.0L/960.0L*pow(x_1, 6) + (41.0L/14400.0L)*pow(x_1, 5) + (83.0L/2304.0L)*pow(x_1, 4) + (205.0L/3456.0L)*pow(x_1, 3) - 101.0L/46080.0L*pow(x_1, 2) + (13841.0L/230400.0L)*x_1 + 111491.0L/552960.0L;
       return __pp_r334___result;
    }
    static inline O __pp_r335__(const I &x_0, const I &x_1) {
       O __pp_r335___result;
       __pp_r335___result = -1.0L/4320.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 71.0L/7200.0L*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/180.0L)*pow(x_0, 4)*x_1 + (5.0L/96.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 7.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 1.0L/24.0L*pow(x_0, 3)*x_1 - 121.0L/1728.0L*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (23.0L/288.0L)*pow(x_0, 2)*x_1 - 3.0L/160.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) - 41.0L/1440.0L*x_0*pow(x_1, 4) - 5.0L/96.0L*x_0*pow(x_1, 3) - 31.0L/576.0L*x_0*pow(x_1, 2) - 71.0L/1280.0L*x_0*x_1 - 839.0L/28800.0L*x_0 - 11.0L/8640.0L*pow(x_1, 6) + (1.0L/14400.0L)*pow(x_1, 5) + (17.0L/768.0L)*pow(x_1, 4) + (77.0L/3456.0L)*pow(x_1, 3) - 887.0L/15360.0L*pow(x_1, 2) + (3601.0L/230400.0L)*x_1 + 34433.0L/184320.0L;
       return __pp_r335___result;
    }
    static inline O __pp_r336__(const I &x_0, const I &x_1) {
       O __pp_r336___result;
       __pp_r336___result = -37.0L/21600.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/800.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/45.0L*pow(x_0, 4)*x_1 + (5.0L/288.0L)*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/36.0L)*pow(x_0, 3)*x_1 - 7.0L/576.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/144.0L*pow(x_0, 2)*x_1 - 841.0L/11520.0L*pow(x_0, 2) - 19.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/40.0L*x_0*pow(x_1, 4) - 5.0L/144.0L*x_0*pow(x_1, 3) - 1.0L/96.0L*x_0*pow(x_1, 2) - 7.0L/5760.0L*x_0*x_1 - 77.0L/38400.0L*x_0 - 7.0L/5400.0L*pow(x_1, 6) - 1.0L/3600.0L*pow(x_1, 5) + (23.0L/1152.0L)*pow(x_1, 4) + (13.0L/864.0L)*pow(x_1, 3) - 1643.0L/23040.0L*pow(x_1, 2) + (119.0L/57600.0L)*x_1 + 50087.0L/276480.0L;
       return __pp_r336___result;
    }
    static inline O __pp_r337__(const I &x_0, const I &x_1) {
       O __pp_r337___result;
       __pp_r337___result = -37.0L/21600.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (23.0L/2400.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/45.0L*pow(x_0, 4)*x_1 - 1.0L/288.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/36.0L)*pow(x_0, 3)*x_1 + (5.0L/576.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/144.0L*pow(x_0, 2)*x_1 - 961.0L/11520.0L*pow(x_0, 2) - 19.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/40.0L*x_0*pow(x_1, 4) - 5.0L/144.0L*x_0*pow(x_1, 3) - 1.0L/96.0L*x_0*pow(x_1, 2) - 7.0L/5760.0L*x_0*x_1 + (23.0L/38400.0L)*x_0 - 7.0L/5400.0L*pow(x_1, 6) - 1.0L/3600.0L*pow(x_1, 5) + (23.0L/1152.0L)*pow(x_1, 4) + (13.0L/864.0L)*pow(x_1, 3) - 1643.0L/23040.0L*pow(x_1, 2) + (119.0L/57600.0L)*x_1 + 10003.0L/55296.0L;
       return __pp_r337___result;
    }
    static inline O __pp_r338__(const I &x_0, const I &x_1) {
       O __pp_r338___result;
       __pp_r338___result = (2.0L/675.0L)*pow(x_0, 6) - 7.0L/600.0L*pow(x_0, 5)*x_1 - 7.0L/160.0L*pow(x_0, 5) + (1.0L/90.0L)*pow(x_0, 4)*pow(x_1, 2) + (7.0L/72.0L)*pow(x_0, 4)*x_1 + (97.0L/480.0L)*pow(x_0, 4) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/8.0L*pow(x_0, 3)*pow(x_1, 2) - 133.0L/360.0L*pow(x_0, 3)*x_1 - 1223.0L/2880.0L*pow(x_0, 3) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) + (961.0L/1440.0L)*pow(x_0, 2)*x_1 + (871.0L/1920.0L)*pow(x_0, 2) - 3.0L/400.0L*x_0*pow(x_1, 5) - 5.0L/96.0L*x_0*pow(x_1, 4) - 253.0L/1440.0L*x_0*pow(x_1, 3) - 129.0L/320.0L*x_0*pow(x_1, 2) - 6721.0L/11520.0L*x_0*x_1 - 1403.0L/3840.0L*x_0 - 43.0L/43200.0L*pow(x_1, 6) + (11.0L/2880.0L)*pow(x_1, 5) + (57.0L/1280.0L)*pow(x_1, 4) + (1711.0L/17280.0L)*pow(x_1, 3) + (1567.0L/15360.0L)*pow(x_1, 2) + (9491.0L/46080.0L)*x_1 + 264251.0L/921600.0L;
       return __pp_r338___result;
    }
    static inline O __pp_r339__(const I &x_0, const I &x_1) {
       O __pp_r339___result;
       __pp_r339___result = -1.0L/360.0L*pow(x_0, 5)*x_1 - 91.0L/7200.0L*pow(x_0, 5) + (7.0L/360.0L)*pow(x_0, 4)*x_1 + (19.0L/288.0L)*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) - 17.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 7.0L/72.0L*pow(x_0, 3)*x_1 - 185.0L/1728.0L*pow(x_0, 3) + (11.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (55.0L/288.0L)*pow(x_0, 2)*x_1 + (53.0L/1440.0L)*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) - 61.0L/1440.0L*x_0*pow(x_1, 4) - 31.0L/288.0L*x_0*pow(x_1, 3) - 95.0L/576.0L*x_0*pow(x_1, 2) - 1919.0L/11520.0L*x_0*x_1 - 2119.0L/28800.0L*x_0 - 1.0L/960.0L*pow(x_1, 6) + (41.0L/14400.0L)*pow(x_1, 5) + (83.0L/2304.0L)*pow(x_1, 4) + (205.0L/3456.0L)*pow(x_1, 3) - 101.0L/46080.0L*pow(x_1, 2) + (13841.0L/230400.0L)*x_1 + 111491.0L/552960.0L;
       return __pp_r339___result;
    }
    static inline O __pp_r340__(const I &x_0, const I &x_1) {
       O __pp_r340___result;
       __pp_r340___result = -1.0L/4320.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 71.0L/7200.0L*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/180.0L)*pow(x_0, 4)*x_1 + (5.0L/96.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 7.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 1.0L/24.0L*pow(x_0, 3)*x_1 - 121.0L/1728.0L*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (23.0L/288.0L)*pow(x_0, 2)*x_1 - 3.0L/160.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) - 41.0L/1440.0L*x_0*pow(x_1, 4) - 5.0L/96.0L*x_0*pow(x_1, 3) - 31.0L/576.0L*x_0*pow(x_1, 2) - 71.0L/1280.0L*x_0*x_1 - 839.0L/28800.0L*x_0 - 11.0L/8640.0L*pow(x_1, 6) + (1.0L/14400.0L)*pow(x_1, 5) + (17.0L/768.0L)*pow(x_1, 4) + (77.0L/3456.0L)*pow(x_1, 3) - 887.0L/15360.0L*pow(x_1, 2) + (3601.0L/230400.0L)*x_1 + 34433.0L/184320.0L;
       return __pp_r340___result;
    }
    static inline O __pp_r341__(const I &x_0, const I &x_1) {
       O __pp_r341___result;
       __pp_r341___result = -37.0L/21600.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/800.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/45.0L*pow(x_0, 4)*x_1 + (5.0L/288.0L)*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/36.0L)*pow(x_0, 3)*x_1 - 7.0L/576.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/144.0L*pow(x_0, 2)*x_1 - 841.0L/11520.0L*pow(x_0, 2) - 19.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/40.0L*x_0*pow(x_1, 4) - 5.0L/144.0L*x_0*pow(x_1, 3) - 1.0L/96.0L*x_0*pow(x_1, 2) - 7.0L/5760.0L*x_0*x_1 - 77.0L/38400.0L*x_0 - 7.0L/5400.0L*pow(x_1, 6) - 1.0L/3600.0L*pow(x_1, 5) + (23.0L/1152.0L)*pow(x_1, 4) + (13.0L/864.0L)*pow(x_1, 3) - 1643.0L/23040.0L*pow(x_1, 2) + (119.0L/57600.0L)*x_1 + 50087.0L/276480.0L;
       return __pp_r341___result;
    }
    static inline O __pp_r342__(const I &x_0, const I &x_1) {
       O __pp_r342___result;
       __pp_r342___result = -37.0L/21600.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (23.0L/2400.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/45.0L*pow(x_0, 4)*x_1 - 1.0L/288.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/36.0L)*pow(x_0, 3)*x_1 + (5.0L/576.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/144.0L*pow(x_0, 2)*x_1 - 961.0L/11520.0L*pow(x_0, 2) - 19.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/40.0L*x_0*pow(x_1, 4) - 5.0L/144.0L*x_0*pow(x_1, 3) - 1.0L/96.0L*x_0*pow(x_1, 2) - 7.0L/5760.0L*x_0*x_1 + (23.0L/38400.0L)*x_0 - 7.0L/5400.0L*pow(x_1, 6) - 1.0L/3600.0L*pow(x_1, 5) + (23.0L/1152.0L)*pow(x_1, 4) + (13.0L/864.0L)*pow(x_1, 3) - 1643.0L/23040.0L*pow(x_1, 2) + (119.0L/57600.0L)*x_1 + 10003.0L/55296.0L;
       return __pp_r342___result;
    }
    static inline O __pp_r343__(const I &x_0, const I &x_1) {
       O __pp_r343___result;
       __pp_r343___result = -1.0L/4320.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (7.0L/2400.0L)*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/180.0L*pow(x_0, 4)*x_1 + (13.0L/1440.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/360.0L)*pow(x_0, 3)*x_1 - 11.0L/2880.0L*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/240.0L)*pow(x_0, 2)*pow(x_1, 2) + (17.0L/1440.0L)*pow(x_0, 2)*x_1 - 11.0L/144.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) - 13.0L/480.0L*x_0*pow(x_1, 4) - 59.0L/1440.0L*x_0*pow(x_1, 3) - 19.0L/960.0L*x_0*pow(x_1, 2) - 19.0L/2304.0L*x_0*x_1 - 29.0L/19200.0L*x_0 - 11.0L/8640.0L*pow(x_1, 6) - 1.0L/14400.0L*pow(x_1, 5) + (239.0L/11520.0L)*pow(x_1, 4) + (287.0L/17280.0L)*pow(x_1, 3) - 641.0L/9216.0L*pow(x_1, 2) + (719.0L/230400.0L)*x_1 + 500879.0L/2764800.0L;
       return __pp_r343___result;
    }
    static inline O __pp_r344__(const I &x_0, const I &x_1) {
       O __pp_r344___result;
       __pp_r344___result = -1.0L/2160.0L*pow(x_0, 6) + (31.0L/7200.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/80.0L*pow(x_0, 4)*x_1 + (1.0L/180.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 3)*x_1 + (7.0L/8640.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/480.0L*pow(x_0, 2)*x_1 - 23.0L/288.0L*pow(x_0, 2) - 1.0L/240.0L*x_0*pow(x_1, 5) - 29.0L/1440.0L*x_0*pow(x_1, 4) - 13.0L/480.0L*x_0*pow(x_1, 3) - 17.0L/2880.0L*x_0*pow(x_1, 2) - 1.0L/768.0L*x_0*x_1 - 7.0L/57600.0L*x_0 - 13.0L/8640.0L*pow(x_1, 6) - 7.0L/4800.0L*pow(x_1, 5) + (199.0L/11520.0L)*pow(x_1, 4) + (23.0L/1920.0L)*pow(x_1, 3) - 673.0L/9216.0L*pow(x_1, 2) + (133.0L/76800.0L)*x_1 + 500239.0L/2764800.0L;
       return __pp_r344___result;
    }
    static inline O __pp_r345__(const I &x_0, const I &x_1) {
       O __pp_r345___result;
       __pp_r345___result = (1.0L/400.0L)*pow(x_0, 6) - 2.0L/225.0L*pow(x_0, 5)*x_1 - 1.0L/7200.0L*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/720.0L*pow(x_0, 4)*x_1 + (1.0L/120.0L)*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/90.0L)*pow(x_0, 3)*x_1 - 1.0L/8640.0L*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/1440.0L*pow(x_0, 2)*x_1 - 51.0L/640.0L*pow(x_0, 2) - 17.0L/3600.0L*x_0*pow(x_1, 5) - 31.0L/1440.0L*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) - 19.0L/2880.0L*x_0*pow(x_1, 2) - 17.0L/11520.0L*x_0*x_1 - 1.0L/7200.0L*x_0 - 7.0L/4800.0L*pow(x_1, 6) - 19.0L/14400.0L*pow(x_1, 5) + (67.0L/3840.0L)*pow(x_1, 4) + (209.0L/17280.0L)*pow(x_1, 3) - 1121.0L/15360.0L*pow(x_1, 2) + (401.0L/230400.0L)*x_1 + 166747.0L/921600.0L;
       return __pp_r345___result;
    }
    static inline O __pp_r346__(const I &x_0, const I &x_1) {
       O __pp_r346___result;
       __pp_r346___result = -1.0L/4320.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (7.0L/2400.0L)*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/180.0L*pow(x_0, 4)*x_1 + (13.0L/1440.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/360.0L)*pow(x_0, 3)*x_1 - 11.0L/2880.0L*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/240.0L)*pow(x_0, 2)*pow(x_1, 2) + (17.0L/1440.0L)*pow(x_0, 2)*x_1 - 11.0L/144.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) - 13.0L/480.0L*x_0*pow(x_1, 4) - 59.0L/1440.0L*x_0*pow(x_1, 3) - 19.0L/960.0L*x_0*pow(x_1, 2) - 19.0L/2304.0L*x_0*x_1 - 29.0L/19200.0L*x_0 - 11.0L/8640.0L*pow(x_1, 6) - 1.0L/14400.0L*pow(x_1, 5) + (239.0L/11520.0L)*pow(x_1, 4) + (287.0L/17280.0L)*pow(x_1, 3) - 641.0L/9216.0L*pow(x_1, 2) + (719.0L/230400.0L)*x_1 + 500879.0L/2764800.0L;
       return __pp_r346___result;
    }
    static inline O __pp_r347__(const I &x_0, const I &x_1) {
       O __pp_r347___result;
       __pp_r347___result = -1.0L/2160.0L*pow(x_0, 6) + (31.0L/7200.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/80.0L*pow(x_0, 4)*x_1 + (1.0L/180.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 3)*x_1 + (7.0L/8640.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/480.0L*pow(x_0, 2)*x_1 - 23.0L/288.0L*pow(x_0, 2) - 1.0L/240.0L*x_0*pow(x_1, 5) - 29.0L/1440.0L*x_0*pow(x_1, 4) - 13.0L/480.0L*x_0*pow(x_1, 3) - 17.0L/2880.0L*x_0*pow(x_1, 2) - 1.0L/768.0L*x_0*x_1 - 7.0L/57600.0L*x_0 - 13.0L/8640.0L*pow(x_1, 6) - 7.0L/4800.0L*pow(x_1, 5) + (199.0L/11520.0L)*pow(x_1, 4) + (23.0L/1920.0L)*pow(x_1, 3) - 673.0L/9216.0L*pow(x_1, 2) + (133.0L/76800.0L)*x_1 + 500239.0L/2764800.0L;
       return __pp_r347___result;
    }
    static inline O __pp_r348__(const I &x_0, const I &x_1) {
       O __pp_r348___result;
       __pp_r348___result = (1.0L/400.0L)*pow(x_0, 6) - 2.0L/225.0L*pow(x_0, 5)*x_1 - 1.0L/7200.0L*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/720.0L*pow(x_0, 4)*x_1 + (1.0L/120.0L)*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/90.0L)*pow(x_0, 3)*x_1 - 1.0L/8640.0L*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/1440.0L*pow(x_0, 2)*x_1 - 51.0L/640.0L*pow(x_0, 2) - 17.0L/3600.0L*x_0*pow(x_1, 5) - 31.0L/1440.0L*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) - 19.0L/2880.0L*x_0*pow(x_1, 2) - 17.0L/11520.0L*x_0*x_1 - 1.0L/7200.0L*x_0 - 7.0L/4800.0L*pow(x_1, 6) - 19.0L/14400.0L*pow(x_1, 5) + (67.0L/3840.0L)*pow(x_1, 4) + (209.0L/17280.0L)*pow(x_1, 3) - 1121.0L/15360.0L*pow(x_1, 2) + (401.0L/230400.0L)*x_1 + 166747.0L/921600.0L;
       return __pp_r348___result;
    }
    static inline O __pp_r349__(const I &x_0, const I &x_1) {
       O __pp_r349___result;
       __pp_r349___result = -1.0L/2400.0L*pow(x_0, 6) + (1.0L/1800.0L)*pow(x_0, 5)*x_1 + (1.0L/225.0L)*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/90.0L*pow(x_0, 4)*x_1 + (11.0L/1920.0L)*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (13.0L/720.0L)*pow(x_0, 3)*x_1 + (1.0L/1080.0L)*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/720.0L*pow(x_0, 2)*x_1 - 613.0L/7680.0L*pow(x_0, 2) + (2.0L/225.0L)*x_0*pow(x_1, 5) + (1.0L/720.0L)*x_0*pow(x_1, 4) - 1.0L/90.0L*x_0*pow(x_1, 3) + (1.0L/1440.0L)*x_0*pow(x_1, 2) + (1.0L/5760.0L)*x_0*x_1 + (1.0L/57600.0L)*x_0 + (1.0L/400.0L)*pow(x_1, 6) - 1.0L/7200.0L*pow(x_1, 5) + (1.0L/120.0L)*pow(x_1, 4) - 1.0L/8640.0L*pow(x_1, 3) - 51.0L/640.0L*pow(x_1, 2) - 1.0L/115200.0L*x_1 + 83291.0L/460800.0L;
       return __pp_r349___result;
    }
    static inline O __pp_r350__(const I &x_0, const I &x_1) {
       O __pp_r350___result;
       __pp_r350___result = -1.0L/2160.0L*pow(x_0, 6) + (31.0L/7200.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/80.0L*pow(x_0, 4)*x_1 + (1.0L/180.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 3)*x_1 + (7.0L/8640.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/480.0L*pow(x_0, 2)*x_1 - 23.0L/288.0L*pow(x_0, 2) - 7.0L/720.0L*x_0*pow(x_1, 4) - 1.0L/60.0L*x_0*pow(x_1, 3) - 1.0L/1440.0L*x_0*pow(x_1, 2) + (1.0L/115200.0L)*x_0 - 1.0L/2160.0L*pow(x_1, 6) - 11.0L/2400.0L*pow(x_1, 5) + (1.0L/180.0L)*pow(x_1, 4) - 1.0L/960.0L*pow(x_1, 3) - 23.0L/288.0L*pow(x_1, 2) - 1.0L/38400.0L*x_1 + 15617.0L/86400.0L;
       return __pp_r350___result;
    }
    static inline O __pp_r351__(const I &x_0, const I &x_1) {
       O __pp_r351___result;
       __pp_r351___result = (11.0L/4320.0L)*pow(x_0, 6) - 1.0L/120.0L*pow(x_0, 5)*x_1 + (1.0L/144.0L)*pow(x_0, 4)*pow(x_1, 2) + (49.0L/5760.0L)*pow(x_0, 4) + (1.0L/80.0L)*pow(x_0, 3)*x_1 + (1.0L/144.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) - 367.0L/4608.0L*pow(x_0, 2) + (1.0L/120.0L)*x_0*pow(x_1, 5) - 1.0L/80.0L*x_0*pow(x_1, 3) + (11.0L/4320.0L)*pow(x_1, 6) + (49.0L/5760.0L)*pow(x_1, 4) - 367.0L/4608.0L*pow(x_1, 2) + 124937.0L/691200.0L;
       return __pp_r351___result;
    }
    static inline O __pp_r352__(const I &x_0, const I &x_1) {
       O __pp_r352___result;
       __pp_r352___result = (1.0L/400.0L)*pow(x_0, 6) - 2.0L/225.0L*pow(x_0, 5)*x_1 - 1.0L/7200.0L*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/720.0L*pow(x_0, 4)*x_1 + (1.0L/120.0L)*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/90.0L)*pow(x_0, 3)*x_1 - 1.0L/8640.0L*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/1440.0L*pow(x_0, 2)*x_1 - 51.0L/640.0L*pow(x_0, 2) - 1.0L/1800.0L*x_0*pow(x_1, 5) - 1.0L/90.0L*x_0*pow(x_1, 4) - 13.0L/720.0L*x_0*pow(x_1, 3) - 1.0L/720.0L*x_0*pow(x_1, 2) - 1.0L/5760.0L*x_0*x_1 - 1.0L/115200.0L*x_0 - 1.0L/2400.0L*pow(x_1, 6) - 1.0L/225.0L*pow(x_1, 5) + (11.0L/1920.0L)*pow(x_1, 4) - 1.0L/1080.0L*pow(x_1, 3) - 613.0L/7680.0L*pow(x_1, 2) - 1.0L/57600.0L*x_1 + 83291.0L/460800.0L;
       return __pp_r352___result;
    }
    static inline O __pp_r353__(const I &x_0, const I &x_1) {
       O __pp_r353___result;
       __pp_r353___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 1.0L/144.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) - 25.0L/72.0L*pow(x_0, 3)*x_1 - 125.0L/216.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (125.0L/72.0L)*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) - 5.0L/144.0L*x_0*pow(x_1, 4) - 25.0L/72.0L*x_0*pow(x_1, 3) - 125.0L/72.0L*x_0*pow(x_1, 2) - 625.0L/144.0L*x_0*x_1 - 625.0L/144.0L*x_0 + (1.0L/4320.0L)*pow(x_1, 6) + (1.0L/144.0L)*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) + (125.0L/216.0L)*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) + (625.0L/144.0L)*x_1 + 3125.0L/864.0L;
       return __pp_r353___result;
    }
    static inline O __pp_r354__(const I &x_0, const I &x_1) {
       O __pp_r354___result;
       __pp_r354___result = -1.0L/800.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (79.0L/3600.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 3.0L/80.0L*pow(x_0, 4)*x_1 - 71.0L/480.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (11.0L/90.0L)*pow(x_0, 3)*x_1 + (947.0L/2160.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (101.0L/480.0L)*pow(x_0, 2)*x_1 - 1187.0L/3840.0L*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) - 37.0L/1440.0L*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) - 2803.0L/2880.0L*x_0*pow(x_1, 2) - 21439.0L/11520.0L*x_0*x_1 - 128707.0L/115200.0L*x_0 + (1.0L/4800.0L)*pow(x_1, 6) + (29.0L/4800.0L)*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) + (289.0L/640.0L)*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) + (209569.0L/76800.0L)*x_1 + 574799.0L/307200.0L;
       return __pp_r354___result;
    }
    static inline O __pp_r355__(const I &x_0, const I &x_1) {
       O __pp_r355___result;
       __pp_r355___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 29.0L/14400.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/160.0L*pow(x_0, 4)*x_1 - 13.0L/3840.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 + (691.0L/17280.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (337.0L/960.0L)*pow(x_0, 2)*x_1 + (4027.0L/15360.0L)*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) - 37.0L/1440.0L*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) - 2803.0L/2880.0L*x_0*pow(x_1, 2) - 11327.0L/5760.0L*x_0*x_1 - 353399.0L/230400.0L*x_0 + (1.0L/4800.0L)*pow(x_1, 6) + (29.0L/4800.0L)*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) + (289.0L/640.0L)*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) + (211999.0L/76800.0L)*x_1 + 38279.0L/19200.0L;
       return __pp_r355___result;
    }
    static inline O __pp_r356__(const I &x_0, const I &x_1) {
       O __pp_r356___result;
       __pp_r356___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 1.0L/144.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 5.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) - 25.0L/72.0L*pow(x_0, 3)*x_1 - 125.0L/216.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (125.0L/72.0L)*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) - 5.0L/144.0L*x_0*pow(x_1, 4) - 25.0L/72.0L*x_0*pow(x_1, 3) - 125.0L/72.0L*x_0*pow(x_1, 2) - 625.0L/144.0L*x_0*x_1 - 625.0L/144.0L*x_0 + (1.0L/4320.0L)*pow(x_1, 6) + (1.0L/144.0L)*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) + (125.0L/216.0L)*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) + (625.0L/144.0L)*x_1 + 3125.0L/864.0L;
       return __pp_r356___result;
    }
    static inline O __pp_r357__(const I &x_0, const I &x_1) {
       O __pp_r357___result;
       __pp_r357___result = -1.0L/800.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (79.0L/3600.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 3.0L/80.0L*pow(x_0, 4)*x_1 - 71.0L/480.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (11.0L/90.0L)*pow(x_0, 3)*x_1 + (947.0L/2160.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (101.0L/480.0L)*pow(x_0, 2)*x_1 - 1187.0L/3840.0L*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) - 37.0L/1440.0L*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) - 2803.0L/2880.0L*x_0*pow(x_1, 2) - 21439.0L/11520.0L*x_0*x_1 - 128707.0L/115200.0L*x_0 + (1.0L/4800.0L)*pow(x_1, 6) + (29.0L/4800.0L)*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) + (289.0L/640.0L)*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) + (209569.0L/76800.0L)*x_1 + 574799.0L/307200.0L;
       return __pp_r357___result;
    }
    static inline O __pp_r358__(const I &x_0, const I &x_1) {
       O __pp_r358___result;
       __pp_r358___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 29.0L/14400.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/160.0L*pow(x_0, 4)*x_1 - 13.0L/3840.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 + (691.0L/17280.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (337.0L/960.0L)*pow(x_0, 2)*x_1 + (4027.0L/15360.0L)*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) - 37.0L/1440.0L*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) - 2803.0L/2880.0L*x_0*pow(x_1, 2) - 11327.0L/5760.0L*x_0*x_1 - 353399.0L/230400.0L*x_0 + (1.0L/4800.0L)*pow(x_1, 6) + (29.0L/4800.0L)*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) + (289.0L/640.0L)*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) + (211999.0L/76800.0L)*x_1 + 38279.0L/19200.0L;
       return __pp_r358___result;
    }
    static inline O __pp_r359__(const I &x_0, const I &x_1) {
       O __pp_r359___result;
       __pp_r359___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 29.0L/14400.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/160.0L*pow(x_0, 4)*x_1 - 13.0L/3840.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 + (691.0L/17280.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (337.0L/960.0L)*pow(x_0, 2)*x_1 + (4027.0L/15360.0L)*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) - 37.0L/1440.0L*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) - 2803.0L/2880.0L*x_0*pow(x_1, 2) - 11327.0L/5760.0L*x_0*x_1 - 353399.0L/230400.0L*x_0 + (1.0L/4800.0L)*pow(x_1, 6) + (29.0L/4800.0L)*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) + (289.0L/640.0L)*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) + (211999.0L/76800.0L)*x_1 + 38279.0L/19200.0L;
       return __pp_r359___result;
    }
    static inline O __pp_r360__(const I &x_0, const I &x_1) {
       O __pp_r360___result;
       __pp_r360___result = -29.0L/43200.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 + (131.0L/14400.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 89.0L/1440.0L*pow(x_0, 4)*x_1 - 1319.0L/11520.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (41.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (227.0L/480.0L)*pow(x_0, 3)*x_1 + (10931.0L/17280.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 4109.0L/2880.0L*pow(x_0, 2)*x_1 - 69839.0L/46080.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (43.0L/1440.0L)*x_0*pow(x_1, 4) + (103.0L/480.0L)*x_0*pow(x_1, 3) + (2317.0L/2880.0L)*x_0*pow(x_1, 2) + (1017.0L/640.0L)*x_0*x_1 + (301961.0L/230400.0L)*x_0 - 11.0L/43200.0L*pow(x_1, 6) - 73.0L/14400.0L*pow(x_1, 5) - 449.0L/11520.0L*pow(x_1, 4) - 2437.0L/17280.0L*pow(x_1, 3) - 10481.0L/46080.0L*pow(x_1, 2) - 19363.0L/230400.0L*x_1 + 16831.0L/172800.0L;
       return __pp_r360___result;
    }
    static inline O __pp_r361__(const I &x_0, const I &x_1) {
       O __pp_r361___result;
       __pp_r361___result = (11.0L/4800.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 - 89.0L/2880.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/288.0L)*pow(x_0, 4)*x_1 + (1273.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 3)*x_1 - 733.0L/17280.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/36.0L*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 1193.0L/2880.0L*pow(x_0, 2)*x_1 - 17351.0L/46080.0L*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (5.0L/288.0L)*x_0*pow(x_1, 4) + (49.0L/480.0L)*x_0*pow(x_1, 3) + (859.0L/2880.0L)*x_0*pow(x_1, 2) + (9.0L/20.0L)*x_0*x_1 + (13153.0L/46080.0L)*x_0 - 1.0L/4800.0L*pow(x_1, 6) - 11.0L/2880.0L*pow(x_1, 5) - 287.0L/11520.0L*pow(x_1, 4) - 979.0L/17280.0L*pow(x_1, 3) + (2641.0L/46080.0L)*pow(x_1, 2) + (19747.0L/46080.0L)*x_1 + 666089.0L/1382400.0L;
       return __pp_r361___result;
    }
    static inline O __pp_r362__(const I &x_0, const I &x_1) {
       O __pp_r362___result;
       __pp_r362___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 29.0L/14400.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/160.0L*pow(x_0, 4)*x_1 - 13.0L/3840.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 + (691.0L/17280.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (337.0L/960.0L)*pow(x_0, 2)*x_1 + (4027.0L/15360.0L)*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) - 37.0L/1440.0L*x_0*pow(x_1, 4) - 331.0L/1440.0L*x_0*pow(x_1, 3) - 2803.0L/2880.0L*x_0*pow(x_1, 2) - 11327.0L/5760.0L*x_0*x_1 - 353399.0L/230400.0L*x_0 + (1.0L/4800.0L)*pow(x_1, 6) + (29.0L/4800.0L)*pow(x_1, 5) + (277.0L/3840.0L)*pow(x_1, 4) + (289.0L/640.0L)*pow(x_1, 3) + (23813.0L/15360.0L)*pow(x_1, 2) + (211999.0L/76800.0L)*x_1 + 38279.0L/19200.0L;
       return __pp_r362___result;
    }
    static inline O __pp_r363__(const I &x_0, const I &x_1) {
       O __pp_r363___result;
       __pp_r363___result = -29.0L/43200.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 + (131.0L/14400.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 89.0L/1440.0L*pow(x_0, 4)*x_1 - 1319.0L/11520.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (41.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (227.0L/480.0L)*pow(x_0, 3)*x_1 + (10931.0L/17280.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 4109.0L/2880.0L*pow(x_0, 2)*x_1 - 69839.0L/46080.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (43.0L/1440.0L)*x_0*pow(x_1, 4) + (103.0L/480.0L)*x_0*pow(x_1, 3) + (2317.0L/2880.0L)*x_0*pow(x_1, 2) + (1017.0L/640.0L)*x_0*x_1 + (301961.0L/230400.0L)*x_0 - 11.0L/43200.0L*pow(x_1, 6) - 73.0L/14400.0L*pow(x_1, 5) - 449.0L/11520.0L*pow(x_1, 4) - 2437.0L/17280.0L*pow(x_1, 3) - 10481.0L/46080.0L*pow(x_1, 2) - 19363.0L/230400.0L*x_1 + 16831.0L/172800.0L;
       return __pp_r363___result;
    }
    static inline O __pp_r364__(const I &x_0, const I &x_1) {
       O __pp_r364___result;
       __pp_r364___result = (11.0L/4800.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 - 89.0L/2880.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/288.0L)*pow(x_0, 4)*x_1 + (1273.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 3)*x_1 - 733.0L/17280.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/36.0L*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 1193.0L/2880.0L*pow(x_0, 2)*x_1 - 17351.0L/46080.0L*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (5.0L/288.0L)*x_0*pow(x_1, 4) + (49.0L/480.0L)*x_0*pow(x_1, 3) + (859.0L/2880.0L)*x_0*pow(x_1, 2) + (9.0L/20.0L)*x_0*x_1 + (13153.0L/46080.0L)*x_0 - 1.0L/4800.0L*pow(x_1, 6) - 11.0L/2880.0L*pow(x_1, 5) - 287.0L/11520.0L*pow(x_1, 4) - 979.0L/17280.0L*pow(x_1, 3) + (2641.0L/46080.0L)*pow(x_1, 2) + (19747.0L/46080.0L)*x_1 + 666089.0L/1382400.0L;
       return __pp_r364___result;
    }
    static inline O __pp_r365__(const I &x_0, const I &x_1) {
       O __pp_r365___result;
       __pp_r365___result = -29.0L/43200.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 + (131.0L/14400.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 89.0L/1440.0L*pow(x_0, 4)*x_1 - 1319.0L/11520.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (41.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (227.0L/480.0L)*pow(x_0, 3)*x_1 + (10931.0L/17280.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 4109.0L/2880.0L*pow(x_0, 2)*x_1 - 69839.0L/46080.0L*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) + (59.0L/720.0L)*x_0*pow(x_1, 4) + (19.0L/40.0L)*x_0*pow(x_1, 3) + (131.0L/90.0L)*x_0*pow(x_1, 2) + (9227.0L/3840.0L)*x_0*x_1 + (395711.0L/230400.0L)*x_0 + (17.0L/21600.0L)*pow(x_1, 6) + (1.0L/7200.0L)*pow(x_1, 5) - 103.0L/1440.0L*pow(x_1, 4) - 4031.0L/8640.0L*pow(x_1, 3) - 14339.0L/11520.0L*pow(x_1, 2) - 10859.0L/7200.0L*x_1 - 1840079.0L/2764800.0L;
       return __pp_r365___result;
    }
    static inline O __pp_r366__(const I &x_0, const I &x_1) {
       O __pp_r366___result;
       __pp_r366___result = -29.0L/43200.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 + (131.0L/14400.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 89.0L/1440.0L*pow(x_0, 4)*x_1 - 1319.0L/11520.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (41.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (227.0L/480.0L)*pow(x_0, 3)*x_1 + (10931.0L/17280.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 4109.0L/2880.0L*pow(x_0, 2)*x_1 - 69839.0L/46080.0L*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) + (59.0L/720.0L)*x_0*pow(x_1, 4) + (19.0L/40.0L)*x_0*pow(x_1, 3) + (131.0L/90.0L)*x_0*pow(x_1, 2) + (9227.0L/3840.0L)*x_0*x_1 + (395711.0L/230400.0L)*x_0 + (17.0L/21600.0L)*pow(x_1, 6) + (1.0L/7200.0L)*pow(x_1, 5) - 103.0L/1440.0L*pow(x_1, 4) - 4031.0L/8640.0L*pow(x_1, 3) - 14339.0L/11520.0L*pow(x_1, 2) - 10859.0L/7200.0L*x_1 - 1840079.0L/2764800.0L;
       return __pp_r366___result;
    }
    static inline O __pp_r367__(const I &x_0, const I &x_1) {
       O __pp_r367___result;
       __pp_r367___result = (11.0L/4800.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 - 89.0L/2880.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/288.0L)*pow(x_0, 4)*x_1 + (1273.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 3)*x_1 - 733.0L/17280.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/36.0L*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 1193.0L/2880.0L*pow(x_0, 2)*x_1 - 17351.0L/46080.0L*pow(x_0, 2) + (19.0L/3600.0L)*x_0*pow(x_1, 5) + (5.0L/72.0L)*x_0*pow(x_1, 4) + (29.0L/80.0L)*x_0*pow(x_1, 3) + (1367.0L/1440.0L)*x_0*pow(x_1, 2) + (4853.0L/3840.0L)*x_0*x_1 + (31903.0L/46080.0L)*x_0 + (1.0L/1200.0L)*pow(x_1, 6) + (1.0L/720.0L)*pow(x_1, 5) - 331.0L/5760.0L*pow(x_1, 4) - 1651.0L/4320.0L*pow(x_1, 3) - 22117.0L/23040.0L*pow(x_1, 2) - 22939.0L/23040.0L*x_1 - 777197.0L/2764800.0L;
       return __pp_r367___result;
    }
    static inline O __pp_r368__(const I &x_0, const I &x_1) {
       O __pp_r368___result;
       __pp_r368___result = (11.0L/4800.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 - 89.0L/2880.0L*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/288.0L)*pow(x_0, 4)*x_1 + (1273.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 3)*x_1 - 733.0L/17280.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/36.0L*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 1193.0L/2880.0L*pow(x_0, 2)*x_1 - 17351.0L/46080.0L*pow(x_0, 2) + (19.0L/3600.0L)*x_0*pow(x_1, 5) + (5.0L/72.0L)*x_0*pow(x_1, 4) + (29.0L/80.0L)*x_0*pow(x_1, 3) + (1367.0L/1440.0L)*x_0*pow(x_1, 2) + (4853.0L/3840.0L)*x_0*x_1 + (31903.0L/46080.0L)*x_0 + (1.0L/1200.0L)*pow(x_1, 6) + (1.0L/720.0L)*pow(x_1, 5) - 331.0L/5760.0L*pow(x_1, 4) - 1651.0L/4320.0L*pow(x_1, 3) - 22117.0L/23040.0L*pow(x_1, 2) - 22939.0L/23040.0L*x_1 - 777197.0L/2764800.0L;
       return __pp_r368___result;
    }
    static inline O __pp_r369__(const I &x_0, const I &x_1) {
       O __pp_r369___result;
       __pp_r369___result = (43.0L/14400.0L)*pow(x_0, 6) - 41.0L/3600.0L*pow(x_0, 5)*x_1 - 25.0L/576.0L*pow(x_0, 5) + (1.0L/80.0L)*pow(x_0, 4)*pow(x_1, 2) + (29.0L/288.0L)*pow(x_0, 4)*x_1 + (2353.0L/11520.0L)*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/9.0L*pow(x_0, 3)*pow(x_1, 2) - 169.0L/480.0L*pow(x_0, 3)*x_1 - 7213.0L/17280.0L*pow(x_0, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 4) + (7.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (193.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (2047.0L/2880.0L)*pow(x_0, 2)*x_1 + (21529.0L/46080.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (1.0L/144.0L)*x_0*pow(x_1, 4) - 1.0L/80.0L*x_0*pow(x_1, 3) - 253.0L/1440.0L*x_0*pow(x_1, 2) - 1627.0L/3840.0L*x_0*x_1 - 14753.0L/46080.0L*x_0 + (11.0L/7200.0L)*pow(x_1, 6) + (1.0L/72.0L)*pow(x_1, 5) + (209.0L/5760.0L)*pow(x_1, 4) - 31.0L/4320.0L*pow(x_1, 3) - 2677.0L/23040.0L*pow(x_1, 2) + (389.0L/23040.0L)*x_1 + 622483.0L/2764800.0L;
       return __pp_r369___result;
    }
    static inline O __pp_r370__(const I &x_0, const I &x_1) {
       O __pp_r370___result;
       __pp_r370___result = (43.0L/14400.0L)*pow(x_0, 6) - 41.0L/3600.0L*pow(x_0, 5)*x_1 - 25.0L/576.0L*pow(x_0, 5) + (1.0L/80.0L)*pow(x_0, 4)*pow(x_1, 2) + (29.0L/288.0L)*pow(x_0, 4)*x_1 + (2353.0L/11520.0L)*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/9.0L*pow(x_0, 3)*pow(x_1, 2) - 169.0L/480.0L*pow(x_0, 3)*x_1 - 7213.0L/17280.0L*pow(x_0, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 4) + (7.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (193.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (2047.0L/2880.0L)*pow(x_0, 2)*x_1 + (21529.0L/46080.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (1.0L/144.0L)*x_0*pow(x_1, 4) - 1.0L/80.0L*x_0*pow(x_1, 3) - 253.0L/1440.0L*x_0*pow(x_1, 2) - 1627.0L/3840.0L*x_0*x_1 - 14753.0L/46080.0L*x_0 + (11.0L/7200.0L)*pow(x_1, 6) + (1.0L/72.0L)*pow(x_1, 5) + (209.0L/5760.0L)*pow(x_1, 4) - 31.0L/4320.0L*pow(x_1, 3) - 2677.0L/23040.0L*pow(x_1, 2) + (389.0L/23040.0L)*x_1 + 622483.0L/2764800.0L;
       return __pp_r370___result;
    }
    static inline O __pp_r371__(const I &x_0, const I &x_1) {
       O __pp_r371___result;
       __pp_r371___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 - 59.0L/4800.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 4)*x_1 + (157.0L/2304.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) - 23.0L/288.0L*pow(x_0, 3)*x_1 - 115.0L/1152.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (7.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (15.0L/64.0L)*pow(x_0, 2)*x_1 + (2321.0L/46080.0L)*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (1.0L/60.0L)*x_0*pow(x_1, 4) + (1.0L/18.0L)*x_0*pow(x_1, 3) + (1.0L/16.0L)*x_0*pow(x_1, 2) - 79.0L/11520.0L*x_0*x_1 - 2179.0L/76800.0L*x_0 + (1.0L/675.0L)*pow(x_1, 6) + (31.0L/2400.0L)*pow(x_1, 5) + (1.0L/36.0L)*pow(x_1, 4) - 3.0L/64.0L*pow(x_1, 3) - 2539.0L/11520.0L*pow(x_1, 2) - 2477.0L/19200.0L*x_1 + 77437.0L/552960.0L;
       return __pp_r371___result;
    }
    static inline O __pp_r372__(const I &x_0, const I &x_1) {
       O __pp_r372___result;
       __pp_r372___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 - 59.0L/4800.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 4)*x_1 + (157.0L/2304.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) - 23.0L/288.0L*pow(x_0, 3)*x_1 - 115.0L/1152.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (7.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (15.0L/64.0L)*pow(x_0, 2)*x_1 + (2321.0L/46080.0L)*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (1.0L/60.0L)*x_0*pow(x_1, 4) + (1.0L/18.0L)*x_0*pow(x_1, 3) + (1.0L/16.0L)*x_0*pow(x_1, 2) - 79.0L/11520.0L*x_0*x_1 - 2179.0L/76800.0L*x_0 + (1.0L/675.0L)*pow(x_1, 6) + (31.0L/2400.0L)*pow(x_1, 5) + (1.0L/36.0L)*pow(x_1, 4) - 3.0L/64.0L*pow(x_1, 3) - 2539.0L/11520.0L*pow(x_1, 2) - 2477.0L/19200.0L*x_1 + 77437.0L/552960.0L;
       return __pp_r372___result;
    }
    static inline O __pp_r373__(const I &x_0, const I &x_1) {
       O __pp_r373___result;
       __pp_r373___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 - 19.0L/4800.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 4)*x_1 + (109.0L/2304.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) - 23.0L/288.0L*pow(x_0, 3)*x_1 - 91.0L/1152.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (7.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (15.0L/64.0L)*pow(x_0, 2)*x_1 + (1841.0L/46080.0L)*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (1.0L/60.0L)*x_0*pow(x_1, 4) + (1.0L/18.0L)*x_0*pow(x_1, 3) + (1.0L/16.0L)*x_0*pow(x_1, 2) - 79.0L/11520.0L*x_0*x_1 - 1979.0L/76800.0L*x_0 + (1.0L/675.0L)*pow(x_1, 6) + (31.0L/2400.0L)*pow(x_1, 5) + (1.0L/36.0L)*pow(x_1, 4) - 3.0L/64.0L*pow(x_1, 3) - 2539.0L/11520.0L*pow(x_1, 2) - 2477.0L/19200.0L*x_1 + 77293.0L/552960.0L;
       return __pp_r373___result;
    }
    static inline O __pp_r374__(const I &x_0, const I &x_1) {
       O __pp_r374___result;
       __pp_r374___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 - 19.0L/4800.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 4)*x_1 + (109.0L/2304.0L)*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) - 23.0L/288.0L*pow(x_0, 3)*x_1 - 91.0L/1152.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (7.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (15.0L/64.0L)*pow(x_0, 2)*x_1 + (1841.0L/46080.0L)*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (1.0L/60.0L)*x_0*pow(x_1, 4) + (1.0L/18.0L)*x_0*pow(x_1, 3) + (1.0L/16.0L)*x_0*pow(x_1, 2) - 79.0L/11520.0L*x_0*x_1 - 1979.0L/76800.0L*x_0 + (1.0L/675.0L)*pow(x_1, 6) + (31.0L/2400.0L)*pow(x_1, 5) + (1.0L/36.0L)*pow(x_1, 4) - 3.0L/64.0L*pow(x_1, 3) - 2539.0L/11520.0L*pow(x_1, 2) - 2477.0L/19200.0L*x_1 + 77293.0L/552960.0L;
       return __pp_r374___result;
    }
    static inline O __pp_r375__(const I &x_0, const I &x_1) {
       O __pp_r375___result;
       __pp_r375___result = -7.0L/4800.0L*pow(x_0, 6) + (7.0L/3600.0L)*pow(x_0, 5)*x_1 + (103.0L/14400.0L)*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) - 7.0L/1440.0L*pow(x_0, 4)*x_1 + (29.0L/2304.0L)*pow(x_0, 4) - 1.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 3)*x_1 - 73.0L/3456.0L*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) + (2.0L/45.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (85.0L/576.0L)*pow(x_0, 2)*x_1 - 659.0L/46080.0L*pow(x_0, 2) + (7.0L/3600.0L)*x_0*pow(x_1, 5) + (29.0L/1440.0L)*x_0*pow(x_1, 4) + (7.0L/96.0L)*x_0*pow(x_1, 3) + (61.0L/576.0L)*x_0*pow(x_1, 2) + (91.0L/1920.0L)*x_0*x_1 + (313.0L/230400.0L)*x_0 + (7.0L/4800.0L)*pow(x_1, 6) + (181.0L/14400.0L)*pow(x_1, 5) + (59.0L/2304.0L)*pow(x_1, 4) - 187.0L/3456.0L*pow(x_1, 3) - 10781.0L/46080.0L*pow(x_1, 2) - 32849.0L/230400.0L*x_1 + 9271.0L/69120.0L;
       return __pp_r375___result;
    }
    static inline O __pp_r376__(const I &x_0, const I &x_1) {
       O __pp_r376___result;
       __pp_r376___result = -7.0L/4800.0L*pow(x_0, 6) + (7.0L/3600.0L)*pow(x_0, 5)*x_1 + (103.0L/14400.0L)*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) - 7.0L/1440.0L*pow(x_0, 4)*x_1 + (29.0L/2304.0L)*pow(x_0, 4) - 1.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 3)*x_1 - 73.0L/3456.0L*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) + (2.0L/45.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (85.0L/576.0L)*pow(x_0, 2)*x_1 - 659.0L/46080.0L*pow(x_0, 2) + (7.0L/3600.0L)*x_0*pow(x_1, 5) + (29.0L/1440.0L)*x_0*pow(x_1, 4) + (7.0L/96.0L)*x_0*pow(x_1, 3) + (61.0L/576.0L)*x_0*pow(x_1, 2) + (91.0L/1920.0L)*x_0*x_1 + (313.0L/230400.0L)*x_0 + (7.0L/4800.0L)*pow(x_1, 6) + (181.0L/14400.0L)*pow(x_1, 5) + (59.0L/2304.0L)*pow(x_1, 4) - 187.0L/3456.0L*pow(x_1, 3) - 10781.0L/46080.0L*pow(x_1, 2) - 32849.0L/230400.0L*x_1 + 9271.0L/69120.0L;
       return __pp_r376___result;
    }
    static inline O __pp_r377__(const I &x_0, const I &x_1) {
       O __pp_r377___result;
       __pp_r377___result = -73.0L/43200.0L*pow(x_0, 6) + (1.0L/300.0L)*pow(x_0, 5)*x_1 + (143.0L/14400.0L)*pow(x_0, 5) - 11.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 3.0L/160.0L*pow(x_0, 4)*x_1 - 1.0L/768.0L*pow(x_0, 4) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/45.0L)*pow(x_0, 3)*pow(x_1, 2) + (13.0L/288.0L)*pow(x_0, 3)*x_1 + (55.0L/3456.0L)*pow(x_0, 3) + (1.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (7.0L/192.0L)*pow(x_0, 2)*x_1 - 1073.0L/15360.0L*pow(x_0, 2) + (1.0L/300.0L)*x_0*pow(x_1, 5) + (49.0L/1440.0L)*x_0*pow(x_1, 4) + (37.0L/288.0L)*x_0*pow(x_1, 3) + (125.0L/576.0L)*x_0*pow(x_1, 2) + (913.0L/5760.0L)*x_0*x_1 + (10553.0L/230400.0L)*x_0 + (53.0L/43200.0L)*pow(x_1, 6) + (47.0L/4800.0L)*pow(x_1, 5) + (3.0L/256.0L)*pow(x_1, 4) - 35.0L/384.0L*pow(x_1, 3) - 4447.0L/15360.0L*pow(x_1, 2) - 14363.0L/76800.0L*x_1 + 2749.0L/23040.0L;
       return __pp_r377___result;
    }
    static inline O __pp_r378__(const I &x_0, const I &x_1) {
       O __pp_r378___result;
       __pp_r378___result = -73.0L/43200.0L*pow(x_0, 6) + (1.0L/300.0L)*pow(x_0, 5)*x_1 + (143.0L/14400.0L)*pow(x_0, 5) - 11.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 3.0L/160.0L*pow(x_0, 4)*x_1 - 1.0L/768.0L*pow(x_0, 4) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/45.0L)*pow(x_0, 3)*pow(x_1, 2) + (13.0L/288.0L)*pow(x_0, 3)*x_1 + (55.0L/3456.0L)*pow(x_0, 3) + (1.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (7.0L/192.0L)*pow(x_0, 2)*x_1 - 1073.0L/15360.0L*pow(x_0, 2) + (1.0L/300.0L)*x_0*pow(x_1, 5) + (49.0L/1440.0L)*x_0*pow(x_1, 4) + (37.0L/288.0L)*x_0*pow(x_1, 3) + (125.0L/576.0L)*x_0*pow(x_1, 2) + (913.0L/5760.0L)*x_0*x_1 + (10553.0L/230400.0L)*x_0 + (53.0L/43200.0L)*pow(x_1, 6) + (47.0L/4800.0L)*pow(x_1, 5) + (3.0L/256.0L)*pow(x_1, 4) - 35.0L/384.0L*pow(x_1, 3) - 4447.0L/15360.0L*pow(x_1, 2) - 14363.0L/76800.0L*x_1 + 2749.0L/23040.0L;
       return __pp_r378___result;
    }
    static inline O __pp_r379__(const I &x_0, const I &x_1) {
       O __pp_r379___result;
       __pp_r379___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (47.0L/14400.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/480.0L*pow(x_0, 4)*x_1 + (43.0L/3840.0L)*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (29.0L/1440.0L)*pow(x_0, 3)*x_1 + (59.0L/17280.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/40.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (53.0L/960.0L)*pow(x_0, 2)*x_1 - 193.0L/3072.0L*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) + (23.0L/720.0L)*x_0*pow(x_1, 4) + (11.0L/90.0L)*x_0*pow(x_1, 3) + (299.0L/1440.0L)*x_0*pow(x_1, 2) + (349.0L/2304.0L)*x_0*x_1 + (10067.0L/230400.0L)*x_0 + (1.0L/800.0L)*pow(x_1, 6) + (1.0L/100.0L)*pow(x_1, 5) + (1.0L/80.0L)*pow(x_1, 4) - 43.0L/480.0L*pow(x_1, 3) - 221.0L/768.0L*pow(x_1, 2) - 7141.0L/38400.0L*x_1 + 110203.0L/921600.0L;
       return __pp_r379___result;
    }
    static inline O __pp_r380__(const I &x_0, const I &x_1) {
       O __pp_r380___result;
       __pp_r380___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (47.0L/14400.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/480.0L*pow(x_0, 4)*x_1 + (43.0L/3840.0L)*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (29.0L/1440.0L)*pow(x_0, 3)*x_1 + (59.0L/17280.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/40.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (53.0L/960.0L)*pow(x_0, 2)*x_1 - 193.0L/3072.0L*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) + (23.0L/720.0L)*x_0*pow(x_1, 4) + (11.0L/90.0L)*x_0*pow(x_1, 3) + (299.0L/1440.0L)*x_0*pow(x_1, 2) + (349.0L/2304.0L)*x_0*x_1 + (10067.0L/230400.0L)*x_0 + (1.0L/800.0L)*pow(x_1, 6) + (1.0L/100.0L)*pow(x_1, 5) + (1.0L/80.0L)*pow(x_1, 4) - 43.0L/480.0L*pow(x_1, 3) - 221.0L/768.0L*pow(x_1, 2) - 7141.0L/38400.0L*x_1 + 110203.0L/921600.0L;
       return __pp_r380___result;
    }
    static inline O __pp_r381__(const I &x_0, const I &x_1) {
       O __pp_r381___result;
       __pp_r381___result = -37.0L/21600.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (23.0L/2400.0L)*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/45.0L*pow(x_0, 4)*x_1 - 1.0L/288.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/36.0L)*pow(x_0, 3)*x_1 + (5.0L/576.0L)*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/144.0L*pow(x_0, 2)*x_1 - 961.0L/11520.0L*pow(x_0, 2) - 19.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/40.0L*x_0*pow(x_1, 4) - 5.0L/144.0L*x_0*pow(x_1, 3) - 1.0L/96.0L*x_0*pow(x_1, 2) - 7.0L/5760.0L*x_0*x_1 + (23.0L/38400.0L)*x_0 - 7.0L/5400.0L*pow(x_1, 6) - 1.0L/3600.0L*pow(x_1, 5) + (23.0L/1152.0L)*pow(x_1, 4) + (13.0L/864.0L)*pow(x_1, 3) - 1643.0L/23040.0L*pow(x_1, 2) + (119.0L/57600.0L)*x_1 + 10003.0L/55296.0L;
       return __pp_r381___result;
    }
    static inline O __pp_r382__(const I &x_0, const I &x_1) {
       O __pp_r382___result;
       __pp_r382___result = -1.0L/4320.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (7.0L/2400.0L)*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/180.0L*pow(x_0, 4)*x_1 + (13.0L/1440.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/360.0L)*pow(x_0, 3)*x_1 - 11.0L/2880.0L*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/240.0L)*pow(x_0, 2)*pow(x_1, 2) + (17.0L/1440.0L)*pow(x_0, 2)*x_1 - 11.0L/144.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) - 13.0L/480.0L*x_0*pow(x_1, 4) - 59.0L/1440.0L*x_0*pow(x_1, 3) - 19.0L/960.0L*x_0*pow(x_1, 2) - 19.0L/2304.0L*x_0*x_1 - 29.0L/19200.0L*x_0 - 11.0L/8640.0L*pow(x_1, 6) - 1.0L/14400.0L*pow(x_1, 5) + (239.0L/11520.0L)*pow(x_1, 4) + (287.0L/17280.0L)*pow(x_1, 3) - 641.0L/9216.0L*pow(x_1, 2) + (719.0L/230400.0L)*x_1 + 500879.0L/2764800.0L;
       return __pp_r382___result;
    }
    static inline O __pp_r383__(const I &x_0, const I &x_1) {
       O __pp_r383___result;
       __pp_r383___result = -1.0L/2160.0L*pow(x_0, 6) + (31.0L/7200.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/80.0L*pow(x_0, 4)*x_1 + (1.0L/180.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 3)*x_1 + (7.0L/8640.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/480.0L*pow(x_0, 2)*x_1 - 23.0L/288.0L*pow(x_0, 2) - 1.0L/240.0L*x_0*pow(x_1, 5) - 29.0L/1440.0L*x_0*pow(x_1, 4) - 13.0L/480.0L*x_0*pow(x_1, 3) - 17.0L/2880.0L*x_0*pow(x_1, 2) - 1.0L/768.0L*x_0*x_1 - 7.0L/57600.0L*x_0 - 13.0L/8640.0L*pow(x_1, 6) - 7.0L/4800.0L*pow(x_1, 5) + (199.0L/11520.0L)*pow(x_1, 4) + (23.0L/1920.0L)*pow(x_1, 3) - 673.0L/9216.0L*pow(x_1, 2) + (133.0L/76800.0L)*x_1 + 500239.0L/2764800.0L;
       return __pp_r383___result;
    }
    static inline O __pp_r384__(const I &x_0, const I &x_1) {
       O __pp_r384___result;
       __pp_r384___result = -73.0L/43200.0L*pow(x_0, 6) + (1.0L/300.0L)*pow(x_0, 5)*x_1 + (143.0L/14400.0L)*pow(x_0, 5) - 11.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 3.0L/160.0L*pow(x_0, 4)*x_1 - 1.0L/768.0L*pow(x_0, 4) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/45.0L)*pow(x_0, 3)*pow(x_1, 2) + (13.0L/288.0L)*pow(x_0, 3)*x_1 + (55.0L/3456.0L)*pow(x_0, 3) + (1.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (7.0L/192.0L)*pow(x_0, 2)*x_1 - 1073.0L/15360.0L*pow(x_0, 2) - 1.0L/1200.0L*x_0*pow(x_1, 5) + (1.0L/360.0L)*x_0*pow(x_1, 4) + (5.0L/144.0L)*x_0*pow(x_1, 3) + (11.0L/144.0L)*x_0*pow(x_1, 2) + (611.0L/11520.0L)*x_0*x_1 + (3263.0L/230400.0L)*x_0 + (1.0L/5400.0L)*pow(x_1, 6) + (13.0L/1200.0L)*pow(x_1, 5) + (7.0L/128.0L)*pow(x_1, 4) + (7.0L/96.0L)*pow(x_1, 3) - 131.0L/7680.0L*pow(x_1, 2) + (1121.0L/38400.0L)*x_1 + 6877.0L/36864.0L;
       return __pp_r384___result;
    }
    static inline O __pp_r385__(const I &x_0, const I &x_1) {
       O __pp_r385___result;
       __pp_r385___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (47.0L/14400.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/480.0L*pow(x_0, 4)*x_1 + (43.0L/3840.0L)*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (29.0L/1440.0L)*pow(x_0, 3)*x_1 + (59.0L/17280.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/40.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (53.0L/960.0L)*pow(x_0, 2)*x_1 - 193.0L/3072.0L*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) + (1.0L/1440.0L)*x_0*pow(x_1, 4) + (41.0L/1440.0L)*x_0*pow(x_1, 3) + (193.0L/2880.0L)*x_0*pow(x_1, 2) + (53.0L/1152.0L)*x_0*x_1 + (2777.0L/230400.0L)*x_0 + (1.0L/4800.0L)*pow(x_1, 6) + (53.0L/4800.0L)*pow(x_1, 5) + (71.0L/1280.0L)*pow(x_1, 4) + (143.0L/1920.0L)*pow(x_1, 3) - 47.0L/3072.0L*pow(x_1, 2) + (2323.0L/76800.0L)*x_1 + 21521.0L/115200.0L;
       return __pp_r385___result;
    }
    static inline O __pp_r386__(const I &x_0, const I &x_1) {
       O __pp_r386___result;
       __pp_r386___result = -19.0L/43200.0L*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (67.0L/14400.0L)*pow(x_0, 5) - 1.0L/180.0L*pow(x_0, 4)*pow(x_1, 2) - 13.0L/1440.0L*pow(x_0, 4)*x_1 + (89.0L/11520.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (49.0L/1440.0L)*pow(x_0, 3)*x_1 + (139.0L/17280.0L)*pow(x_0, 3) - 1.0L/720.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (29.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (119.0L/2880.0L)*pow(x_0, 2)*x_1 - 611.0L/9216.0L*pow(x_0, 2) + (1.0L/3600.0L)*x_0*pow(x_1, 5) + (11.0L/1440.0L)*x_0*pow(x_1, 4) + (61.0L/1440.0L)*x_0*pow(x_1, 3) + (233.0L/2880.0L)*x_0*pow(x_1, 2) + (61.0L/1152.0L)*x_0*x_1 + (3097.0L/230400.0L)*x_0 - 1.0L/43200.0L*pow(x_1, 6) + (139.0L/14400.0L)*pow(x_1, 5) + (599.0L/11520.0L)*pow(x_1, 4) + (1207.0L/17280.0L)*pow(x_1, 3) - 173.0L/9216.0L*pow(x_1, 2) + (6649.0L/230400.0L)*x_1 + 64483.0L/345600.0L;
       return __pp_r386___result;
    }
    static inline O __pp_r387__(const I &x_0, const I &x_1) {
       O __pp_r387___result;
       __pp_r387___result = (1.0L/400.0L)*pow(x_0, 6) - 2.0L/225.0L*pow(x_0, 5)*x_1 - 1.0L/7200.0L*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/720.0L*pow(x_0, 4)*x_1 + (1.0L/120.0L)*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/90.0L)*pow(x_0, 3)*x_1 - 1.0L/8640.0L*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/1440.0L*pow(x_0, 2)*x_1 - 51.0L/640.0L*pow(x_0, 2) - 17.0L/3600.0L*x_0*pow(x_1, 5) - 31.0L/1440.0L*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) - 19.0L/2880.0L*x_0*pow(x_1, 2) - 17.0L/11520.0L*x_0*x_1 - 1.0L/7200.0L*x_0 - 7.0L/4800.0L*pow(x_1, 6) - 19.0L/14400.0L*pow(x_1, 5) + (67.0L/3840.0L)*pow(x_1, 4) + (209.0L/17280.0L)*pow(x_1, 3) - 1121.0L/15360.0L*pow(x_1, 2) + (401.0L/230400.0L)*x_1 + 166747.0L/921600.0L;
       return __pp_r387___result;
    }
    static inline O __pp_r388__(const I &x_0, const I &x_1) {
       O __pp_r388___result;
       __pp_r388___result = (1.0L/2400.0L)*pow(x_0, 6) - 1.0L/1800.0L*pow(x_0, 5)*x_1 + (7.0L/3600.0L)*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) + (7.0L/360.0L)*pow(x_0, 4)*x_1 + (41.0L/1920.0L)*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (23.0L/720.0L)*pow(x_0, 3)*x_1 + (67.0L/4320.0L)*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) + (7.0L/720.0L)*pow(x_0, 2)*x_1 - 547.0L/7680.0L*pow(x_0, 2) - 17.0L/3600.0L*x_0*pow(x_1, 5) - 31.0L/1440.0L*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) - 19.0L/2880.0L*x_0*pow(x_1, 2) + (13.0L/11520.0L)*x_0*x_1 + (239.0L/115200.0L)*x_0 - 7.0L/4800.0L*pow(x_1, 6) - 19.0L/14400.0L*pow(x_1, 5) + (67.0L/3840.0L)*pow(x_1, 4) + (209.0L/17280.0L)*pow(x_1, 3) - 1121.0L/15360.0L*pow(x_1, 2) + (461.0L/230400.0L)*x_1 + 166957.0L/921600.0L;
       return __pp_r388___result;
    }
    static inline O __pp_r389__(const I &x_0, const I &x_1) {
       O __pp_r389___result;
       __pp_r389___result = -23.0L/21600.0L*pow(x_0, 6) + (7.0L/1800.0L)*pow(x_0, 5)*x_1 - 1.0L/3600.0L*pow(x_0, 5) - 1.0L/720.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/40.0L)*pow(x_0, 4)*x_1 + (23.0L/1152.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 3)*x_1 + (13.0L/864.0L)*pow(x_0, 3) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 2)*x_1 - 1643.0L/23040.0L*pow(x_0, 2) - 1.0L/225.0L*x_0*pow(x_1, 5) - 1.0L/45.0L*x_0*pow(x_1, 4) - 1.0L/36.0L*x_0*pow(x_1, 3) - 1.0L/144.0L*x_0*pow(x_1, 2) + (7.0L/5760.0L)*x_0*x_1 + (119.0L/57600.0L)*x_0 - 1.0L/675.0L*pow(x_1, 6) - 1.0L/800.0L*pow(x_1, 5) + (5.0L/288.0L)*pow(x_1, 4) + (7.0L/576.0L)*pow(x_1, 3) - 841.0L/11520.0L*pow(x_1, 2) + (77.0L/38400.0L)*x_1 + 50087.0L/276480.0L;
       return __pp_r389___result;
    }
    static inline O __pp_r390__(const I &x_0, const I &x_1) {
       O __pp_r390___result;
       __pp_r390___result = (109.0L/43200.0L)*pow(x_0, 6) - 31.0L/3600.0L*pow(x_0, 5)*x_1 + (1.0L/4800.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/480.0L)*pow(x_0, 4)*x_1 + (121.0L/11520.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 + (41.0L/5760.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) + (31.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (41.0L/960.0L)*pow(x_0, 2)*x_1 - 3047.0L/46080.0L*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/160.0L)*x_0*pow(x_1, 4) + (59.0L/1440.0L)*x_0*pow(x_1, 3) + (77.0L/960.0L)*x_0*pow(x_1, 2) + (19.0L/360.0L)*x_0*x_1 + (1031.0L/76800.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (47.0L/4800.0L)*pow(x_1, 5) + (601.0L/11520.0L)*pow(x_1, 4) + (403.0L/5760.0L)*pow(x_1, 3) - 863.0L/46080.0L*pow(x_1, 2) + (739.0L/25600.0L)*x_1 + 257933.0L/1382400.0L;
       return __pp_r390___result;
    }
    static inline O __pp_r391__(const I &x_0, const I &x_1) {
       O __pp_r391___result;
       __pp_r391___result = (19.0L/43200.0L)*pow(x_0, 6) - 1.0L/3600.0L*pow(x_0, 5)*x_1 + (11.0L/4800.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 4)*x_1 + (271.0L/11520.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (71.0L/1440.0L)*pow(x_0, 3)*x_1 + (131.0L/5760.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/60.0L)*pow(x_0, 2)*pow(x_1, 3) + (31.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (17.0L/320.0L)*pow(x_0, 2)*x_1 - 2657.0L/46080.0L*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/160.0L)*x_0*pow(x_1, 4) + (59.0L/1440.0L)*x_0*pow(x_1, 3) + (77.0L/960.0L)*x_0*pow(x_1, 2) + (319.0L/5760.0L)*x_0*x_1 + (1201.0L/76800.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (47.0L/4800.0L)*pow(x_1, 5) + (601.0L/11520.0L)*pow(x_1, 4) + (403.0L/5760.0L)*pow(x_1, 3) - 863.0L/46080.0L*pow(x_1, 2) + (2237.0L/76800.0L)*x_1 + 32281.0L/172800.0L;
       return __pp_r391___result;
    }
    static inline O __pp_r392__(const I &x_0, const I &x_1) {
       O __pp_r392___result;
       __pp_r392___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 + (1.0L/14400.0L)*pow(x_0, 5) + (41.0L/1440.0L)*pow(x_0, 4)*x_1 + (17.0L/768.0L)*pow(x_0, 4) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/96.0L)*pow(x_0, 3)*x_1 + (77.0L/3456.0L)*pow(x_0, 3) + (7.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (31.0L/576.0L)*pow(x_0, 2)*x_1 - 887.0L/15360.0L*pow(x_0, 2) + (1.0L/180.0L)*x_0*pow(x_1, 4) + (1.0L/24.0L)*x_0*pow(x_1, 3) + (23.0L/288.0L)*x_0*pow(x_1, 2) + (71.0L/1280.0L)*x_0*x_1 + (3601.0L/230400.0L)*x_0 + (71.0L/7200.0L)*pow(x_1, 5) + (5.0L/96.0L)*pow(x_1, 4) + (121.0L/1728.0L)*pow(x_1, 3) - 3.0L/160.0L*pow(x_1, 2) + (839.0L/28800.0L)*x_1 + 34433.0L/184320.0L;
       return __pp_r392___result;
    }
    static inline O __pp_r393__(const I &x_0, const I &x_1) {
       O __pp_r393___result;
       __pp_r393___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (11.0L/14400.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/1440.0L)*pow(x_0, 4)*x_1 + (121.0L/11520.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (11.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (121.0L/1440.0L)*pow(x_0, 3)*x_1 + (1331.0L/17280.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (11.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (121.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (1331.0L/2880.0L)*pow(x_0, 2)*x_1 + (14641.0L/46080.0L)*pow(x_0, 2) + (1.0L/3600.0L)*x_0*pow(x_1, 5) - 17.0L/1440.0L*x_0*pow(x_1, 4) - 251.0L/1440.0L*x_0*pow(x_1, 3) - 2483.0L/2880.0L*x_0*pow(x_1, 2) - 10687.0L/5760.0L*x_0*x_1 - 343159.0L/230400.0L*x_0 + (19.0L/43200.0L)*pow(x_1, 6) + (127.0L/14400.0L)*pow(x_1, 5) + (991.0L/11520.0L)*pow(x_1, 4) + (8443.0L/17280.0L)*pow(x_1, 3) + (73999.0L/46080.0L)*pow(x_1, 2) + (646237.0L/230400.0L)*x_1 + 347071.0L/172800.0L;
       return __pp_r393___result;
    }
    static inline O __pp_r394__(const I &x_0, const I &x_1) {
       O __pp_r394___result;
       __pp_r394___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (11.0L/14400.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/1440.0L)*pow(x_0, 4)*x_1 + (121.0L/11520.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (11.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (121.0L/1440.0L)*pow(x_0, 3)*x_1 + (1331.0L/17280.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (11.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (121.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (1331.0L/2880.0L)*pow(x_0, 2)*x_1 + (14641.0L/46080.0L)*pow(x_0, 2) + (1.0L/3600.0L)*x_0*pow(x_1, 5) - 17.0L/1440.0L*x_0*pow(x_1, 4) - 251.0L/1440.0L*x_0*pow(x_1, 3) - 2483.0L/2880.0L*x_0*pow(x_1, 2) - 10687.0L/5760.0L*x_0*x_1 - 343159.0L/230400.0L*x_0 + (19.0L/43200.0L)*pow(x_1, 6) + (127.0L/14400.0L)*pow(x_1, 5) + (991.0L/11520.0L)*pow(x_1, 4) + (8443.0L/17280.0L)*pow(x_1, 3) + (73999.0L/46080.0L)*pow(x_1, 2) + (646237.0L/230400.0L)*x_1 + 347071.0L/172800.0L;
       return __pp_r394___result;
    }
    static inline O __pp_r395__(const I &x_0, const I &x_1) {
       O __pp_r395___result;
       __pp_r395___result = -19.0L/43200.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (19.0L/1600.0L)*pow(x_0, 5) - 1.0L/180.0L*pow(x_0, 4)*pow(x_1, 2) - 23.0L/480.0L*pow(x_0, 4)*x_1 - 1159.0L/11520.0L*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) + (17.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (761.0L/1440.0L)*pow(x_0, 3)*x_1 + (3857.0L/5760.0L)*pow(x_0, 3) - 1.0L/720.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/20.0L*pow(x_0, 2)*pow(x_1, 3) - 199.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 421.0L/320.0L*pow(x_0, 2)*x_1 - 67279.0L/46080.0L*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) + (7.0L/160.0L)*x_0*pow(x_1, 4) + (389.0L/1440.0L)*x_0*pow(x_1, 3) + (293.0L/320.0L)*x_0*pow(x_1, 2) + (9793.0L/5760.0L)*x_0*x_1 + (34689.0L/25600.0L)*x_0 - 1.0L/43200.0L*pow(x_1, 6) - 11.0L/4800.0L*pow(x_1, 5) - 289.0L/11520.0L*pow(x_1, 4) - 599.0L/5760.0L*pow(x_1, 3) - 7921.0L/46080.0L*pow(x_1, 2) - 3041.0L/76800.0L*x_1 + 19391.0L/172800.0L;
       return __pp_r395___result;
    }
    static inline O __pp_r396__(const I &x_0, const I &x_1) {
       O __pp_r396___result;
       __pp_r396___result = -1.0L/240.0L*x_0*pow(x_1, 5) - 7.0L/96.0L*x_0*pow(x_1, 4) - 49.0L/96.0L*x_0*pow(x_1, 3) - 343.0L/192.0L*x_0*pow(x_1, 2) - 2401.0L/768.0L*x_0*x_1 - 16807.0L/7680.0L*x_0 - 1.0L/960.0L*pow(x_1, 6) - 1.0L/64.0L*pow(x_1, 5) - 21.0L/256.0L*pow(x_1, 4) - 49.0L/384.0L*pow(x_1, 3) + (343.0L/1024.0L)*pow(x_1, 2) + (7203.0L/5120.0L)*x_1 + 16807.0L/12288.0L;
       return __pp_r396___result;
    }
    static inline O __pp_r397__(const I &x_0, const I &x_1) {
       O __pp_r397___result;
       __pp_r397___result = -1.0L/240.0L*x_0*pow(x_1, 5) - 7.0L/96.0L*x_0*pow(x_1, 4) - 49.0L/96.0L*x_0*pow(x_1, 3) - 343.0L/192.0L*x_0*pow(x_1, 2) - 2401.0L/768.0L*x_0*x_1 - 16807.0L/7680.0L*x_0 - 1.0L/960.0L*pow(x_1, 6) - 1.0L/64.0L*pow(x_1, 5) - 21.0L/256.0L*pow(x_1, 4) - 49.0L/384.0L*pow(x_1, 3) + (343.0L/1024.0L)*pow(x_1, 2) + (7203.0L/5120.0L)*x_1 + 16807.0L/12288.0L;
       return __pp_r397___result;
    }
    static inline O __pp_r398__(const I &x_0, const I &x_1) {
       O __pp_r398___result;
       __pp_r398___result = -1.0L/2160.0L*pow(x_0, 6) + (1.0L/360.0L)*pow(x_0, 5)*x_1 + (1.0L/90.0L)*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 4)*x_1 - 1.0L/9.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/9.0L)*pow(x_0, 3)*pow(x_1, 2) + (4.0L/9.0L)*pow(x_0, 3)*x_1 + (16.0L/27.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/9.0L*pow(x_0, 2)*pow(x_1, 3) - 2.0L/3.0L*pow(x_0, 2)*pow(x_1, 2) - 16.0L/9.0L*pow(x_0, 2)*x_1 - 16.0L/9.0L*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) - 5.0L/288.0L*x_0*pow(x_1, 4) - 19.0L/288.0L*x_0*pow(x_1, 3) - 5.0L/576.0L*x_0*pow(x_1, 2) + (989.0L/2304.0L)*x_0*x_1 + (3023.0L/4608.0L)*x_0 - 13.0L/8640.0L*pow(x_1, 6) - 77.0L/2880.0L*pow(x_1, 5) - 445.0L/2304.0L*pow(x_1, 4) - 2489.0L/3456.0L*pow(x_1, 3) - 13297.0L/9216.0L*pow(x_1, 2) - 13249.0L/9216.0L*x_1 - 292261.0L/552960.0L;
       return __pp_r398___result;
    }
    static inline O __pp_r399__(const I &x_0, const I &x_1) {
       O __pp_r399___result;
       __pp_r399___result = (109.0L/43200.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 - 9.0L/320.0L*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/96.0L)*pow(x_0, 4)*x_1 + (1433.0L/11520.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (113.0L/1440.0L)*pow(x_0, 3)*x_1 - 31.0L/5760.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 37.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 97.0L/320.0L*pow(x_0, 2)*x_1 - 14791.0L/46080.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) + (1.0L/32.0L)*x_0*pow(x_1, 4) + (227.0L/1440.0L)*x_0*pow(x_1, 3) + (131.0L/320.0L)*x_0*pow(x_1, 2) + (101.0L/180.0L)*x_0*x_1 + (1689.0L/5120.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 1.0L/960.0L*pow(x_1, 5) - 127.0L/11520.0L*pow(x_1, 4) - 113.0L/5760.0L*pow(x_1, 3) + (5201.0L/46080.0L)*pow(x_1, 2) + (1453.0L/3072.0L)*x_1 + 686569.0L/1382400.0L;
       return __pp_r399___result;
    }
    static inline O __pp_r400__(const I &x_0, const I &x_1) {
       O __pp_r400___result;
       __pp_r400___result = (1.0L/400.0L)*pow(x_0, 6) - 11.0L/1800.0L*pow(x_0, 5)*x_1 - 13.0L/450.0L*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) + (2.0L/45.0L)*pow(x_0, 4)*x_1 + (41.0L/360.0L)*pow(x_0, 4) + (1.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) - 1.0L/180.0L*pow(x_0, 3)*x_1 - 89.0L/1080.0L*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) - 11.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 79.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) - 551.0L/720.0L*pow(x_0, 2)*x_1 - 3679.0L/5760.0L*pow(x_0, 2) - 7.0L/3600.0L*x_0*pow(x_1, 5) - 43.0L/1440.0L*x_0*pow(x_1, 4) - 257.0L/1440.0L*x_0*pow(x_1, 3) - 1483.0L/2880.0L*x_0*pow(x_1, 2) - 8177.0L/11520.0L*x_0*x_1 - 42523.0L/115200.0L*x_0 - 7.0L/4800.0L*pow(x_1, 6) - 367.0L/14400.0L*pow(x_1, 5) - 2063.0L/11520.0L*pow(x_1, 4) - 10987.0L/17280.0L*pow(x_1, 3) - 53363.0L/46080.0L*pow(x_1, 2) - 213127.0L/230400.0L*x_1 - 398423.0L/2764800.0L;
       return __pp_r400___result;
    }
    static inline O __pp_r401__(const I &x_0, const I &x_1) {
       O __pp_r401___result;
       __pp_r401___result = (19.0L/43200.0L)*pow(x_0, 6) + (1.0L/400.0L)*pow(x_0, 5)*x_1 + (7.0L/960.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/32.0L)*pow(x_0, 4)*x_1 + (503.0L/11520.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (143.0L/1440.0L)*pow(x_0, 3)*x_1 + (419.0L/5760.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 37.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 301.0L/960.0L*pow(x_0, 2)*x_1 - 16561.0L/46080.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) + (1.0L/32.0L)*x_0*pow(x_1, 4) + (227.0L/1440.0L)*x_0*pow(x_1, 3) + (131.0L/320.0L)*x_0*pow(x_1, 2) + (3247.0L/5760.0L)*x_0*x_1 + (5213.0L/15360.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 1.0L/960.0L*pow(x_1, 5) - 127.0L/11520.0L*pow(x_1, 4) - 113.0L/5760.0L*pow(x_1, 3) + (5201.0L/46080.0L)*pow(x_1, 2) + (7261.0L/15360.0L)*x_1 + 42829.0L/86400.0L;
       return __pp_r401___result;
    }
    static inline O __pp_r402__(const I &x_0, const I &x_1) {
       O __pp_r402___result;
       __pp_r402___result = (1.0L/2400.0L)*pow(x_0, 6) + (1.0L/450.0L)*pow(x_0, 5)*x_1 + (47.0L/7200.0L)*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) + (17.0L/720.0L)*pow(x_0, 4)*x_1 + (191.0L/5760.0L)*pow(x_0, 4) + (1.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (11.0L/720.0L)*pow(x_0, 3)*x_1 - 37.0L/8640.0L*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) - 11.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 79.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) - 1117.0L/1440.0L*pow(x_0, 2)*x_1 - 15601.0L/23040.0L*pow(x_0, 2) - 7.0L/3600.0L*x_0*pow(x_1, 5) - 43.0L/1440.0L*x_0*pow(x_1, 4) - 257.0L/1440.0L*x_0*pow(x_1, 3) - 1483.0L/2880.0L*x_0*pow(x_1, 2) - 8147.0L/11520.0L*x_0*x_1 - 10357.0L/28800.0L*x_0 - 7.0L/4800.0L*pow(x_1, 6) - 367.0L/14400.0L*pow(x_1, 5) - 2063.0L/11520.0L*pow(x_1, 4) - 10987.0L/17280.0L*pow(x_1, 3) - 53363.0L/46080.0L*pow(x_1, 2) - 213187.0L/230400.0L*x_1 - 401033.0L/2764800.0L;
       return __pp_r402___result;
    }
    static inline O __pp_r403__(const I &x_0, const I &x_1) {
       O __pp_r403___result;
       __pp_r403___result = -1.0L/21600.0L*pow(x_0, 6) - 1.0L/1800.0L*pow(x_0, 5)*x_1 - 13.0L/7200.0L*pow(x_0, 5) - 1.0L/360.0L*pow(x_0, 4)*pow(x_1, 2) - 13.0L/720.0L*pow(x_0, 4)*x_1 - 169.0L/5760.0L*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) - 13.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) - 169.0L/720.0L*pow(x_0, 3)*x_1 - 2197.0L/8640.0L*pow(x_0, 3) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 4) - 13.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) - 169.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) - 2197.0L/1440.0L*pow(x_0, 2)*x_1 - 28561.0L/23040.0L*pow(x_0, 2) - 17.0L/3600.0L*x_0*pow(x_1, 5) - 103.0L/1440.0L*x_0*pow(x_1, 4) - 617.0L/1440.0L*x_0*pow(x_1, 3) - 3643.0L/2880.0L*x_0*pow(x_1, 2) - 21107.0L/11520.0L*x_0*x_1 - 29797.0L/28800.0L*x_0 - 83.0L/43200.0L*pow(x_1, 6) - 487.0L/14400.0L*pow(x_1, 5) - 2783.0L/11520.0L*pow(x_1, 4) - 15307.0L/17280.0L*pow(x_1, 3) - 79283.0L/46080.0L*pow(x_1, 2) - 368707.0L/230400.0L*x_1 - 1334153.0L/2764800.0L;
       return __pp_r403___result;
    }
    static inline O __pp_r404__(const I &x_0, const I &x_1) {
       O __pp_r404___result;
       __pp_r404___result = (1.0L/240.0L)*x_0*pow(x_1, 5) + (7.0L/96.0L)*x_0*pow(x_1, 4) + (49.0L/96.0L)*x_0*pow(x_1, 3) + (343.0L/192.0L)*x_0*pow(x_1, 2) + (2401.0L/768.0L)*x_0*x_1 + (16807.0L/7680.0L)*x_0 + (1.0L/960.0L)*pow(x_1, 6) + (23.0L/960.0L)*pow(x_1, 5) + (175.0L/768.0L)*pow(x_1, 4) + (147.0L/128.0L)*pow(x_1, 3) + (9947.0L/3072.0L)*pow(x_1, 2) + (74431.0L/15360.0L)*x_1 + 184877.0L/61440.0L;
       return __pp_r404___result;
    }
    static inline O __pp_r405__(const I &x_0, const I &x_1) {
       O __pp_r405___result;
       __pp_r405___result = (19.0L/43200.0L)*pow(x_0, 6) + (1.0L/400.0L)*pow(x_0, 5)*x_1 + (7.0L/960.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/32.0L)*pow(x_0, 4)*x_1 + (503.0L/11520.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (143.0L/1440.0L)*pow(x_0, 3)*x_1 + (419.0L/5760.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 37.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 301.0L/960.0L*pow(x_0, 2)*x_1 - 16561.0L/46080.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) + (1.0L/32.0L)*x_0*pow(x_1, 4) + (227.0L/1440.0L)*x_0*pow(x_1, 3) + (131.0L/320.0L)*x_0*pow(x_1, 2) + (3247.0L/5760.0L)*x_0*x_1 + (5213.0L/15360.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 1.0L/960.0L*pow(x_1, 5) - 127.0L/11520.0L*pow(x_1, 4) - 113.0L/5760.0L*pow(x_1, 3) + (5201.0L/46080.0L)*pow(x_1, 2) + (7261.0L/15360.0L)*x_1 + 42829.0L/86400.0L;
       return __pp_r405___result;
    }
    static inline O __pp_r406__(const I &x_0, const I &x_1) {
       O __pp_r406___result;
       __pp_r406___result = (1.0L/2400.0L)*pow(x_0, 6) + (1.0L/450.0L)*pow(x_0, 5)*x_1 + (47.0L/7200.0L)*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) + (17.0L/720.0L)*pow(x_0, 4)*x_1 + (191.0L/5760.0L)*pow(x_0, 4) + (1.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (11.0L/720.0L)*pow(x_0, 3)*x_1 - 37.0L/8640.0L*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) - 11.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 79.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) - 1117.0L/1440.0L*pow(x_0, 2)*x_1 - 15601.0L/23040.0L*pow(x_0, 2) - 7.0L/3600.0L*x_0*pow(x_1, 5) - 43.0L/1440.0L*x_0*pow(x_1, 4) - 257.0L/1440.0L*x_0*pow(x_1, 3) - 1483.0L/2880.0L*x_0*pow(x_1, 2) - 8147.0L/11520.0L*x_0*x_1 - 10357.0L/28800.0L*x_0 - 7.0L/4800.0L*pow(x_1, 6) - 367.0L/14400.0L*pow(x_1, 5) - 2063.0L/11520.0L*pow(x_1, 4) - 10987.0L/17280.0L*pow(x_1, 3) - 53363.0L/46080.0L*pow(x_1, 2) - 213187.0L/230400.0L*x_1 - 401033.0L/2764800.0L;
       return __pp_r406___result;
    }
    static inline O __pp_r407__(const I &x_0, const I &x_1) {
       O __pp_r407___result;
       __pp_r407___result = (29.0L/43200.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 + (1.0L/320.0L)*pow(x_0, 5) + (13.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/96.0L)*pow(x_0, 4)*x_1 + (863.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 37.0L/1440.0L*pow(x_0, 3)*x_1 - 301.0L/5760.0L*pow(x_0, 3) + (7.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (53.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (59.0L/960.0L)*pow(x_0, 2)*x_1 - 3601.0L/46080.0L*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (1.0L/96.0L)*x_0*pow(x_1, 4) + (47.0L/1440.0L)*x_0*pow(x_1, 3) + (11.0L/320.0L)*x_0*pow(x_1, 2) + (7.0L/5760.0L)*x_0*x_1 + (29.0L/15360.0L)*x_0 + (11.0L/43200.0L)*pow(x_1, 6) + (1.0L/320.0L)*pow(x_1, 5) + (233.0L/11520.0L)*pow(x_1, 4) + (607.0L/5760.0L)*pow(x_1, 3) + (18161.0L/46080.0L)*pow(x_1, 2) + (2489.0L/3072.0L)*x_1 + 57409.0L/86400.0L;
       return __pp_r407___result;
    }
    static inline O __pp_r408__(const I &x_0, const I &x_1) {
       O __pp_r408___result;
       __pp_r408___result = (7.0L/10800.0L)*pow(x_0, 6) + (1.0L/1200.0L)*pow(x_0, 5)*x_1 + (17.0L/7200.0L)*pow(x_0, 5) + (11.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (2.0L/45.0L)*pow(x_0, 4)*x_1 + (371.0L/5760.0L)*pow(x_0, 4) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 3) - 11.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 79.0L/720.0L*pow(x_0, 3)*x_1 - 1117.0L/8640.0L*pow(x_0, 3) - 1.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 17.0L/120.0L*pow(x_0, 2)*pow(x_1, 2) - 577.0L/1440.0L*pow(x_0, 2)*x_1 - 9121.0L/23040.0L*pow(x_0, 2) - 1.0L/300.0L*x_0*pow(x_1, 5) - 73.0L/1440.0L*x_0*pow(x_1, 4) - 437.0L/1440.0L*x_0*pow(x_1, 3) - 2563.0L/2880.0L*x_0*pow(x_1, 2) - 14627.0L/11520.0L*x_0*x_1 - 20077.0L/28800.0L*x_0 - 53.0L/43200.0L*pow(x_1, 6) - 307.0L/14400.0L*pow(x_1, 5) - 1703.0L/11520.0L*pow(x_1, 4) - 8827.0L/17280.0L*pow(x_1, 3) - 40403.0L/46080.0L*pow(x_1, 2) - 135427.0L/230400.0L*x_1 + 65527.0L/2764800.0L;
       return __pp_r408___result;
    }
    static inline O __pp_r409__(const I &x_0, const I &x_1) {
       O __pp_r409___result;
       __pp_r409___result = -1.0L/21600.0L*pow(x_0, 6) - 1.0L/1800.0L*pow(x_0, 5)*x_1 - 13.0L/7200.0L*pow(x_0, 5) - 1.0L/360.0L*pow(x_0, 4)*pow(x_1, 2) - 13.0L/720.0L*pow(x_0, 4)*x_1 - 169.0L/5760.0L*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) - 13.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) - 169.0L/720.0L*pow(x_0, 3)*x_1 - 2197.0L/8640.0L*pow(x_0, 3) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 4) - 13.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) - 169.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) - 2197.0L/1440.0L*pow(x_0, 2)*x_1 - 28561.0L/23040.0L*pow(x_0, 2) - 17.0L/3600.0L*x_0*pow(x_1, 5) - 103.0L/1440.0L*x_0*pow(x_1, 4) - 617.0L/1440.0L*x_0*pow(x_1, 3) - 3643.0L/2880.0L*x_0*pow(x_1, 2) - 21107.0L/11520.0L*x_0*x_1 - 29797.0L/28800.0L*x_0 - 83.0L/43200.0L*pow(x_1, 6) - 487.0L/14400.0L*pow(x_1, 5) - 2783.0L/11520.0L*pow(x_1, 4) - 15307.0L/17280.0L*pow(x_1, 3) - 79283.0L/46080.0L*pow(x_1, 2) - 368707.0L/230400.0L*x_1 - 1334153.0L/2764800.0L;
       return __pp_r409___result;
    }
    static inline O __pp_r410__(const I &x_0, const I &x_1) {
       O __pp_r410___result;
       __pp_r410___result = (1.0L/240.0L)*x_0*pow(x_1, 5) + (7.0L/96.0L)*x_0*pow(x_1, 4) + (49.0L/96.0L)*x_0*pow(x_1, 3) + (343.0L/192.0L)*x_0*pow(x_1, 2) + (2401.0L/768.0L)*x_0*x_1 + (16807.0L/7680.0L)*x_0 + (1.0L/960.0L)*pow(x_1, 6) + (23.0L/960.0L)*pow(x_1, 5) + (175.0L/768.0L)*pow(x_1, 4) + (147.0L/128.0L)*pow(x_1, 3) + (9947.0L/3072.0L)*pow(x_1, 2) + (74431.0L/15360.0L)*x_1 + 184877.0L/61440.0L;
       return __pp_r410___result;
    }
    static inline O __pp_r411__(const I &x_0, const I &x_1) {
       O __pp_r411___result;
       __pp_r411___result = (1.0L/5400.0L)*pow(x_0, 6) - 7.0L/3600.0L*pow(x_0, 5)*x_1 - 43.0L/7200.0L*pow(x_0, 5) + (1.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/360.0L)*pow(x_0, 4)*x_1 + (11.0L/5760.0L)*pow(x_0, 4) - 13.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 41.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 259.0L/720.0L*pow(x_0, 3)*x_1 - 3277.0L/8640.0L*pow(x_0, 3) - 11.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 37.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 31.0L/60.0L*pow(x_0, 2)*pow(x_1, 2) - 1657.0L/1440.0L*pow(x_0, 2)*x_1 - 22081.0L/23040.0L*pow(x_0, 2) - 11.0L/1800.0L*x_0*pow(x_1, 5) - 133.0L/1440.0L*x_0*pow(x_1, 4) - 797.0L/1440.0L*x_0*pow(x_1, 3) - 4723.0L/2880.0L*x_0*pow(x_1, 2) - 27587.0L/11520.0L*x_0*x_1 - 39517.0L/28800.0L*x_0 - 73.0L/43200.0L*pow(x_1, 6) - 427.0L/14400.0L*pow(x_1, 5) - 2423.0L/11520.0L*pow(x_1, 4) - 13147.0L/17280.0L*pow(x_1, 3) - 66323.0L/46080.0L*pow(x_1, 2) - 290947.0L/230400.0L*x_1 - 867593.0L/2764800.0L;
       return __pp_r411___result;
    }
    static inline O __pp_r412__(const I &x_0, const I &x_1) {
       O __pp_r412___result;
       __pp_r412___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 1.0L/240.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/48.0L)*pow(x_0, 4)*x_1 + (1.0L/32.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/24.0L*pow(x_0, 3)*pow(x_1, 2) - 1.0L/8.0L*pow(x_0, 3)*x_1 - 1.0L/8.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (3.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (3.0L/8.0L)*pow(x_0, 2)*x_1 + (9.0L/32.0L)*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) + (5.0L/96.0L)*x_0*pow(x_1, 4) + (37.0L/96.0L)*x_0*pow(x_1, 3) + (271.0L/192.0L)*x_0*pow(x_1, 2) + (1969.0L/768.0L)*x_0*x_1 + (2843.0L/1536.0L)*x_0 + (11.0L/8640.0L)*pow(x_1, 6) + (9.0L/320.0L)*pow(x_1, 5) + (199.0L/768.0L)*pow(x_1, 4) + (163.0L/128.0L)*pow(x_1, 3) + (10811.0L/3072.0L)*pow(x_1, 2) + (15923.0L/3072.0L)*x_1 + 39049.0L/12288.0L;
       return __pp_r412___result;
    }
    static inline O __pp_r413__(const I &x_0, const I &x_1) {
       O __pp_r413___result;
       __pp_r413___result = -7.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 + (41.0L/2880.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (7.0L/288.0L)*pow(x_0, 4)*x_1 + (463.0L/11520.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/36.0L)*pow(x_0, 3)*pow(x_1, 2) + (7.0L/160.0L)*pow(x_0, 3)*x_1 + (97.0L/17280.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/36.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) - 73.0L/2880.0L*pow(x_0, 2)*x_1 - 6101.0L/46080.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (1.0L/72.0L)*x_0*pow(x_1, 4) + (1.0L/20.0L)*x_0*pow(x_1, 3) + (7.0L/90.0L)*x_0*pow(x_1, 2) + (71.0L/1280.0L)*x_0*x_1 + (1337.0L/46080.0L)*x_0 + (1.0L/4320.0L)*pow(x_1, 6) + (1.0L/360.0L)*pow(x_1, 5) + (13.0L/720.0L)*pow(x_1, 4) + (53.0L/540.0L)*pow(x_1, 3) + (137.0L/360.0L)*pow(x_1, 2) + (3671.0L/4608.0L)*x_1 + 1821463.0L/2764800.0L;
       return __pp_r413___result;
    }
    static inline O __pp_r414__(const I &x_0, const I &x_1) {
       O __pp_r414___result;
       __pp_r414___result = -1.0L/1200.0L*pow(x_0, 6) + (19.0L/3600.0L)*pow(x_0, 5)*x_1 + (97.0L/7200.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 4)*x_1 + (19.0L/640.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 29.0L/720.0L*pow(x_0, 3)*x_1 - 617.0L/8640.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/30.0L*pow(x_0, 2)*pow(x_1, 3) - 31.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 39.0L/80.0L*pow(x_0, 2)*x_1 - 3457.0L/7680.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 17.0L/360.0L*x_0*pow(x_1, 4) - 103.0L/360.0L*x_0*pow(x_1, 3) - 1219.0L/1440.0L*x_0*pow(x_1, 2) - 7001.0L/5760.0L*x_0*x_1 - 77183.0L/115200.0L*x_0 - 1.0L/800.0L*pow(x_1, 6) - 13.0L/600.0L*pow(x_1, 5) - 3.0L/20.0L*pow(x_1, 4) - 373.0L/720.0L*pow(x_1, 3) - 3419.0L/3840.0L*pow(x_1, 2) - 5773.0L/9600.0L*x_1 + 8317.0L/460800.0L;
       return __pp_r414___result;
    }
    static inline O __pp_r415__(const I &x_0, const I &x_1) {
       O __pp_r415___result;
       __pp_r415___result = -7.0L/5400.0L*pow(x_0, 6) + (1.0L/400.0L)*pow(x_0, 5)*x_1 + (37.0L/7200.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/40.0L*pow(x_0, 4)*x_1 - 21.0L/640.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 31.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 209.0L/720.0L*pow(x_0, 3)*x_1 - 2777.0L/8640.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 99.0L/80.0L*pow(x_0, 2)*x_1 - 7777.0L/7680.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) - 4.0L/45.0L*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) - 2299.0L/1440.0L*x_0*pow(x_1, 2) - 13481.0L/5760.0L*x_0*x_1 - 154943.0L/115200.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 3.0L/100.0L*pow(x_1, 5) - 17.0L/80.0L*pow(x_1, 4) - 553.0L/720.0L*pow(x_1, 3) - 5579.0L/3840.0L*pow(x_1, 2) - 12253.0L/9600.0L*x_1 - 147203.0L/460800.0L;
       return __pp_r415___result;
    }
    static inline O __pp_r416__(const I &x_0, const I &x_1) {
       O __pp_r416___result;
       __pp_r416___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 91.0L/14400.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 17.0L/480.0L*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 31.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 - 6139.0L/17280.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 1193.0L/960.0L*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) - 4.0L/45.0L*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) - 2299.0L/1440.0L*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 - 310891.0L/230400.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 3.0L/100.0L*pow(x_1, 5) - 17.0L/80.0L*pow(x_1, 4) - 553.0L/720.0L*pow(x_1, 3) - 5579.0L/3840.0L*pow(x_1, 2) - 16339.0L/12800.0L*x_1 - 294811.0L/921600.0L;
       return __pp_r416___result;
    }
    static inline O __pp_r417__(const I &x_0, const I &x_1) {
       O __pp_r417___result;
       __pp_r417___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 91.0L/14400.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 17.0L/480.0L*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 31.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 - 6139.0L/17280.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 1193.0L/960.0L*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) - 4.0L/45.0L*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) - 2299.0L/1440.0L*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 - 310891.0L/230400.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 3.0L/100.0L*pow(x_1, 5) - 17.0L/80.0L*pow(x_1, 4) - 553.0L/720.0L*pow(x_1, 3) - 5579.0L/3840.0L*pow(x_1, 2) - 16339.0L/12800.0L*x_1 - 294811.0L/921600.0L;
       return __pp_r417___result;
    }
    static inline O __pp_r418__(const I &x_0, const I &x_1) {
       O __pp_r418___result;
       __pp_r418___result = -1.0L/800.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/144.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/144.0L*pow(x_0, 4)*x_1 - 1.0L/288.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 3)*x_1 - 29.0L/432.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/36.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (83.0L/288.0L)*pow(x_0, 2)*x_1 + (523.0L/2304.0L)*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) + (1.0L/18.0L)*x_0*pow(x_1, 4) + (29.0L/72.0L)*x_0*pow(x_1, 3) + (419.0L/288.0L)*x_0*pow(x_1, 2) + (377.0L/144.0L)*x_0*x_1 + (4327.0L/2304.0L)*x_0 + (1.0L/800.0L)*pow(x_1, 6) + (1.0L/36.0L)*pow(x_1, 5) + (37.0L/144.0L)*pow(x_1, 4) + (547.0L/432.0L)*pow(x_1, 3) + (8077.0L/2304.0L)*pow(x_1, 2) + (11911.0L/2304.0L)*x_1 + 10963.0L/3456.0L;
       return __pp_r418___result;
    }
    static inline O __pp_r419__(const I &x_0, const I &x_1) {
       O __pp_r419___result;
       __pp_r419___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 13.0L/2880.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 5.0L/288.0L*pow(x_0, 4)*x_1 - 83.0L/2304.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) - 19.0L/288.0L*pow(x_0, 3)*x_1 - 349.0L/3456.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/36.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (163.0L/576.0L)*pow(x_0, 2)*x_1 + (1933.0L/9216.0L)*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) + (1.0L/18.0L)*x_0*pow(x_1, 4) + (29.0L/72.0L)*x_0*pow(x_1, 3) + (419.0L/288.0L)*x_0*pow(x_1, 2) + (6029.0L/2304.0L)*x_0*x_1 + (86339.0L/46080.0L)*x_0 + (1.0L/800.0L)*pow(x_1, 6) + (1.0L/36.0L)*pow(x_1, 5) + (37.0L/144.0L)*pow(x_1, 4) + (547.0L/432.0L)*pow(x_1, 3) + (8077.0L/2304.0L)*pow(x_1, 2) + (119107.0L/23040.0L)*x_1 + 1753837.0L/552960.0L;
       return __pp_r419___result;
    }
    static inline O __pp_r420__(const I &x_0, const I &x_1) {
       O __pp_r420___result;
       __pp_r420___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/960.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 + (25.0L/128.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/12.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/64.0L)*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) + (1.0L/12.0L)*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) + (75.0L/32.0L)*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 + (3375.0L/1024.0L)*x_0 + (1.0L/675.0L)*pow(x_1, 6) + (1.0L/30.0L)*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) + (25.0L/16.0L)*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) + (3375.0L/512.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r420___result;
    }
    static inline O __pp_r421__(const I &x_0, const I &x_1) {
       O __pp_r421___result;
       __pp_r421___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 13.0L/2880.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 5.0L/288.0L*pow(x_0, 4)*x_1 - 83.0L/2304.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) - 19.0L/288.0L*pow(x_0, 3)*x_1 - 349.0L/3456.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/36.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (163.0L/576.0L)*pow(x_0, 2)*x_1 + (1933.0L/9216.0L)*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) + (1.0L/18.0L)*x_0*pow(x_1, 4) + (29.0L/72.0L)*x_0*pow(x_1, 3) + (419.0L/288.0L)*x_0*pow(x_1, 2) + (6029.0L/2304.0L)*x_0*x_1 + (86339.0L/46080.0L)*x_0 + (1.0L/800.0L)*pow(x_1, 6) + (1.0L/36.0L)*pow(x_1, 5) + (37.0L/144.0L)*pow(x_1, 4) + (547.0L/432.0L)*pow(x_1, 3) + (8077.0L/2304.0L)*pow(x_1, 2) + (119107.0L/23040.0L)*x_1 + 1753837.0L/552960.0L;
       return __pp_r421___result;
    }
    static inline O __pp_r422__(const I &x_0, const I &x_1) {
       O __pp_r422___result;
       __pp_r422___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/960.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 + (25.0L/128.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/12.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/64.0L)*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) + (1.0L/12.0L)*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) + (75.0L/32.0L)*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 + (3375.0L/1024.0L)*x_0 + (1.0L/675.0L)*pow(x_1, 6) + (1.0L/30.0L)*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) + (25.0L/16.0L)*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) + (3375.0L/512.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r422___result;
    }
    static inline O __pp_r423__(const I &x_0, const I &x_1) {
       O __pp_r423___result;
       __pp_r423___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 91.0L/14400.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 17.0L/480.0L*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 31.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 - 6139.0L/17280.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 1193.0L/960.0L*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) - 4.0L/45.0L*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) - 2299.0L/1440.0L*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 - 310891.0L/230400.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 3.0L/100.0L*pow(x_1, 5) - 17.0L/80.0L*pow(x_1, 4) - 553.0L/720.0L*pow(x_1, 3) - 5579.0L/3840.0L*pow(x_1, 2) - 16339.0L/12800.0L*x_1 - 294811.0L/921600.0L;
       return __pp_r423___result;
    }
    static inline O __pp_r424__(const I &x_0, const I &x_1) {
       O __pp_r424___result;
       __pp_r424___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 13.0L/2880.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 5.0L/288.0L*pow(x_0, 4)*x_1 - 83.0L/2304.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) - 19.0L/288.0L*pow(x_0, 3)*x_1 - 349.0L/3456.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/36.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (163.0L/576.0L)*pow(x_0, 2)*x_1 + (1933.0L/9216.0L)*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) + (1.0L/18.0L)*x_0*pow(x_1, 4) + (29.0L/72.0L)*x_0*pow(x_1, 3) + (419.0L/288.0L)*x_0*pow(x_1, 2) + (6029.0L/2304.0L)*x_0*x_1 + (86339.0L/46080.0L)*x_0 + (1.0L/800.0L)*pow(x_1, 6) + (1.0L/36.0L)*pow(x_1, 5) + (37.0L/144.0L)*pow(x_1, 4) + (547.0L/432.0L)*pow(x_1, 3) + (8077.0L/2304.0L)*pow(x_1, 2) + (119107.0L/23040.0L)*x_1 + 1753837.0L/552960.0L;
       return __pp_r424___result;
    }
    static inline O __pp_r425__(const I &x_0, const I &x_1) {
       O __pp_r425___result;
       __pp_r425___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/960.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 + (25.0L/128.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/12.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/64.0L)*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) + (1.0L/12.0L)*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) + (75.0L/32.0L)*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 + (3375.0L/1024.0L)*x_0 + (1.0L/675.0L)*pow(x_1, 6) + (1.0L/30.0L)*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) + (25.0L/16.0L)*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) + (3375.0L/512.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r425___result;
    }
    static inline O __pp_r426__(const I &x_0, const I &x_1) {
       O __pp_r426___result;
       __pp_r426___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 91.0L/14400.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 17.0L/480.0L*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 31.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 - 6139.0L/17280.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 1193.0L/960.0L*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) - 4.0L/45.0L*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) - 2299.0L/1440.0L*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 - 310891.0L/230400.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 3.0L/100.0L*pow(x_1, 5) - 17.0L/80.0L*pow(x_1, 4) - 553.0L/720.0L*pow(x_1, 3) - 5579.0L/3840.0L*pow(x_1, 2) - 16339.0L/12800.0L*x_1 - 294811.0L/921600.0L;
       return __pp_r426___result;
    }
    static inline O __pp_r427__(const I &x_0, const I &x_1) {
       O __pp_r427___result;
       __pp_r427___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 13.0L/2880.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 5.0L/288.0L*pow(x_0, 4)*x_1 - 83.0L/2304.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/72.0L*pow(x_0, 3)*pow(x_1, 2) - 19.0L/288.0L*pow(x_0, 3)*x_1 - 349.0L/3456.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/36.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (163.0L/576.0L)*pow(x_0, 2)*x_1 + (1933.0L/9216.0L)*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) + (1.0L/18.0L)*x_0*pow(x_1, 4) + (29.0L/72.0L)*x_0*pow(x_1, 3) + (419.0L/288.0L)*x_0*pow(x_1, 2) + (6029.0L/2304.0L)*x_0*x_1 + (86339.0L/46080.0L)*x_0 + (1.0L/800.0L)*pow(x_1, 6) + (1.0L/36.0L)*pow(x_1, 5) + (37.0L/144.0L)*pow(x_1, 4) + (547.0L/432.0L)*pow(x_1, 3) + (8077.0L/2304.0L)*pow(x_1, 2) + (119107.0L/23040.0L)*x_1 + 1753837.0L/552960.0L;
       return __pp_r427___result;
    }
    static inline O __pp_r428__(const I &x_0, const I &x_1) {
       O __pp_r428___result;
       __pp_r428___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/960.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 + (25.0L/128.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/12.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/64.0L)*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) + (1.0L/12.0L)*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) + (75.0L/32.0L)*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 + (3375.0L/1024.0L)*x_0 + (1.0L/675.0L)*pow(x_1, 6) + (1.0L/30.0L)*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) + (25.0L/16.0L)*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) + (3375.0L/512.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r428___result;
    }
    static inline O __pp_r429__(const I &x_0, const I &x_1) {
       O __pp_r429___result;
       __pp_r429___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/960.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 + (25.0L/128.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/12.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/64.0L)*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) + (1.0L/12.0L)*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) + (75.0L/32.0L)*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 + (3375.0L/1024.0L)*x_0 + (1.0L/675.0L)*pow(x_1, 6) + (1.0L/30.0L)*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) + (25.0L/16.0L)*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) + (3375.0L/512.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r429___result;
    }
    static inline O __pp_r430__(const I &x_0, const I &x_1) {
       O __pp_r430___result;
       __pp_r430___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/960.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 + (25.0L/128.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/12.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/64.0L)*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) + (1.0L/12.0L)*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) + (75.0L/32.0L)*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 + (3375.0L/1024.0L)*x_0 + (1.0L/675.0L)*pow(x_1, 6) + (1.0L/30.0L)*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) + (25.0L/16.0L)*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) + (3375.0L/512.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r430___result;
    }
    static inline O __pp_r431__(const I &x_0, const I &x_1) {
       O __pp_r431___result;
       __pp_r431___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/960.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 + (25.0L/128.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/12.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/64.0L)*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) + (1.0L/12.0L)*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) + (75.0L/32.0L)*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 + (3375.0L/1024.0L)*x_0 + (1.0L/675.0L)*pow(x_1, 6) + (1.0L/30.0L)*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) + (25.0L/16.0L)*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) + (3375.0L/512.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r431___result;
    }
    static inline O __pp_r432__(const I &x_0, const I &x_1) {
       O __pp_r432___result;
       __pp_r432___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/960.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 + (25.0L/128.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/12.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/64.0L)*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/225.0L)*x_0*pow(x_1, 5) + (1.0L/12.0L)*x_0*pow(x_1, 4) + (5.0L/8.0L)*x_0*pow(x_1, 3) + (75.0L/32.0L)*x_0*pow(x_1, 2) + (1125.0L/256.0L)*x_0*x_1 + (3375.0L/1024.0L)*x_0 + (1.0L/675.0L)*pow(x_1, 6) + (1.0L/30.0L)*pow(x_1, 5) + (5.0L/16.0L)*pow(x_1, 4) + (25.0L/16.0L)*pow(x_1, 3) + (1125.0L/256.0L)*pow(x_1, 2) + (3375.0L/512.0L)*x_1 + 16875.0L/4096.0L;
       return __pp_r432___result;
    }
    static inline O __pp_r433__(const I &x_0, const I &x_1) {
       O __pp_r433___result;
       __pp_r433___result = (31.0L/43200.0L)*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 + (7.0L/1600.0L)*pow(x_0, 5) + (17.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (31.0L/480.0L)*pow(x_0, 4)*x_1 + (205.0L/2304.0L)*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/20.0L)*pow(x_0, 3)*pow(x_1, 2) + (25.0L/288.0L)*pow(x_0, 3)*x_1 + (37.0L/1152.0L)*pow(x_0, 3) + (23.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (17.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (43.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (109.0L/192.0L)*pow(x_0, 2)*x_1 + (9521.0L/46080.0L)*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) + (7.0L/120.0L)*x_0*pow(x_1, 4) + (2.0L/9.0L)*x_0*pow(x_1, 3) + (19.0L/48.0L)*x_0*pow(x_1, 2) + (3761.0L/11520.0L)*x_0*x_1 + (8261.0L/76800.0L)*x_0 + (47.0L/21600.0L)*pow(x_1, 6) + (17.0L/800.0L)*pow(x_1, 5) + (5.0L/72.0L)*pow(x_1, 4) + (37.0L/576.0L)*pow(x_1, 3) - 619.0L/11520.0L*pow(x_1, 2) + (83.0L/19200.0L)*x_1 + 101869.0L/552960.0L;
       return __pp_r433___result;
    }
    static inline O __pp_r434__(const I &x_0, const I &x_1) {
       O __pp_r434___result;
       __pp_r434___result = -11.0L/14400.0L*pow(x_0, 6) + (11.0L/1800.0L)*pow(x_0, 5)*x_1 + (223.0L/14400.0L)*pow(x_0, 5) + (1.0L/160.0L)*pow(x_0, 4)*pow(x_1, 2) + (53.0L/1440.0L)*pow(x_0, 4)*x_1 + (125.0L/2304.0L)*pow(x_0, 4) + (13.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 + (311.0L/3456.0L)*pow(x_0, 3) + (7.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (23.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (19.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (277.0L/576.0L)*pow(x_0, 2)*x_1 + (7021.0L/46080.0L)*pow(x_0, 2) + (11.0L/1800.0L)*x_0*pow(x_1, 5) + (89.0L/1440.0L)*x_0*pow(x_1, 4) + (23.0L/96.0L)*x_0*pow(x_1, 3) + (253.0L/576.0L)*x_0*pow(x_1, 2) + (731.0L/1920.0L)*x_0*x_1 + (31033.0L/230400.0L)*x_0 + (31.0L/14400.0L)*pow(x_1, 6) + (301.0L/14400.0L)*pow(x_1, 5) + (155.0L/2304.0L)*pow(x_1, 4) + (197.0L/3456.0L)*pow(x_1, 3) - 3101.0L/46080.0L*pow(x_1, 2) - 2129.0L/230400.0L*x_1 + 12343.0L/69120.0L;
       return __pp_r434___result;
    }
    static inline O __pp_r435__(const I &x_0, const I &x_1) {
       O __pp_r435___result;
       __pp_r435___result = -43.0L/43200.0L*pow(x_0, 6) + (3.0L/400.0L)*pow(x_0, 5)*x_1 + (263.0L/14400.0L)*pow(x_0, 5) + (1.0L/360.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/480.0L)*pow(x_0, 4)*x_1 + (31.0L/768.0L)*pow(x_0, 4) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 3) + (19.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (61.0L/288.0L)*pow(x_0, 3)*x_1 + (439.0L/3456.0L)*pow(x_0, 3) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/10.0L)*pow(x_0, 2)*pow(x_1, 3) + (5.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (71.0L/192.0L)*pow(x_0, 2)*x_1 + (1487.0L/15360.0L)*pow(x_0, 2) + (3.0L/400.0L)*x_0*pow(x_1, 5) + (109.0L/1440.0L)*x_0*pow(x_1, 4) + (85.0L/288.0L)*x_0*pow(x_1, 3) + (317.0L/576.0L)*x_0*pow(x_1, 2) + (2833.0L/5760.0L)*x_0*x_1 + (41273.0L/230400.0L)*x_0 + (83.0L/43200.0L)*pow(x_1, 6) + (29.0L/1600.0L)*pow(x_1, 5) + (41.0L/768.0L)*pow(x_1, 4) + (23.0L/1152.0L)*pow(x_1, 3) - 629.0L/5120.0L*pow(x_1, 2) - 4123.0L/76800.0L*x_1 + 3773.0L/23040.0L;
       return __pp_r435___result;
    }
    static inline O __pp_r436__(const I &x_0, const I &x_1) {
       O __pp_r436___result;
       __pp_r436___result = (7.0L/14400.0L)*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (167.0L/14400.0L)*pow(x_0, 5) + (1.0L/120.0L)*pow(x_0, 4)*pow(x_1, 2) + (19.0L/480.0L)*pow(x_0, 4)*x_1 + (203.0L/3840.0L)*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) + (4.0L/45.0L)*pow(x_0, 3)*pow(x_1, 2) + (269.0L/1440.0L)*pow(x_0, 3)*x_1 + (1979.0L/17280.0L)*pow(x_0, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 4) + (13.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (53.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (373.0L/960.0L)*pow(x_0, 2)*x_1 + (319.0L/3072.0L)*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) + (53.0L/720.0L)*x_0*pow(x_1, 4) + (13.0L/45.0L)*x_0*pow(x_1, 3) + (779.0L/1440.0L)*x_0*pow(x_1, 2) + (1117.0L/2304.0L)*x_0*x_1 + (40787.0L/230400.0L)*x_0 + (7.0L/3600.0L)*pow(x_1, 6) + (11.0L/600.0L)*pow(x_1, 5) + (13.0L/240.0L)*pow(x_1, 4) + (31.0L/1440.0L)*pow(x_1, 3) - 31.0L/256.0L*pow(x_1, 2) - 2021.0L/38400.0L*x_1 + 151163.0L/921600.0L;
       return __pp_r436___result;
    }
    static inline O __pp_r437__(const I &x_0, const I &x_1) {
       O __pp_r437___result;
       __pp_r437___result = (7.0L/14400.0L)*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (107.0L/14400.0L)*pow(x_0, 5) + (1.0L/120.0L)*pow(x_0, 4)*pow(x_1, 2) + (19.0L/480.0L)*pow(x_0, 4)*x_1 + (163.0L/3840.0L)*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) + (4.0L/45.0L)*pow(x_0, 3)*pow(x_1, 2) + (269.0L/1440.0L)*pow(x_0, 3)*x_1 + (1799.0L/17280.0L)*pow(x_0, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 4) + (13.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (53.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (373.0L/960.0L)*pow(x_0, 2)*x_1 + (101.0L/1024.0L)*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) + (53.0L/720.0L)*x_0*pow(x_1, 4) + (13.0L/45.0L)*x_0*pow(x_1, 3) + (779.0L/1440.0L)*x_0*pow(x_1, 2) + (1117.0L/2304.0L)*x_0*x_1 + (40487.0L/230400.0L)*x_0 + (7.0L/3600.0L)*pow(x_1, 6) + (11.0L/600.0L)*pow(x_1, 5) + (13.0L/240.0L)*pow(x_1, 4) + (31.0L/1440.0L)*pow(x_1, 3) - 31.0L/256.0L*pow(x_1, 2) - 2021.0L/38400.0L*x_1 + 151043.0L/921600.0L;
       return __pp_r437___result;
    }
    static inline O __pp_r438__(const I &x_0, const I &x_1) {
       O __pp_r438___result;
       __pp_r438___result = (29.0L/43200.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 + (1.0L/320.0L)*pow(x_0, 5) + (13.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/96.0L)*pow(x_0, 4)*x_1 + (863.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 37.0L/1440.0L*pow(x_0, 3)*x_1 - 301.0L/5760.0L*pow(x_0, 3) + (7.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (53.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (59.0L/960.0L)*pow(x_0, 2)*x_1 - 3601.0L/46080.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/24.0L*x_0*pow(x_1, 4) - 41.0L/180.0L*x_0*pow(x_1, 3) - 37.0L/60.0L*x_0*pow(x_1, 2) - 9361.0L/11520.0L*x_0*x_1 - 6221.0L/15360.0L*x_0 - 17.0L/21600.0L*pow(x_1, 6) - 3.0L/160.0L*pow(x_1, 5) - 7.0L/45.0L*pow(x_1, 4) - 1759.0L/2880.0L*pow(x_1, 3) - 13741.0L/11520.0L*pow(x_1, 2) - 49.0L/48.0L*x_1 - 553537.0L/2764800.0L;
       return __pp_r438___result;
    }
    static inline O __pp_r439__(const I &x_0, const I &x_1) {
       O __pp_r439___result;
       __pp_r439___result = -7.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 + (41.0L/2880.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (7.0L/288.0L)*pow(x_0, 4)*x_1 + (463.0L/11520.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/36.0L)*pow(x_0, 3)*pow(x_1, 2) + (7.0L/160.0L)*pow(x_0, 3)*x_1 + (97.0L/17280.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/36.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) - 73.0L/2880.0L*pow(x_0, 2)*x_1 - 6101.0L/46080.0L*pow(x_0, 2) - 1.0L/360.0L*x_0*pow(x_1, 5) - 11.0L/288.0L*x_0*pow(x_1, 4) - 101.0L/480.0L*x_0*pow(x_1, 3) - 1651.0L/2880.0L*x_0*pow(x_1, 2) - 91.0L/120.0L*x_0*x_1 - 17413.0L/46080.0L*x_0 - 7.0L/8640.0L*pow(x_1, 6) - 11.0L/576.0L*pow(x_1, 5) - 1817.0L/11520.0L*pow(x_1, 4) - 10679.0L/17280.0L*pow(x_1, 3) - 55589.0L/46080.0L*pow(x_1, 2) - 9533.0L/9216.0L*x_1 - 284581.0L/1382400.0L;
       return __pp_r439___result;
    }
    static inline O __pp_r440__(const I &x_0, const I &x_1) {
       O __pp_r440___result;
       __pp_r440___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 + (49.0L/2880.0L)*pow(x_0, 5) + (1.0L/96.0L)*pow(x_0, 4)*x_1 + (101.0L/3840.0L)*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/18.0L)*pow(x_0, 3)*pow(x_1, 2) + (143.0L/1440.0L)*pow(x_0, 3)*x_1 + (737.0L/17280.0L)*pow(x_0, 3) - 1.0L/40.0L*pow(x_0, 2)*pow(x_1, 2) - 131.0L/960.0L*pow(x_0, 2)*x_1 - 2887.0L/15360.0L*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) - 7.0L/288.0L*x_0*pow(x_1, 4) - 223.0L/1440.0L*x_0*pow(x_1, 3) - 1331.0L/2880.0L*x_0*pow(x_1, 2) - 233.0L/360.0L*x_0*x_1 - 3073.0L/9216.0L*x_0 - 1.0L/960.0L*pow(x_1, 6) - 7.0L/320.0L*pow(x_1, 5) - 659.0L/3840.0L*pow(x_1, 4) - 3773.0L/5760.0L*pow(x_1, 3) - 6461.0L/5120.0L*pow(x_1, 2) - 16571.0L/15360.0L*x_1 - 101687.0L/460800.0L;
       return __pp_r440___result;
    }
    static inline O __pp_r441__(const I &x_0, const I &x_1) {
       O __pp_r441___result;
       __pp_r441___result = (19.0L/43200.0L)*pow(x_0, 6) + (1.0L/400.0L)*pow(x_0, 5)*x_1 + (149.0L/14400.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (13.0L/480.0L)*pow(x_0, 4)*x_1 + (149.0L/3840.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (107.0L/1440.0L)*pow(x_0, 3)*x_1 + (521.0L/17280.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) - 1.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 113.0L/960.0L*pow(x_0, 2)*x_1 - 2779.0L/15360.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) - 19.0L/720.0L*x_0*pow(x_1, 4) - 29.0L/180.0L*x_0*pow(x_1, 3) - 679.0L/1440.0L*x_0*pow(x_1, 2) - 7537.0L/11520.0L*x_0*x_1 - 77311.0L/230400.0L*x_0 - 11.0L/10800.0L*pow(x_1, 6) - 13.0L/600.0L*pow(x_1, 5) - 41.0L/240.0L*pow(x_1, 4) - 941.0L/1440.0L*pow(x_1, 3) - 1613.0L/1280.0L*pow(x_1, 2) - 41387.0L/38400.0L*x_1 - 203131.0L/921600.0L;
       return __pp_r441___result;
    }
    static inline O __pp_r442__(const I &x_0, const I &x_1) {
       O __pp_r442___result;
       __pp_r442___result = (19.0L/43200.0L)*pow(x_0, 6) + (1.0L/400.0L)*pow(x_0, 5)*x_1 + (89.0L/14400.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (13.0L/480.0L)*pow(x_0, 4)*x_1 + (109.0L/3840.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (107.0L/1440.0L)*pow(x_0, 3)*x_1 + (341.0L/17280.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) - 1.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 113.0L/960.0L*pow(x_0, 2)*x_1 - 953.0L/5120.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) - 19.0L/720.0L*x_0*pow(x_1, 4) - 29.0L/180.0L*x_0*pow(x_1, 3) - 679.0L/1440.0L*x_0*pow(x_1, 2) - 7537.0L/11520.0L*x_0*x_1 - 77611.0L/230400.0L*x_0 - 11.0L/10800.0L*pow(x_1, 6) - 13.0L/600.0L*pow(x_1, 5) - 41.0L/240.0L*pow(x_1, 4) - 941.0L/1440.0L*pow(x_1, 3) - 1613.0L/1280.0L*pow(x_1, 2) - 41387.0L/38400.0L*x_1 - 203251.0L/921600.0L;
       return __pp_r442___result;
    }
    static inline O __pp_r443__(const I &x_0, const I &x_1) {
       O __pp_r443___result;
       __pp_r443___result = (7.0L/14400.0L)*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 + (107.0L/14400.0L)*pow(x_0, 5) + (1.0L/120.0L)*pow(x_0, 4)*pow(x_1, 2) + (19.0L/480.0L)*pow(x_0, 4)*x_1 + (163.0L/3840.0L)*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) + (4.0L/45.0L)*pow(x_0, 3)*pow(x_1, 2) + (269.0L/1440.0L)*pow(x_0, 3)*x_1 + (1799.0L/17280.0L)*pow(x_0, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 4) + (13.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (53.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (373.0L/960.0L)*pow(x_0, 2)*x_1 + (101.0L/1024.0L)*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) + (53.0L/720.0L)*x_0*pow(x_1, 4) + (13.0L/45.0L)*x_0*pow(x_1, 3) + (779.0L/1440.0L)*x_0*pow(x_1, 2) + (1117.0L/2304.0L)*x_0*x_1 + (40487.0L/230400.0L)*x_0 + (7.0L/3600.0L)*pow(x_1, 6) + (11.0L/600.0L)*pow(x_1, 5) + (13.0L/240.0L)*pow(x_1, 4) + (31.0L/1440.0L)*pow(x_1, 3) - 31.0L/256.0L*pow(x_1, 2) - 2021.0L/38400.0L*x_1 + 151043.0L/921600.0L;
       return __pp_r443___result;
    }
    static inline O __pp_r444__(const I &x_0, const I &x_1) {
       O __pp_r444___result;
       __pp_r444___result = (31.0L/43200.0L)*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 + (29.0L/4800.0L)*pow(x_0, 5) + (17.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (67.0L/1440.0L)*pow(x_0, 4)*x_1 + (529.0L/11520.0L)*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (3.0L/40.0L)*pow(x_0, 3)*pow(x_1, 2) + (83.0L/480.0L)*pow(x_0, 3)*x_1 + (191.0L/1920.0L)*pow(x_0, 3) + (23.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (11.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (169.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (1159.0L/2880.0L)*pow(x_0, 2)*x_1 + (941.0L/9216.0L)*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) + (1.0L/15.0L)*x_0*pow(x_1, 4) + (11.0L/40.0L)*x_0*pow(x_1, 3) + (253.0L/480.0L)*x_0*pow(x_1, 2) + (367.0L/768.0L)*x_0*x_1 + (4463.0L/25600.0L)*x_0 + (47.0L/21600.0L)*pow(x_1, 6) + (71.0L/3600.0L)*pow(x_1, 5) + (83.0L/1440.0L)*pow(x_1, 4) + (113.0L/4320.0L)*pow(x_1, 3) - 271.0L/2304.0L*pow(x_1, 2) - 5903.0L/115200.0L*x_1 + 453769.0L/2764800.0L;
       return __pp_r444___result;
    }
    static inline O __pp_r445__(const I &x_0, const I &x_1) {
       O __pp_r445___result;
       __pp_r445___result = -11.0L/14400.0L*pow(x_0, 6) + (11.0L/1800.0L)*pow(x_0, 5)*x_1 + (11.0L/2880.0L)*pow(x_0, 5) + (1.0L/160.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/96.0L)*pow(x_0, 4)*x_1 + (57.0L/1280.0L)*pow(x_0, 4) + (13.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (253.0L/1440.0L)*pow(x_0, 3)*x_1 + (1711.0L/17280.0L)*pow(x_0, 3) + (7.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/8.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) + (129.0L/320.0L)*pow(x_0, 2)*x_1 + (1567.0L/15360.0L)*pow(x_0, 2) + (11.0L/1800.0L)*x_0*pow(x_1, 5) + (19.0L/288.0L)*x_0*pow(x_1, 4) + (397.0L/1440.0L)*x_0*pow(x_1, 3) + (1517.0L/2880.0L)*x_0*pow(x_1, 2) + (2753.0L/5760.0L)*x_0*x_1 + (8033.0L/46080.0L)*x_0 + (31.0L/14400.0L)*pow(x_1, 6) + (19.0L/960.0L)*pow(x_1, 5) + (221.0L/3840.0L)*pow(x_1, 4) + (151.0L/5760.0L)*pow(x_1, 3) - 1807.0L/15360.0L*pow(x_1, 2) - 787.0L/15360.0L*x_1 + 18907.0L/115200.0L;
       return __pp_r445___result;
    }
    static inline O __pp_r446__(const I &x_0, const I &x_1) {
       O __pp_r446___result;
       __pp_r446___result = (19.0L/43200.0L)*pow(x_0, 6) + (1.0L/400.0L)*pow(x_0, 5)*x_1 + (89.0L/14400.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (13.0L/480.0L)*pow(x_0, 4)*x_1 + (109.0L/3840.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (107.0L/1440.0L)*pow(x_0, 3)*x_1 + (341.0L/17280.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) - 1.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 113.0L/960.0L*pow(x_0, 2)*x_1 - 953.0L/5120.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) - 19.0L/720.0L*x_0*pow(x_1, 4) - 29.0L/180.0L*x_0*pow(x_1, 3) - 679.0L/1440.0L*x_0*pow(x_1, 2) - 7537.0L/11520.0L*x_0*x_1 - 77611.0L/230400.0L*x_0 - 11.0L/10800.0L*pow(x_1, 6) - 13.0L/600.0L*pow(x_1, 5) - 41.0L/240.0L*pow(x_1, 4) - 941.0L/1440.0L*pow(x_1, 3) - 1613.0L/1280.0L*pow(x_1, 2) - 41387.0L/38400.0L*x_1 - 203251.0L/921600.0L;
       return __pp_r446___result;
    }
    static inline O __pp_r447__(const I &x_0, const I &x_1) {
       O __pp_r447___result;
       __pp_r447___result = (29.0L/43200.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 + (23.0L/4800.0L)*pow(x_0, 5) + (13.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (49.0L/1440.0L)*pow(x_0, 4)*x_1 + (367.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/40.0L)*pow(x_0, 3)*pow(x_1, 2) + (29.0L/480.0L)*pow(x_0, 3)*x_1 + (29.0L/1920.0L)*pow(x_0, 3) + (7.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/45.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 299.0L/2880.0L*pow(x_0, 2)*x_1 - 8417.0L/46080.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/30.0L*x_0*pow(x_1, 4) - 7.0L/40.0L*x_0*pow(x_1, 3) - 233.0L/480.0L*x_0*pow(x_1, 2) - 2539.0L/3840.0L*x_0*x_1 - 8659.0L/25600.0L*x_0 - 17.0L/21600.0L*pow(x_1, 6) - 73.0L/3600.0L*pow(x_1, 5) - 241.0L/1440.0L*pow(x_1, 4) - 2803.0L/4320.0L*pow(x_1, 3) - 14477.0L/11520.0L*pow(x_1, 2) - 124001.0L/115200.0L*x_1 - 609113.0L/2764800.0L;
       return __pp_r447___result;
    }
    static inline O __pp_r448__(const I &x_0, const I &x_1) {
       O __pp_r448___result;
       __pp_r448___result = -7.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 + (37.0L/14400.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (19.0L/480.0L)*pow(x_0, 4)*x_1 + (39.0L/1280.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (91.0L/1440.0L)*pow(x_0, 3)*x_1 + (253.0L/17280.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/40.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) - 33.0L/320.0L*pow(x_0, 2)*x_1 - 2807.0L/15360.0L*pow(x_0, 2) - 1.0L/360.0L*x_0*pow(x_1, 5) - 49.0L/1440.0L*x_0*pow(x_1, 4) - 251.0L/1440.0L*x_0*pow(x_1, 3) - 1399.0L/2880.0L*x_0*pow(x_1, 2) - 119.0L/180.0L*x_0*x_1 - 77933.0L/230400.0L*x_0 - 7.0L/8640.0L*pow(x_1, 6) - 97.0L/4800.0L*pow(x_1, 5) - 643.0L/3840.0L*pow(x_1, 4) - 3737.0L/5760.0L*pow(x_1, 3) - 19303.0L/15360.0L*pow(x_1, 2) - 82667.0L/76800.0L*x_1 - 101519.0L/460800.0L;
       return __pp_r448___result;
    }
    static inline O __pp_r449__(const I &x_0, const I &x_1) {
       O __pp_r449___result;
       __pp_r449___result = (29.0L/43200.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 + (73.0L/14400.0L)*pow(x_0, 5) + (13.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (53.0L/1440.0L)*pow(x_0, 4)*x_1 + (431.0L/11520.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (13.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (151.0L/1440.0L)*pow(x_0, 3)*x_1 + (1033.0L/17280.0L)*pow(x_0, 3) + (7.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (2.0L/45.0L)*pow(x_0, 2)*pow(x_1, 3) + (71.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (473.0L/2880.0L)*pow(x_0, 2)*x_1 - 97.0L/46080.0L*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) + (29.0L/1440.0L)*x_0*pow(x_1, 4) + (139.0L/1440.0L)*x_0*pow(x_1, 3) + (551.0L/2880.0L)*x_0*pow(x_1, 2) + (959.0L/5760.0L)*x_0*x_1 + (13843.0L/230400.0L)*x_0 + (11.0L/43200.0L)*pow(x_1, 6) + (181.0L/14400.0L)*pow(x_1, 5) + (761.0L/11520.0L)*pow(x_1, 4) + (1849.0L/17280.0L)*pow(x_1, 3) + (1697.0L/46080.0L)*pow(x_1, 2) + (16951.0L/230400.0L)*x_1 + 34841.0L/172800.0L;
       return __pp_r449___result;
    }
    static inline O __pp_r450__(const I &x_0, const I &x_1) {
       O __pp_r450___result;
       __pp_r450___result = (31.0L/43200.0L)*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 + (29.0L/4800.0L)*pow(x_0, 5) + (17.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (67.0L/1440.0L)*pow(x_0, 4)*x_1 + (529.0L/11520.0L)*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (3.0L/40.0L)*pow(x_0, 3)*pow(x_1, 2) + (83.0L/480.0L)*pow(x_0, 3)*x_1 + (191.0L/1920.0L)*pow(x_0, 3) + (23.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (11.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (169.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (1159.0L/2880.0L)*pow(x_0, 2)*x_1 + (941.0L/9216.0L)*pow(x_0, 2) + (1.0L/100.0L)*x_0*pow(x_1, 5) + (47.0L/480.0L)*x_0*pow(x_1, 4) + (59.0L/160.0L)*x_0*pow(x_1, 3) + (641.0L/960.0L)*x_0*pow(x_1, 2) + (7.0L/12.0L)*x_0*x_1 + (5273.0L/25600.0L)*x_0 + (139.0L/43200.0L)*pow(x_1, 6) + (629.0L/14400.0L)*pow(x_1, 5) + (2329.0L/11520.0L)*pow(x_1, 4) + (7337.0L/17280.0L)*pow(x_1, 3) + (4181.0L/9216.0L)*pow(x_1, 2) + (84179.0L/230400.0L)*x_1 + 396377.0L/1382400.0L;
       return __pp_r450___result;
    }
    static inline O __pp_r451__(const I &x_0, const I &x_1) {
       O __pp_r451___result;
       __pp_r451___result = -7.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 + (41.0L/14400.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (61.0L/1440.0L)*pow(x_0, 4)*x_1 + (83.0L/2304.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (11.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (31.0L/288.0L)*pow(x_0, 3)*x_1 + (205.0L/3456.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (17.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (95.0L/576.0L)*pow(x_0, 2)*x_1 - 101.0L/46080.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (7.0L/360.0L)*x_0*pow(x_1, 4) + (7.0L/72.0L)*x_0*pow(x_1, 3) + (55.0L/288.0L)*x_0*pow(x_1, 2) + (1919.0L/11520.0L)*x_0*x_1 + (13841.0L/230400.0L)*x_0 + (1.0L/4320.0L)*pow(x_1, 6) + (91.0L/7200.0L)*pow(x_1, 5) + (19.0L/288.0L)*pow(x_1, 4) + (185.0L/1728.0L)*pow(x_1, 3) + (53.0L/1440.0L)*pow(x_1, 2) + (2119.0L/28800.0L)*x_1 + 111491.0L/552960.0L;
       return __pp_r451___result;
    }
    static inline O __pp_r452__(const I &x_0, const I &x_1) {
       O __pp_r452___result;
       __pp_r452___result = -11.0L/14400.0L*pow(x_0, 6) + (11.0L/1800.0L)*pow(x_0, 5)*x_1 + (11.0L/2880.0L)*pow(x_0, 5) + (1.0L/160.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/96.0L)*pow(x_0, 4)*x_1 + (57.0L/1280.0L)*pow(x_0, 4) + (13.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (253.0L/1440.0L)*pow(x_0, 3)*x_1 + (1711.0L/17280.0L)*pow(x_0, 3) + (7.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/8.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) + (129.0L/320.0L)*pow(x_0, 2)*x_1 + (1567.0L/15360.0L)*pow(x_0, 2) + (37.0L/3600.0L)*x_0*pow(x_1, 5) + (7.0L/72.0L)*x_0*pow(x_1, 4) + (133.0L/360.0L)*x_0*pow(x_1, 3) + (961.0L/1440.0L)*x_0*pow(x_1, 2) + (6721.0L/11520.0L)*x_0*x_1 + (9491.0L/46080.0L)*x_0 + (23.0L/7200.0L)*pow(x_1, 6) + (7.0L/160.0L)*pow(x_1, 5) + (97.0L/480.0L)*pow(x_1, 4) + (1223.0L/2880.0L)*pow(x_1, 3) + (871.0L/1920.0L)*pow(x_1, 2) + (1403.0L/3840.0L)*x_1 + 264251.0L/921600.0L;
       return __pp_r452___result;
    }
    static inline O __pp_r453__(const I &x_0, const I &x_1) {
       O __pp_r453___result;
       __pp_r453___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 91.0L/14400.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 17.0L/480.0L*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 31.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 - 6139.0L/17280.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 1193.0L/960.0L*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) - 4.0L/45.0L*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) - 2299.0L/1440.0L*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 - 310891.0L/230400.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 41.0L/1200.0L*pow(x_1, 5) - 127.0L/480.0L*pow(x_1, 4) - 1481.0L/1440.0L*pow(x_1, 3) - 2693.0L/1280.0L*pow(x_1, 2) - 80267.0L/38400.0L*x_1 - 669811.0L/921600.0L;
       return __pp_r453___result;
    }
    static inline O __pp_r454__(const I &x_0, const I &x_1) {
       O __pp_r454___result;
       __pp_r454___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 91.0L/14400.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 17.0L/480.0L*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 31.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 - 6139.0L/17280.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 1193.0L/960.0L*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) - 4.0L/45.0L*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) - 2299.0L/1440.0L*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 - 310891.0L/230400.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 41.0L/1200.0L*pow(x_1, 5) - 127.0L/480.0L*pow(x_1, 4) - 1481.0L/1440.0L*pow(x_1, 3) - 2693.0L/1280.0L*pow(x_1, 2) - 80267.0L/38400.0L*x_1 - 669811.0L/921600.0L;
       return __pp_r454___result;
    }
    static inline O __pp_r455__(const I &x_0, const I &x_1) {
       O __pp_r455___result;
       __pp_r455___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 91.0L/14400.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 17.0L/480.0L*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 31.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 - 6139.0L/17280.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 1193.0L/960.0L*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) - 4.0L/45.0L*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) - 2299.0L/1440.0L*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 - 310891.0L/230400.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 41.0L/1200.0L*pow(x_1, 5) - 127.0L/480.0L*pow(x_1, 4) - 1481.0L/1440.0L*pow(x_1, 3) - 2693.0L/1280.0L*pow(x_1, 2) - 80267.0L/38400.0L*x_1 - 669811.0L/921600.0L;
       return __pp_r455___result;
    }
    static inline O __pp_r456__(const I &x_0, const I &x_1) {
       O __pp_r456___result;
       __pp_r456___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 91.0L/14400.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 17.0L/480.0L*pow(x_0, 4)*x_1 - 251.0L/3840.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 31.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) - 433.0L/1440.0L*pow(x_0, 3)*x_1 - 6139.0L/17280.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/60.0L*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) - 1193.0L/960.0L*pow(x_0, 2)*x_1 - 5273.0L/5120.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) - 4.0L/45.0L*x_0*pow(x_1, 4) - 193.0L/360.0L*x_0*pow(x_1, 3) - 2299.0L/1440.0L*x_0*pow(x_1, 2) - 26977.0L/11520.0L*x_0*x_1 - 310891.0L/230400.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 41.0L/1200.0L*pow(x_1, 5) - 127.0L/480.0L*pow(x_1, 4) - 1481.0L/1440.0L*pow(x_1, 3) - 2693.0L/1280.0L*pow(x_1, 2) - 80267.0L/38400.0L*x_1 - 669811.0L/921600.0L;
       return __pp_r456___result;
    }
    static inline O __pp_r457__(const I &x_0, const I &x_1) {
       O __pp_r457___result;
       __pp_r457___result = -1.0L/43200.0L*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 - 37.0L/4800.0L*pow(x_0, 5) - 1.0L/720.0L*pow(x_0, 4)*pow(x_1, 2) - 41.0L/1440.0L*pow(x_0, 4)*x_1 - 713.0L/11520.0L*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/10.0L*pow(x_0, 3)*pow(x_1, 2) - 151.0L/480.0L*pow(x_0, 3)*x_1 - 691.0L/1920.0L*pow(x_0, 3) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 4) - 37.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 263.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 3539.0L/2880.0L*pow(x_0, 2)*x_1 - 47297.0L/46080.0L*pow(x_0, 2) - 13.0L/1800.0L*x_0*pow(x_1, 5) - 23.0L/240.0L*x_0*pow(x_1, 4) - 11.0L/20.0L*x_0*pow(x_1, 3) - 773.0L/480.0L*x_0*pow(x_1, 2) - 9019.0L/3840.0L*x_0*x_1 - 34579.0L/25600.0L*x_0 - 1.0L/675.0L*pow(x_1, 6) - 59.0L/1800.0L*pow(x_1, 5) - 47.0L/180.0L*pow(x_1, 4) - 4423.0L/4320.0L*pow(x_1, 3) - 24197.0L/11520.0L*pow(x_1, 2) - 240641.0L/115200.0L*x_1 - 2008793.0L/2764800.0L;
       return __pp_r457___result;
    }
    static inline O __pp_r458__(const I &x_0, const I &x_1) {
       O __pp_r458___result;
       __pp_r458___result = -1.0L/43200.0L*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 - 37.0L/4800.0L*pow(x_0, 5) - 1.0L/720.0L*pow(x_0, 4)*pow(x_1, 2) - 41.0L/1440.0L*pow(x_0, 4)*x_1 - 713.0L/11520.0L*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/10.0L*pow(x_0, 3)*pow(x_1, 2) - 151.0L/480.0L*pow(x_0, 3)*x_1 - 691.0L/1920.0L*pow(x_0, 3) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 4) - 37.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 263.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 3539.0L/2880.0L*pow(x_0, 2)*x_1 - 47297.0L/46080.0L*pow(x_0, 2) - 13.0L/1800.0L*x_0*pow(x_1, 5) - 23.0L/240.0L*x_0*pow(x_1, 4) - 11.0L/20.0L*x_0*pow(x_1, 3) - 773.0L/480.0L*x_0*pow(x_1, 2) - 9019.0L/3840.0L*x_0*x_1 - 34579.0L/25600.0L*x_0 - 1.0L/675.0L*pow(x_1, 6) - 59.0L/1800.0L*pow(x_1, 5) - 47.0L/180.0L*pow(x_1, 4) - 4423.0L/4320.0L*pow(x_1, 3) - 24197.0L/11520.0L*pow(x_1, 2) - 240641.0L/115200.0L*x_1 - 2008793.0L/2764800.0L;
       return __pp_r458___result;
    }
    static inline O __pp_r459__(const I &x_0, const I &x_1) {
       O __pp_r459___result;
       __pp_r459___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 - 143.0L/14400.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 11.0L/480.0L*pow(x_0, 4)*x_1 - 81.0L/1280.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) - 19.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) - 449.0L/1440.0L*pow(x_0, 3)*x_1 - 6227.0L/17280.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/10.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 393.0L/320.0L*pow(x_0, 2)*x_1 - 15767.0L/15360.0L*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) - 139.0L/1440.0L*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) - 4639.0L/2880.0L*x_0*pow(x_1, 2) - 1691.0L/720.0L*x_0*x_1 - 311213.0L/230400.0L*x_0 - 13.0L/8640.0L*pow(x_1, 6) - 157.0L/4800.0L*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) - 5897.0L/5760.0L*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) - 160427.0L/76800.0L*x_1 - 334799.0L/460800.0L;
       return __pp_r459___result;
    }
    static inline O __pp_r460__(const I &x_0, const I &x_1) {
       O __pp_r460___result;
       __pp_r460___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 - 143.0L/14400.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 11.0L/480.0L*pow(x_0, 4)*x_1 - 81.0L/1280.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) - 19.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) - 449.0L/1440.0L*pow(x_0, 3)*x_1 - 6227.0L/17280.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/10.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 393.0L/320.0L*pow(x_0, 2)*x_1 - 15767.0L/15360.0L*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) - 139.0L/1440.0L*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) - 4639.0L/2880.0L*x_0*pow(x_1, 2) - 1691.0L/720.0L*x_0*x_1 - 311213.0L/230400.0L*x_0 - 13.0L/8640.0L*pow(x_1, 6) - 157.0L/4800.0L*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) - 5897.0L/5760.0L*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) - 160427.0L/76800.0L*x_1 - 334799.0L/460800.0L;
       return __pp_r460___result;
    }
    static inline O __pp_r461__(const I &x_0, const I &x_1) {
       O __pp_r461___result;
       __pp_r461___result = -1.0L/2160.0L*pow(x_0, 6) - 1.0L/360.0L*pow(x_0, 5)*x_1 - 79.0L/7200.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 13.0L/240.0L*pow(x_0, 4)*x_1 - 17.0L/160.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) - 19.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) - 73.0L/180.0L*pow(x_0, 3)*x_1 - 4531.0L/8640.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/10.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 219.0L/160.0L*pow(x_0, 2)*x_1 - 1247.0L/960.0L*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) - 139.0L/1440.0L*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) - 4639.0L/2880.0L*x_0*pow(x_1, 2) - 28271.0L/11520.0L*x_0*x_1 - 90257.0L/57600.0L*x_0 - 13.0L/8640.0L*pow(x_1, 6) - 157.0L/4800.0L*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) - 5897.0L/5760.0L*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) - 162857.0L/76800.0L*x_1 - 731563.0L/921600.0L;
       return __pp_r461___result;
    }
    static inline O __pp_r462__(const I &x_0, const I &x_1) {
       O __pp_r462___result;
       __pp_r462___result = -1.0L/2160.0L*pow(x_0, 6) - 1.0L/360.0L*pow(x_0, 5)*x_1 - 79.0L/7200.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 13.0L/240.0L*pow(x_0, 4)*x_1 - 17.0L/160.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) - 19.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) - 73.0L/180.0L*pow(x_0, 3)*x_1 - 4531.0L/8640.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/10.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 219.0L/160.0L*pow(x_0, 2)*x_1 - 1247.0L/960.0L*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) - 139.0L/1440.0L*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) - 4639.0L/2880.0L*x_0*pow(x_1, 2) - 28271.0L/11520.0L*x_0*x_1 - 90257.0L/57600.0L*x_0 - 13.0L/8640.0L*pow(x_1, 6) - 157.0L/4800.0L*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) - 5897.0L/5760.0L*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) - 162857.0L/76800.0L*x_1 - 731563.0L/921600.0L;
       return __pp_r462___result;
    }
    static inline O __pp_r463__(const I &x_0, const I &x_1) {
       O __pp_r463___result;
       __pp_r463___result = (1.0L/7200.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*x_1 + (7.0L/1440.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (7.0L/180.0L)*pow(x_0, 3)*x_1 + (589.0L/8640.0L)*pow(x_0, 3) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) + (589.0L/1440.0L)*pow(x_0, 2)*x_1 + (1379.0L/2880.0L)*pow(x_0, 2) - 1.0L/240.0L*x_0*pow(x_1, 5) - 59.0L/1440.0L*x_0*pow(x_1, 4) - 151.0L/1440.0L*x_0*pow(x_1, 3) + (481.0L/2880.0L)*x_0*pow(x_1, 2) + (12689.0L/11520.0L)*x_0*x_1 + (73583.0L/57600.0L)*x_0 - 1.0L/960.0L*pow(x_1, 6) - 311.0L/14400.0L*pow(x_1, 5) - 1729.0L/11520.0L*pow(x_1, 4) - 7451.0L/17280.0L*pow(x_1, 3) - 14869.0L/46080.0L*pow(x_1, 2) + (166789.0L/230400.0L)*x_1 + 3048191.0L/2764800.0L;
       return __pp_r463___result;
    }
    static inline O __pp_r464__(const I &x_0, const I &x_1) {
       O __pp_r464___result;
       __pp_r464___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/960.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 + (25.0L/128.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/12.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/64.0L)*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/3600.0L)*x_0*pow(x_1, 5) + (1.0L/32.0L)*x_0*pow(x_1, 4) + (35.0L/96.0L)*x_0*pow(x_1, 3) + (325.0L/192.0L)*x_0*pow(x_1, 2) + (1375.0L/384.0L)*x_0*x_1 + (8875.0L/3072.0L)*x_0 + (19.0L/43200.0L)*pow(x_1, 6) + (7.0L/960.0L)*pow(x_1, 5) + (65.0L/768.0L)*pow(x_1, 4) + (75.0L/128.0L)*pow(x_1, 3) + (6625.0L/3072.0L)*pow(x_1, 2) + (12125.0L/3072.0L)*x_1 + 4375.0L/1536.0L;
       return __pp_r464___result;
    }
    static inline O __pp_r465__(const I &x_0, const I &x_1) {
       O __pp_r465___result;
       __pp_r465___result = -1.0L/2160.0L*pow(x_0, 6) - 1.0L/360.0L*pow(x_0, 5)*x_1 - 79.0L/7200.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 13.0L/240.0L*pow(x_0, 4)*x_1 - 17.0L/160.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) - 19.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) - 73.0L/180.0L*pow(x_0, 3)*x_1 - 4531.0L/8640.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/10.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 219.0L/160.0L*pow(x_0, 2)*x_1 - 1247.0L/960.0L*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) - 139.0L/1440.0L*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) - 4639.0L/2880.0L*x_0*pow(x_1, 2) - 28271.0L/11520.0L*x_0*x_1 - 90257.0L/57600.0L*x_0 - 13.0L/8640.0L*pow(x_1, 6) - 157.0L/4800.0L*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) - 5897.0L/5760.0L*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) - 162857.0L/76800.0L*x_1 - 731563.0L/921600.0L;
       return __pp_r465___result;
    }
    static inline O __pp_r466__(const I &x_0, const I &x_1) {
       O __pp_r466___result;
       __pp_r466___result = -1.0L/2160.0L*pow(x_0, 6) - 1.0L/360.0L*pow(x_0, 5)*x_1 - 79.0L/7200.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) - 13.0L/240.0L*pow(x_0, 4)*x_1 - 17.0L/160.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) - 19.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) - 73.0L/180.0L*pow(x_0, 3)*x_1 - 4531.0L/8640.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/10.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 219.0L/160.0L*pow(x_0, 2)*x_1 - 1247.0L/960.0L*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) - 139.0L/1440.0L*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) - 4639.0L/2880.0L*x_0*pow(x_1, 2) - 28271.0L/11520.0L*x_0*x_1 - 90257.0L/57600.0L*x_0 - 13.0L/8640.0L*pow(x_1, 6) - 157.0L/4800.0L*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) - 5897.0L/5760.0L*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) - 162857.0L/76800.0L*x_1 - 731563.0L/921600.0L;
       return __pp_r466___result;
    }
    static inline O __pp_r467__(const I &x_0, const I &x_1) {
       O __pp_r467___result;
       __pp_r467___result = (1.0L/7200.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*x_1 + (7.0L/1440.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (7.0L/180.0L)*pow(x_0, 3)*x_1 + (589.0L/8640.0L)*pow(x_0, 3) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) + (589.0L/1440.0L)*pow(x_0, 2)*x_1 + (1379.0L/2880.0L)*pow(x_0, 2) - 1.0L/240.0L*x_0*pow(x_1, 5) - 59.0L/1440.0L*x_0*pow(x_1, 4) - 151.0L/1440.0L*x_0*pow(x_1, 3) + (481.0L/2880.0L)*x_0*pow(x_1, 2) + (12689.0L/11520.0L)*x_0*x_1 + (73583.0L/57600.0L)*x_0 - 1.0L/960.0L*pow(x_1, 6) - 311.0L/14400.0L*pow(x_1, 5) - 1729.0L/11520.0L*pow(x_1, 4) - 7451.0L/17280.0L*pow(x_1, 3) - 14869.0L/46080.0L*pow(x_1, 2) + (166789.0L/230400.0L)*x_1 + 3048191.0L/2764800.0L;
       return __pp_r467___result;
    }
    static inline O __pp_r468__(const I &x_0, const I &x_1) {
       O __pp_r468___result;
       __pp_r468___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (1.0L/960.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 + (25.0L/128.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/12.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/64.0L)*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/3600.0L)*x_0*pow(x_1, 5) + (1.0L/32.0L)*x_0*pow(x_1, 4) + (35.0L/96.0L)*x_0*pow(x_1, 3) + (325.0L/192.0L)*x_0*pow(x_1, 2) + (1375.0L/384.0L)*x_0*x_1 + (8875.0L/3072.0L)*x_0 + (19.0L/43200.0L)*pow(x_1, 6) + (7.0L/960.0L)*pow(x_1, 5) + (65.0L/768.0L)*pow(x_1, 4) + (75.0L/128.0L)*pow(x_1, 3) + (6625.0L/3072.0L)*pow(x_1, 2) + (12125.0L/3072.0L)*x_1 + 4375.0L/1536.0L;
       return __pp_r468___result;
    }
    static inline O __pp_r469__(const I &x_0, const I &x_1) {
       O __pp_r469___result;
       __pp_r469___result = (11.0L/4320.0L)*pow(x_0, 6) - 1.0L/120.0L*pow(x_0, 5)*x_1 + (1.0L/144.0L)*pow(x_0, 4)*pow(x_1, 2) + (49.0L/5760.0L)*pow(x_0, 4) + (1.0L/80.0L)*pow(x_0, 3)*x_1 + (1.0L/144.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) - 367.0L/4608.0L*pow(x_0, 2) + (1.0L/120.0L)*x_0*pow(x_1, 5) - 1.0L/80.0L*x_0*pow(x_1, 3) + (11.0L/4320.0L)*pow(x_1, 6) + (49.0L/5760.0L)*pow(x_1, 4) - 367.0L/4608.0L*pow(x_1, 2) + 124937.0L/691200.0L;
       return __pp_r469___result;
    }
    static inline O __pp_r470__(const I &x_0, const I &x_1) {
       O __pp_r470___result;
       __pp_r470___result = (1.0L/400.0L)*pow(x_0, 6) - 2.0L/225.0L*pow(x_0, 5)*x_1 - 1.0L/7200.0L*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/720.0L*pow(x_0, 4)*x_1 + (1.0L/120.0L)*pow(x_0, 4) - 1.0L/135.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/90.0L)*pow(x_0, 3)*x_1 - 1.0L/8640.0L*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/1440.0L*pow(x_0, 2)*x_1 - 51.0L/640.0L*pow(x_0, 2) - 1.0L/1800.0L*x_0*pow(x_1, 5) - 1.0L/90.0L*x_0*pow(x_1, 4) - 13.0L/720.0L*x_0*pow(x_1, 3) - 1.0L/720.0L*x_0*pow(x_1, 2) - 1.0L/5760.0L*x_0*x_1 - 1.0L/115200.0L*x_0 - 1.0L/2400.0L*pow(x_1, 6) - 1.0L/225.0L*pow(x_1, 5) + (11.0L/1920.0L)*pow(x_1, 4) - 1.0L/1080.0L*pow(x_1, 3) - 613.0L/7680.0L*pow(x_1, 2) - 1.0L/57600.0L*x_1 + 83291.0L/460800.0L;
       return __pp_r470___result;
    }
    static inline O __pp_r471__(const I &x_0, const I &x_1) {
       O __pp_r471___result;
       __pp_r471___result = -1.0L/2400.0L*pow(x_0, 6) + (1.0L/1800.0L)*pow(x_0, 5)*x_1 - 1.0L/225.0L*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/90.0L)*pow(x_0, 4)*x_1 + (11.0L/1920.0L)*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (13.0L/720.0L)*pow(x_0, 3)*x_1 - 1.0L/1080.0L*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/720.0L)*pow(x_0, 2)*x_1 - 613.0L/7680.0L*pow(x_0, 2) + (2.0L/225.0L)*x_0*pow(x_1, 5) - 1.0L/720.0L*x_0*pow(x_1, 4) - 1.0L/90.0L*x_0*pow(x_1, 3) - 1.0L/1440.0L*x_0*pow(x_1, 2) + (1.0L/5760.0L)*x_0*x_1 - 1.0L/57600.0L*x_0 + (1.0L/400.0L)*pow(x_1, 6) + (1.0L/7200.0L)*pow(x_1, 5) + (1.0L/120.0L)*pow(x_1, 4) + (1.0L/8640.0L)*pow(x_1, 3) - 51.0L/640.0L*pow(x_1, 2) + (1.0L/115200.0L)*x_1 + 83291.0L/460800.0L;
       return __pp_r471___result;
    }
    static inline O __pp_r472__(const I &x_0, const I &x_1) {
       O __pp_r472___result;
       __pp_r472___result = -1.0L/2160.0L*pow(x_0, 6) - 11.0L/2400.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (7.0L/720.0L)*pow(x_0, 4)*x_1 + (1.0L/180.0L)*pow(x_0, 4) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/60.0L)*pow(x_0, 3)*x_1 - 1.0L/960.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/1440.0L)*pow(x_0, 2)*x_1 - 23.0L/288.0L*pow(x_0, 2) - 1.0L/80.0L*x_0*pow(x_1, 4) - 1.0L/60.0L*x_0*pow(x_1, 3) - 1.0L/480.0L*x_0*pow(x_1, 2) - 1.0L/38400.0L*x_0 - 1.0L/2160.0L*pow(x_1, 6) - 31.0L/7200.0L*pow(x_1, 5) + (1.0L/180.0L)*pow(x_1, 4) - 7.0L/8640.0L*pow(x_1, 3) - 23.0L/288.0L*pow(x_1, 2) - 1.0L/115200.0L*x_1 + 15617.0L/86400.0L;
       return __pp_r472___result;
    }
    static inline O __pp_r473__(const I &x_0, const I &x_1) {
       O __pp_r473___result;
       __pp_r473___result = -7.0L/4800.0L*pow(x_0, 6) + (17.0L/3600.0L)*pow(x_0, 5)*x_1 - 19.0L/14400.0L*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) + (31.0L/1440.0L)*pow(x_0, 4)*x_1 + (67.0L/3840.0L)*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 + (209.0L/17280.0L)*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) + (19.0L/2880.0L)*pow(x_0, 2)*x_1 - 1121.0L/15360.0L*pow(x_0, 2) + (2.0L/225.0L)*x_0*pow(x_1, 5) - 1.0L/720.0L*x_0*pow(x_1, 4) - 1.0L/90.0L*x_0*pow(x_1, 3) - 1.0L/1440.0L*x_0*pow(x_1, 2) + (17.0L/11520.0L)*x_0*x_1 + (401.0L/230400.0L)*x_0 + (1.0L/400.0L)*pow(x_1, 6) + (1.0L/7200.0L)*pow(x_1, 5) + (1.0L/120.0L)*pow(x_1, 4) + (1.0L/8640.0L)*pow(x_1, 3) - 51.0L/640.0L*pow(x_1, 2) + (1.0L/7200.0L)*x_1 + 166747.0L/921600.0L;
       return __pp_r473___result;
    }
    static inline O __pp_r474__(const I &x_0, const I &x_1) {
       O __pp_r474___result;
       __pp_r474___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 - 7.0L/4800.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (29.0L/1440.0L)*pow(x_0, 4)*x_1 + (199.0L/11520.0L)*pow(x_0, 4) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) + (13.0L/480.0L)*pow(x_0, 3)*x_1 + (23.0L/1920.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) + (17.0L/2880.0L)*pow(x_0, 2)*x_1 - 673.0L/9216.0L*pow(x_0, 2) - 1.0L/80.0L*x_0*pow(x_1, 4) - 1.0L/60.0L*x_0*pow(x_1, 3) - 1.0L/480.0L*x_0*pow(x_1, 2) + (1.0L/768.0L)*x_0*x_1 + (133.0L/76800.0L)*x_0 - 1.0L/2160.0L*pow(x_1, 6) - 31.0L/7200.0L*pow(x_1, 5) + (1.0L/180.0L)*pow(x_1, 4) - 7.0L/8640.0L*pow(x_1, 3) - 23.0L/288.0L*pow(x_1, 2) + (7.0L/57600.0L)*x_1 + 500239.0L/2764800.0L;
       return __pp_r474___result;
    }
    static inline O __pp_r475__(const I &x_0, const I &x_1) {
       O __pp_r475___result;
       __pp_r475___result = -11.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 - 1.0L/14400.0L*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) + (13.0L/480.0L)*pow(x_0, 4)*x_1 + (239.0L/11520.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (59.0L/1440.0L)*pow(x_0, 3)*x_1 + (287.0L/17280.0L)*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/240.0L)*pow(x_0, 2)*pow(x_1, 2) + (19.0L/960.0L)*pow(x_0, 2)*x_1 - 641.0L/9216.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 1.0L/180.0L*x_0*pow(x_1, 4) - 1.0L/360.0L*x_0*pow(x_1, 3) + (17.0L/1440.0L)*x_0*pow(x_1, 2) + (19.0L/2304.0L)*x_0*x_1 + (719.0L/230400.0L)*x_0 - 1.0L/4320.0L*pow(x_1, 6) - 7.0L/2400.0L*pow(x_1, 5) + (13.0L/1440.0L)*pow(x_1, 4) + (11.0L/2880.0L)*pow(x_1, 3) - 11.0L/144.0L*pow(x_1, 2) + (29.0L/19200.0L)*x_1 + 500879.0L/2764800.0L;
       return __pp_r475___result;
    }
    static inline O __pp_r476__(const I &x_0, const I &x_1) {
       O __pp_r476___result;
       __pp_r476___result = -7.0L/4800.0L*pow(x_0, 6) + (17.0L/3600.0L)*pow(x_0, 5)*x_1 - 19.0L/14400.0L*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) + (31.0L/1440.0L)*pow(x_0, 4)*x_1 + (67.0L/3840.0L)*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 + (209.0L/17280.0L)*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) + (19.0L/2880.0L)*pow(x_0, 2)*x_1 - 1121.0L/15360.0L*pow(x_0, 2) + (2.0L/225.0L)*x_0*pow(x_1, 5) - 1.0L/720.0L*x_0*pow(x_1, 4) - 1.0L/90.0L*x_0*pow(x_1, 3) - 1.0L/1440.0L*x_0*pow(x_1, 2) + (17.0L/11520.0L)*x_0*x_1 + (401.0L/230400.0L)*x_0 + (1.0L/400.0L)*pow(x_1, 6) + (1.0L/7200.0L)*pow(x_1, 5) + (1.0L/120.0L)*pow(x_1, 4) + (1.0L/8640.0L)*pow(x_1, 3) - 51.0L/640.0L*pow(x_1, 2) + (1.0L/7200.0L)*x_1 + 166747.0L/921600.0L;
       return __pp_r476___result;
    }
    static inline O __pp_r477__(const I &x_0, const I &x_1) {
       O __pp_r477___result;
       __pp_r477___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 - 7.0L/4800.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (29.0L/1440.0L)*pow(x_0, 4)*x_1 + (199.0L/11520.0L)*pow(x_0, 4) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) + (13.0L/480.0L)*pow(x_0, 3)*x_1 + (23.0L/1920.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) + (17.0L/2880.0L)*pow(x_0, 2)*x_1 - 673.0L/9216.0L*pow(x_0, 2) - 1.0L/80.0L*x_0*pow(x_1, 4) - 1.0L/60.0L*x_0*pow(x_1, 3) - 1.0L/480.0L*x_0*pow(x_1, 2) + (1.0L/768.0L)*x_0*x_1 + (133.0L/76800.0L)*x_0 - 1.0L/2160.0L*pow(x_1, 6) - 31.0L/7200.0L*pow(x_1, 5) + (1.0L/180.0L)*pow(x_1, 4) - 7.0L/8640.0L*pow(x_1, 3) - 23.0L/288.0L*pow(x_1, 2) + (7.0L/57600.0L)*x_1 + 500239.0L/2764800.0L;
       return __pp_r477___result;
    }
    static inline O __pp_r478__(const I &x_0, const I &x_1) {
       O __pp_r478___result;
       __pp_r478___result = -11.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 - 1.0L/14400.0L*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) + (13.0L/480.0L)*pow(x_0, 4)*x_1 + (239.0L/11520.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (59.0L/1440.0L)*pow(x_0, 3)*x_1 + (287.0L/17280.0L)*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/240.0L)*pow(x_0, 2)*pow(x_1, 2) + (19.0L/960.0L)*pow(x_0, 2)*x_1 - 641.0L/9216.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 1.0L/180.0L*x_0*pow(x_1, 4) - 1.0L/360.0L*x_0*pow(x_1, 3) + (17.0L/1440.0L)*x_0*pow(x_1, 2) + (19.0L/2304.0L)*x_0*x_1 + (719.0L/230400.0L)*x_0 - 1.0L/4320.0L*pow(x_1, 6) - 7.0L/2400.0L*pow(x_1, 5) + (13.0L/1440.0L)*pow(x_1, 4) + (11.0L/2880.0L)*pow(x_1, 3) - 11.0L/144.0L*pow(x_1, 2) + (29.0L/19200.0L)*x_1 + 500879.0L/2764800.0L;
       return __pp_r478___result;
    }
    static inline O __pp_r479__(const I &x_0, const I &x_1) {
       O __pp_r479___result;
       __pp_r479___result = -7.0L/5400.0L*pow(x_0, 6) + (19.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/3600.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/40.0L)*pow(x_0, 4)*x_1 + (23.0L/1152.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 3)*x_1 + (13.0L/864.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 2)*x_1 - 1643.0L/23040.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/45.0L*x_0*pow(x_1, 4) - 1.0L/36.0L*x_0*pow(x_1, 3) - 1.0L/144.0L*x_0*pow(x_1, 2) + (7.0L/5760.0L)*x_0*x_1 + (119.0L/57600.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 1.0L/800.0L*pow(x_1, 5) + (5.0L/288.0L)*pow(x_1, 4) + (7.0L/576.0L)*pow(x_1, 3) - 841.0L/11520.0L*pow(x_1, 2) + (77.0L/38400.0L)*x_1 + 50087.0L/276480.0L;
       return __pp_r479___result;
    }
    static inline O __pp_r480__(const I &x_0, const I &x_1) {
       O __pp_r480___result;
       __pp_r480___result = -11.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 + (1.0L/14400.0L)*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 4)*x_1 + (17.0L/768.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/96.0L)*pow(x_0, 3)*x_1 + (77.0L/3456.0L)*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (31.0L/576.0L)*pow(x_0, 2)*x_1 - 887.0L/15360.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (1.0L/180.0L)*x_0*pow(x_1, 4) + (1.0L/24.0L)*x_0*pow(x_1, 3) + (23.0L/288.0L)*x_0*pow(x_1, 2) + (71.0L/1280.0L)*x_0*x_1 + (3601.0L/230400.0L)*x_0 - 1.0L/4320.0L*pow(x_1, 6) + (71.0L/7200.0L)*pow(x_1, 5) + (5.0L/96.0L)*pow(x_1, 4) + (121.0L/1728.0L)*pow(x_1, 3) - 3.0L/160.0L*pow(x_1, 2) + (839.0L/28800.0L)*x_1 + 34433.0L/184320.0L;
       return __pp_r480___result;
    }
    static inline O __pp_r481__(const I &x_0, const I &x_1) {
       O __pp_r481___result;
       __pp_r481___result = -7.0L/5400.0L*pow(x_0, 6) + (19.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/3600.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/40.0L)*pow(x_0, 4)*x_1 + (23.0L/1152.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 3)*x_1 + (13.0L/864.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 2)*x_1 - 1643.0L/23040.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/45.0L*x_0*pow(x_1, 4) - 1.0L/36.0L*x_0*pow(x_1, 3) - 1.0L/144.0L*x_0*pow(x_1, 2) + (7.0L/5760.0L)*x_0*x_1 + (119.0L/57600.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 23.0L/2400.0L*pow(x_1, 5) - 1.0L/288.0L*pow(x_1, 4) - 5.0L/576.0L*pow(x_1, 3) - 961.0L/11520.0L*pow(x_1, 2) - 23.0L/38400.0L*x_1 + 10003.0L/55296.0L;
       return __pp_r481___result;
    }
    static inline O __pp_r482__(const I &x_0, const I &x_1) {
       O __pp_r482___result;
       __pp_r482___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 + (41.0L/14400.0L)*pow(x_0, 5) + (61.0L/1440.0L)*pow(x_0, 4)*x_1 + (83.0L/2304.0L)*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) + (11.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (31.0L/288.0L)*pow(x_0, 3)*x_1 + (205.0L/3456.0L)*pow(x_0, 3) + (17.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (95.0L/576.0L)*pow(x_0, 2)*x_1 - 101.0L/46080.0L*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) + (7.0L/360.0L)*x_0*pow(x_1, 4) + (7.0L/72.0L)*x_0*pow(x_1, 3) + (55.0L/288.0L)*x_0*pow(x_1, 2) + (1919.0L/11520.0L)*x_0*x_1 + (13841.0L/230400.0L)*x_0 + (91.0L/7200.0L)*pow(x_1, 5) + (19.0L/288.0L)*pow(x_1, 4) + (185.0L/1728.0L)*pow(x_1, 3) + (53.0L/1440.0L)*pow(x_1, 2) + (2119.0L/28800.0L)*x_1 + 111491.0L/552960.0L;
       return __pp_r482___result;
    }
    static inline O __pp_r483__(const I &x_0, const I &x_1) {
       O __pp_r483___result;
       __pp_r483___result = -43.0L/43200.0L*pow(x_0, 6) + (3.0L/400.0L)*pow(x_0, 5)*x_1 + (11.0L/2880.0L)*pow(x_0, 5) + (1.0L/360.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/96.0L)*pow(x_0, 4)*x_1 + (57.0L/1280.0L)*pow(x_0, 4) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (253.0L/1440.0L)*pow(x_0, 3)*x_1 + (1711.0L/17280.0L)*pow(x_0, 3) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/8.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) + (129.0L/320.0L)*pow(x_0, 2)*x_1 + (1567.0L/15360.0L)*pow(x_0, 2) + (7.0L/600.0L)*x_0*pow(x_1, 5) + (7.0L/72.0L)*x_0*pow(x_1, 4) + (133.0L/360.0L)*x_0*pow(x_1, 3) + (961.0L/1440.0L)*x_0*pow(x_1, 2) + (6721.0L/11520.0L)*x_0*x_1 + (9491.0L/46080.0L)*x_0 + (2.0L/675.0L)*pow(x_1, 6) + (7.0L/160.0L)*pow(x_1, 5) + (97.0L/480.0L)*pow(x_1, 4) + (1223.0L/2880.0L)*pow(x_1, 3) + (871.0L/1920.0L)*pow(x_1, 2) + (1403.0L/3840.0L)*x_1 + 264251.0L/921600.0L;
       return __pp_r483___result;
    }
    static inline O __pp_r484__(const I &x_0, const I &x_1) {
       O __pp_r484___result;
       __pp_r484___result = -7.0L/5400.0L*pow(x_0, 6) + (19.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/3600.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/40.0L)*pow(x_0, 4)*x_1 + (23.0L/1152.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 3)*x_1 + (13.0L/864.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 2)*x_1 - 1643.0L/23040.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/45.0L*x_0*pow(x_1, 4) - 1.0L/36.0L*x_0*pow(x_1, 3) - 1.0L/144.0L*x_0*pow(x_1, 2) + (7.0L/5760.0L)*x_0*x_1 + (119.0L/57600.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 1.0L/800.0L*pow(x_1, 5) + (5.0L/288.0L)*pow(x_1, 4) + (7.0L/576.0L)*pow(x_1, 3) - 841.0L/11520.0L*pow(x_1, 2) + (77.0L/38400.0L)*x_1 + 50087.0L/276480.0L;
       return __pp_r484___result;
    }
    static inline O __pp_r485__(const I &x_0, const I &x_1) {
       O __pp_r485___result;
       __pp_r485___result = -11.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 + (1.0L/14400.0L)*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 4)*x_1 + (17.0L/768.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/96.0L)*pow(x_0, 3)*x_1 + (77.0L/3456.0L)*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) + (7.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) + (31.0L/576.0L)*pow(x_0, 2)*x_1 - 887.0L/15360.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (1.0L/180.0L)*x_0*pow(x_1, 4) + (1.0L/24.0L)*x_0*pow(x_1, 3) + (23.0L/288.0L)*x_0*pow(x_1, 2) + (71.0L/1280.0L)*x_0*x_1 + (3601.0L/230400.0L)*x_0 - 1.0L/4320.0L*pow(x_1, 6) + (71.0L/7200.0L)*pow(x_1, 5) + (5.0L/96.0L)*pow(x_1, 4) + (121.0L/1728.0L)*pow(x_1, 3) - 3.0L/160.0L*pow(x_1, 2) + (839.0L/28800.0L)*x_1 + 34433.0L/184320.0L;
       return __pp_r485___result;
    }
    static inline O __pp_r486__(const I &x_0, const I &x_1) {
       O __pp_r486___result;
       __pp_r486___result = -7.0L/5400.0L*pow(x_0, 6) + (19.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/3600.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/40.0L)*pow(x_0, 4)*x_1 + (23.0L/1152.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 3)*x_1 + (13.0L/864.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 2)*x_1 - 1643.0L/23040.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/45.0L*x_0*pow(x_1, 4) - 1.0L/36.0L*x_0*pow(x_1, 3) - 1.0L/144.0L*x_0*pow(x_1, 2) + (7.0L/5760.0L)*x_0*x_1 + (119.0L/57600.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 23.0L/2400.0L*pow(x_1, 5) - 1.0L/288.0L*pow(x_1, 4) - 5.0L/576.0L*pow(x_1, 3) - 961.0L/11520.0L*pow(x_1, 2) - 23.0L/38400.0L*x_1 + 10003.0L/55296.0L;
       return __pp_r486___result;
    }
    static inline O __pp_r487__(const I &x_0, const I &x_1) {
       O __pp_r487___result;
       __pp_r487___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 + (41.0L/14400.0L)*pow(x_0, 5) + (61.0L/1440.0L)*pow(x_0, 4)*x_1 + (83.0L/2304.0L)*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) + (11.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (31.0L/288.0L)*pow(x_0, 3)*x_1 + (205.0L/3456.0L)*pow(x_0, 3) + (17.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (95.0L/576.0L)*pow(x_0, 2)*x_1 - 101.0L/46080.0L*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) + (7.0L/360.0L)*x_0*pow(x_1, 4) + (7.0L/72.0L)*x_0*pow(x_1, 3) + (55.0L/288.0L)*x_0*pow(x_1, 2) + (1919.0L/11520.0L)*x_0*x_1 + (13841.0L/230400.0L)*x_0 + (91.0L/7200.0L)*pow(x_1, 5) + (19.0L/288.0L)*pow(x_1, 4) + (185.0L/1728.0L)*pow(x_1, 3) + (53.0L/1440.0L)*pow(x_1, 2) + (2119.0L/28800.0L)*x_1 + 111491.0L/552960.0L;
       return __pp_r487___result;
    }
    static inline O __pp_r488__(const I &x_0, const I &x_1) {
       O __pp_r488___result;
       __pp_r488___result = -43.0L/43200.0L*pow(x_0, 6) + (3.0L/400.0L)*pow(x_0, 5)*x_1 + (11.0L/2880.0L)*pow(x_0, 5) + (1.0L/360.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/96.0L)*pow(x_0, 4)*x_1 + (57.0L/1280.0L)*pow(x_0, 4) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (253.0L/1440.0L)*pow(x_0, 3)*x_1 + (1711.0L/17280.0L)*pow(x_0, 3) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/8.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) + (129.0L/320.0L)*pow(x_0, 2)*x_1 + (1567.0L/15360.0L)*pow(x_0, 2) + (7.0L/600.0L)*x_0*pow(x_1, 5) + (7.0L/72.0L)*x_0*pow(x_1, 4) + (133.0L/360.0L)*x_0*pow(x_1, 3) + (961.0L/1440.0L)*x_0*pow(x_1, 2) + (6721.0L/11520.0L)*x_0*x_1 + (9491.0L/46080.0L)*x_0 + (2.0L/675.0L)*pow(x_1, 6) + (7.0L/160.0L)*pow(x_1, 5) + (97.0L/480.0L)*pow(x_1, 4) + (1223.0L/2880.0L)*pow(x_1, 3) + (871.0L/1920.0L)*pow(x_1, 2) + (1403.0L/3840.0L)*x_1 + 264251.0L/921600.0L;
       return __pp_r488___result;
    }
    static inline O __pp_r489__(const I &x_0, const I &x_1) {
       O __pp_r489___result;
       __pp_r489___result = (1.0L/360.0L)*pow(x_0, 5)*x_1 + (13.0L/7200.0L)*pow(x_0, 5) + (1.0L/90.0L)*pow(x_0, 4)*x_1 - 1.0L/144.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) + (11.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/72.0L)*pow(x_0, 3)*x_1 - 181.0L/1728.0L*pow(x_0, 3) + (17.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (7.0L/288.0L)*pow(x_0, 2)*x_1 - 791.0L/2880.0L*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) + (7.0L/360.0L)*x_0*pow(x_1, 4) + (7.0L/72.0L)*x_0*pow(x_1, 3) + (55.0L/288.0L)*x_0*pow(x_1, 2) + (11.0L/180.0L)*x_0*x_1 - 17987.0L/115200.0L*x_0 + (91.0L/7200.0L)*pow(x_1, 5) + (19.0L/288.0L)*pow(x_1, 4) + (185.0L/1728.0L)*pow(x_1, 3) + (53.0L/1440.0L)*pow(x_1, 2) + (4831.0L/115200.0L)*x_1 + 9289.0L/69120.0L;
       return __pp_r489___result;
    }
    static inline O __pp_r490__(const I &x_0, const I &x_1) {
       O __pp_r490___result;
       __pp_r490___result = (1.0L/21600.0L)*pow(x_0, 6) + (1.0L/300.0L)*pow(x_0, 5)*x_1 + (1.0L/360.0L)*pow(x_0, 5) + (1.0L/360.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/48.0L)*pow(x_0, 4)*x_1 + (1.0L/640.0L)*pow(x_0, 4) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (59.0L/720.0L)*pow(x_0, 3)*x_1 - 281.0L/4320.0L*pow(x_0, 3) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/8.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) + (21.0L/80.0L)*pow(x_0, 2)*x_1 - 1309.0L/7680.0L*pow(x_0, 2) + (7.0L/600.0L)*x_0*pow(x_1, 5) + (7.0L/72.0L)*x_0*pow(x_1, 4) + (133.0L/360.0L)*x_0*pow(x_1, 3) + (961.0L/1440.0L)*x_0*pow(x_1, 2) + (2753.0L/5760.0L)*x_0*x_1 - 59.0L/5760.0L*x_0 + (2.0L/675.0L)*pow(x_1, 6) + (7.0L/160.0L)*pow(x_1, 5) + (97.0L/480.0L)*pow(x_1, 4) + (1223.0L/2880.0L)*pow(x_1, 3) + (871.0L/1920.0L)*pow(x_1, 2) + (2563.0L/7680.0L)*x_1 + 101143.0L/460800.0L;
       return __pp_r490___result;
    }
    static inline O __pp_r491__(const I &x_0, const I &x_1) {
       O __pp_r491___result;
       __pp_r491___result = -1.0L/1440.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 79.0L/7200.0L*pow(x_0, 5) - 1.0L/96.0L*pow(x_0, 4)*pow(x_1, 2) - 13.0L/240.0L*pow(x_0, 4)*x_1 - 17.0L/160.0L*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 19.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) - 73.0L/180.0L*pow(x_0, 3)*x_1 - 4531.0L/8640.0L*pow(x_0, 3) - 1.0L/96.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/10.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 219.0L/160.0L*pow(x_0, 2)*x_1 - 1247.0L/960.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) - 139.0L/1440.0L*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) - 4639.0L/2880.0L*x_0*pow(x_1, 2) - 28271.0L/11520.0L*x_0*x_1 - 90257.0L/57600.0L*x_0 - 1.0L/576.0L*pow(x_1, 6) - 157.0L/4800.0L*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) - 5897.0L/5760.0L*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) - 162857.0L/76800.0L*x_1 - 731563.0L/921600.0L;
       return __pp_r491___result;
    }
    static inline O __pp_r492__(const I &x_0, const I &x_1) {
       O __pp_r492___result;
       __pp_r492___result = -7.0L/10800.0L*pow(x_0, 6) - 1.0L/1200.0L*pow(x_0, 5)*x_1 - 7.0L/720.0L*pow(x_0, 5) - 11.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 1.0L/24.0L*pow(x_0, 4)*x_1 - 59.0L/640.0L*pow(x_0, 4) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/18.0L*pow(x_0, 3)*pow(x_1, 2) - 211.0L/720.0L*pow(x_0, 3)*x_1 - 1901.0L/4320.0L*pow(x_0, 3) + (1.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 17.0L/80.0L*pow(x_0, 2)*pow(x_1, 2) - 69.0L/80.0L*pow(x_0, 2)*x_1 - 7789.0L/7680.0L*pow(x_0, 2) + (3.0L/400.0L)*x_0*pow(x_1, 5) + (5.0L/144.0L)*x_0*pow(x_1, 4) - 1.0L/180.0L*x_0*pow(x_1, 3) - 659.0L/1440.0L*x_0*pow(x_1, 2) - 6967.0L/5760.0L*x_0*x_1 - 5891.0L/5760.0L*x_0 + (49.0L/21600.0L)*pow(x_1, 6) + (1.0L/32.0L)*pow(x_1, 5) + (13.0L/120.0L)*pow(x_1, 4) + (143.0L/2880.0L)*pow(x_1, 3) - 749.0L/1920.0L*pow(x_1, 2) - 5213.0L/7680.0L*x_1 - 132137.0L/460800.0L;
       return __pp_r492___result;
    }
    static inline O __pp_r493__(const I &x_0, const I &x_1) {
       O __pp_r493___result;
       __pp_r493___result = -1.0L/1440.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 79.0L/7200.0L*pow(x_0, 5) - 1.0L/96.0L*pow(x_0, 4)*pow(x_1, 2) - 13.0L/240.0L*pow(x_0, 4)*x_1 - 17.0L/160.0L*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 19.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) - 73.0L/180.0L*pow(x_0, 3)*x_1 - 4531.0L/8640.0L*pow(x_0, 3) - 1.0L/96.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/10.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 219.0L/160.0L*pow(x_0, 2)*x_1 - 1247.0L/960.0L*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) - 47.0L/720.0L*x_0*pow(x_1, 4) - 41.0L/90.0L*x_0*pow(x_1, 3) - 2117.0L/1440.0L*x_0*pow(x_1, 2) - 1691.0L/720.0L*x_0*x_1 - 176869.0L/115200.0L*x_0 - 1.0L/1440.0L*pow(x_1, 6) - 7.0L/800.0L*pow(x_1, 5) - 7.0L/60.0L*pow(x_1, 4) - 1801.0L/2880.0L*pow(x_1, 3) - 367.0L/240.0L*pow(x_1, 2) - 65431.0L/38400.0L*x_1 - 77321.0L/115200.0L;
       return __pp_r493___result;
    }
    static inline O __pp_r494__(const I &x_0, const I &x_1) {
       O __pp_r494___result;
       __pp_r494___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (31.0L/2400.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/60.0L*pow(x_0, 4)*x_1 + (1.0L/36.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 3)*x_1 - 3.0L/64.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/16.0L*pow(x_0, 2)*x_1 - 2539.0L/11520.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) + (11.0L/480.0L)*x_0*pow(x_1, 4) + (23.0L/288.0L)*x_0*pow(x_1, 3) + (15.0L/64.0L)*x_0*pow(x_1, 2) + (79.0L/11520.0L)*x_0*x_1 - 2477.0L/19200.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (59.0L/4800.0L)*pow(x_1, 5) + (157.0L/2304.0L)*pow(x_1, 4) + (115.0L/1152.0L)*pow(x_1, 3) + (2321.0L/46080.0L)*pow(x_1, 2) + (2179.0L/76800.0L)*x_1 + 77437.0L/552960.0L;
       return __pp_r494___result;
    }
    static inline O __pp_r495__(const I &x_0, const I &x_1) {
       O __pp_r495___result;
       __pp_r495___result = (11.0L/7200.0L)*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (1.0L/72.0L)*pow(x_0, 5) + (1.0L/120.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/144.0L*pow(x_0, 4)*x_1 + (209.0L/5760.0L)*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/80.0L)*pow(x_0, 3)*x_1 - 31.0L/4320.0L*pow(x_0, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/9.0L)*pow(x_0, 2)*pow(x_1, 3) + (193.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (253.0L/1440.0L)*pow(x_0, 2)*x_1 - 2677.0L/23040.0L*pow(x_0, 2) + (41.0L/3600.0L)*x_0*pow(x_1, 5) + (29.0L/288.0L)*x_0*pow(x_1, 4) + (169.0L/480.0L)*x_0*pow(x_1, 3) + (2047.0L/2880.0L)*x_0*pow(x_1, 2) + (1627.0L/3840.0L)*x_0*x_1 + (389.0L/23040.0L)*x_0 + (43.0L/14400.0L)*pow(x_1, 6) + (25.0L/576.0L)*pow(x_1, 5) + (2353.0L/11520.0L)*pow(x_1, 4) + (7213.0L/17280.0L)*pow(x_1, 3) + (21529.0L/46080.0L)*pow(x_1, 2) + (14753.0L/46080.0L)*x_1 + 622483.0L/2764800.0L;
       return __pp_r495___result;
    }
    static inline O __pp_r496__(const I &x_0, const I &x_1) {
       O __pp_r496___result;
       __pp_r496___result = (17.0L/21600.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 + (1.0L/7200.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 59.0L/720.0L*pow(x_0, 4)*x_1 - 103.0L/1440.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 7.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) - 19.0L/40.0L*pow(x_0, 3)*x_1 - 4031.0L/8640.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 41.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 131.0L/90.0L*pow(x_0, 2)*x_1 - 14339.0L/11520.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) - 67.0L/720.0L*x_0*pow(x_1, 4) - 17.0L/30.0L*x_0*pow(x_1, 3) - 2257.0L/1440.0L*x_0*pow(x_1, 2) - 301.0L/120.0L*x_0*x_1 - 177389.0L/115200.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 119.0L/3600.0L*pow(x_1, 5) - 373.0L/1440.0L*pow(x_1, 4) - 2227.0L/2160.0L*pow(x_1, 3) - 24041.0L/11520.0L*pow(x_1, 2) - 30731.0L/14400.0L*x_1 - 272383.0L/345600.0L;
       return __pp_r496___result;
    }
    static inline O __pp_r497__(const I &x_0, const I &x_1) {
       O __pp_r497___result;
       __pp_r497___result = (1.0L/1200.0L)*pow(x_0, 6) - 19.0L/3600.0L*pow(x_0, 5)*x_1 + (1.0L/720.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 5.0L/72.0L*pow(x_0, 4)*x_1 - 331.0L/5760.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/36.0L*pow(x_0, 3)*pow(x_1, 2) - 29.0L/80.0L*pow(x_0, 3)*x_1 - 1651.0L/4320.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 1367.0L/1440.0L*pow(x_0, 2)*x_1 - 22117.0L/23040.0L*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) + (11.0L/288.0L)*x_0*pow(x_1, 4) - 11.0L/480.0L*x_0*pow(x_1, 3) - 1193.0L/2880.0L*x_0*pow(x_1, 2) - 4853.0L/3840.0L*x_0*x_1 - 22939.0L/23040.0L*x_0 + (11.0L/4800.0L)*pow(x_1, 6) + (89.0L/2880.0L)*pow(x_1, 5) + (1273.0L/11520.0L)*pow(x_1, 4) + (733.0L/17280.0L)*pow(x_1, 3) - 17351.0L/46080.0L*pow(x_1, 2) - 31903.0L/46080.0L*x_1 - 777197.0L/2764800.0L;
       return __pp_r497___result;
    }
    static inline O __pp_r498__(const I &x_0, const I &x_1) {
       O __pp_r498___result;
       __pp_r498___result = (17.0L/21600.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 + (1.0L/7200.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 59.0L/720.0L*pow(x_0, 4)*x_1 - 103.0L/1440.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 7.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) - 19.0L/40.0L*pow(x_0, 3)*x_1 - 4031.0L/8640.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 41.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 131.0L/90.0L*pow(x_0, 2)*x_1 - 14339.0L/11520.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) - 89.0L/1440.0L*x_0*pow(x_1, 4) - 227.0L/480.0L*x_0*pow(x_1, 3) - 4109.0L/2880.0L*x_0*pow(x_1, 2) - 9227.0L/3840.0L*x_0*x_1 - 10859.0L/7200.0L*x_0 - 29.0L/43200.0L*pow(x_1, 6) - 131.0L/14400.0L*pow(x_1, 5) - 1319.0L/11520.0L*pow(x_1, 4) - 10931.0L/17280.0L*pow(x_1, 3) - 69839.0L/46080.0L*pow(x_1, 2) - 395711.0L/230400.0L*x_1 - 1840079.0L/2764800.0L;
       return __pp_r498___result;
    }
    static inline O __pp_r499__(const I &x_0, const I &x_1) {
       O __pp_r499___result;
       __pp_r499___result = -1.0L/1440.0L*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 - 79.0L/7200.0L*pow(x_0, 5) - 1.0L/96.0L*pow(x_0, 4)*pow(x_1, 2) - 13.0L/240.0L*pow(x_0, 4)*x_1 - 17.0L/160.0L*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) - 19.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) - 73.0L/180.0L*pow(x_0, 3)*x_1 - 4531.0L/8640.0L*pow(x_0, 3) - 1.0L/96.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/10.0L*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) - 219.0L/160.0L*pow(x_0, 2)*x_1 - 1247.0L/960.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) - 139.0L/1440.0L*x_0*pow(x_1, 4) - 791.0L/1440.0L*x_0*pow(x_1, 3) - 4639.0L/2880.0L*x_0*pow(x_1, 2) - 28271.0L/11520.0L*x_0*x_1 - 90257.0L/57600.0L*x_0 - 1.0L/576.0L*pow(x_1, 6) - 157.0L/4800.0L*pow(x_1, 5) - 1003.0L/3840.0L*pow(x_1, 4) - 5897.0L/5760.0L*pow(x_1, 3) - 32263.0L/15360.0L*pow(x_1, 2) - 162857.0L/76800.0L*x_1 - 731563.0L/921600.0L;
       return __pp_r499___result;
    }
    static inline O __pp_r500__(const I &x_0, const I &x_1) {
       O __pp_r500___result;
       __pp_r500___result = -1.0L/4320.0L*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 + (1.0L/7200.0L)*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/720.0L)*pow(x_0, 4)*x_1 + (7.0L/1440.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (7.0L/180.0L)*pow(x_0, 3)*x_1 + (589.0L/8640.0L)*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) + (589.0L/1440.0L)*pow(x_0, 2)*x_1 + (1379.0L/2880.0L)*pow(x_0, 2) - 1.0L/360.0L*x_0*pow(x_1, 5) - 59.0L/1440.0L*x_0*pow(x_1, 4) - 151.0L/1440.0L*x_0*pow(x_1, 3) + (481.0L/2880.0L)*x_0*pow(x_1, 2) + (12689.0L/11520.0L)*x_0*x_1 + (73583.0L/57600.0L)*x_0 - 11.0L/8640.0L*pow(x_1, 6) - 311.0L/14400.0L*pow(x_1, 5) - 1729.0L/11520.0L*pow(x_1, 4) - 7451.0L/17280.0L*pow(x_1, 3) - 14869.0L/46080.0L*pow(x_1, 2) + (166789.0L/230400.0L)*x_1 + 3048191.0L/2764800.0L;
       return __pp_r500___result;
    }
    static inline O __pp_r501__(const I &x_0, const I &x_1) {
       O __pp_r501___result;
       __pp_r501___result = -1.0L/4800.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 + (1.0L/960.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 4)*x_1 + (5.0L/256.0L)*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (5.0L/32.0L)*pow(x_0, 3)*x_1 + (25.0L/128.0L)*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/12.0L)*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) + (75.0L/64.0L)*pow(x_0, 2)*x_1 + (1125.0L/1024.0L)*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) + (1.0L/32.0L)*x_0*pow(x_1, 4) + (35.0L/96.0L)*x_0*pow(x_1, 3) + (325.0L/192.0L)*x_0*pow(x_1, 2) + (1375.0L/384.0L)*x_0*x_1 + (8875.0L/3072.0L)*x_0 + (1.0L/4800.0L)*pow(x_1, 6) + (7.0L/960.0L)*pow(x_1, 5) + (65.0L/768.0L)*pow(x_1, 4) + (75.0L/128.0L)*pow(x_1, 3) + (6625.0L/3072.0L)*pow(x_1, 2) + (12125.0L/3072.0L)*x_1 + 4375.0L/1536.0L;
       return __pp_r501___result;
    }
    static inline O __pp_r502__(const I &x_0, const I &x_1) {
       O __pp_r502___result;
       __pp_r502___result = (17.0L/21600.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 + (1.0L/7200.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 59.0L/720.0L*pow(x_0, 4)*x_1 - 103.0L/1440.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 7.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) - 19.0L/40.0L*pow(x_0, 3)*x_1 - 4031.0L/8640.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 41.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 131.0L/90.0L*pow(x_0, 2)*x_1 - 14339.0L/11520.0L*pow(x_0, 2) - 7.0L/1200.0L*x_0*pow(x_1, 5) - 67.0L/720.0L*x_0*pow(x_1, 4) - 17.0L/30.0L*x_0*pow(x_1, 3) - 2257.0L/1440.0L*x_0*pow(x_1, 2) - 301.0L/120.0L*x_0*x_1 - 177389.0L/115200.0L*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 119.0L/3600.0L*pow(x_1, 5) - 373.0L/1440.0L*pow(x_1, 4) - 2227.0L/2160.0L*pow(x_1, 3) - 24041.0L/11520.0L*pow(x_1, 2) - 30731.0L/14400.0L*x_1 - 272383.0L/345600.0L;
       return __pp_r502___result;
    }
    static inline O __pp_r503__(const I &x_0, const I &x_1) {
       O __pp_r503___result;
       __pp_r503___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 + (9.0L/800.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 19.0L/720.0L*pow(x_0, 4)*x_1 + (19.0L/480.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) - 11.0L/360.0L*pow(x_0, 3)*x_1 + (121.0L/960.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (29.0L/90.0L)*pow(x_0, 2)*x_1 + (2047.0L/3840.0L)*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 3.0L/80.0L*x_0*pow(x_1, 4) - 11.0L/90.0L*x_0*pow(x_1, 3) + (101.0L/480.0L)*x_0*pow(x_1, 2) + (377.0L/360.0L)*x_0*x_1 + (16699.0L/12800.0L)*x_0 - 1.0L/800.0L*pow(x_1, 6) - 79.0L/3600.0L*pow(x_1, 5) - 71.0L/480.0L*pow(x_1, 4) - 947.0L/2160.0L*pow(x_1, 3) - 1187.0L/3840.0L*pow(x_1, 2) + (10229.0L/14400.0L)*x_1 + 42553.0L/38400.0L;
       return __pp_r503___result;
    }
    static inline O __pp_r504__(const I &x_0, const I &x_1) {
       O __pp_r504___result;
       __pp_r504___result = (11.0L/8640.0L)*pow(x_0, 6) - 1.0L/360.0L*pow(x_0, 5)*x_1 + (7.0L/576.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/288.0L*pow(x_0, 4)*x_1 + (125.0L/2304.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (25.0L/288.0L)*pow(x_0, 3)*x_1 + (875.0L/3456.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (625.0L/576.0L)*pow(x_0, 2)*x_1 + (10625.0L/9216.0L)*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (5.0L/144.0L)*x_0*pow(x_1, 4) + (25.0L/72.0L)*x_0*pow(x_1, 3) + (125.0L/72.0L)*x_0*pow(x_1, 2) + (8125.0L/2304.0L)*x_0*x_1 + (26875.0L/9216.0L)*x_0 + (1.0L/4320.0L)*pow(x_1, 6) + (1.0L/144.0L)*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) + (125.0L/216.0L)*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) + (18125.0L/4608.0L)*x_1 + 315625.0L/110592.0L;
       return __pp_r504___result;
    }
    static inline O __pp_r505__(const I &x_0, const I &x_1) {
       O __pp_r505___result;
       __pp_r505___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 + (29.0L/4800.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (37.0L/1440.0L)*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 + (289.0L/640.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (2803.0L/2880.0L)*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 3.0L/80.0L*x_0*pow(x_1, 4) - 11.0L/90.0L*x_0*pow(x_1, 3) + (101.0L/480.0L)*x_0*pow(x_1, 2) + (21439.0L/11520.0L)*x_0*x_1 + (209569.0L/76800.0L)*x_0 - 1.0L/800.0L*pow(x_1, 6) - 79.0L/3600.0L*pow(x_1, 5) - 71.0L/480.0L*pow(x_1, 4) - 947.0L/2160.0L*pow(x_1, 3) - 1187.0L/3840.0L*pow(x_1, 2) + (128707.0L/115200.0L)*x_1 + 574799.0L/307200.0L;
       return __pp_r505___result;
    }
    static inline O __pp_r506__(const I &x_0, const I &x_1) {
       O __pp_r506___result;
       __pp_r506___result = (1.0L/4320.0L)*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 + (1.0L/144.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (25.0L/72.0L)*pow(x_0, 3)*x_1 + (125.0L/216.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (125.0L/72.0L)*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (5.0L/144.0L)*x_0*pow(x_1, 4) + (25.0L/72.0L)*x_0*pow(x_1, 3) + (125.0L/72.0L)*x_0*pow(x_1, 2) + (625.0L/144.0L)*x_0*x_1 + (625.0L/144.0L)*x_0 + (1.0L/4320.0L)*pow(x_1, 6) + (1.0L/144.0L)*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) + (125.0L/216.0L)*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) + (625.0L/144.0L)*x_1 + 3125.0L/864.0L;
       return __pp_r506___result;
    }
    static inline O __pp_r507__(const I &x_0, const I &x_1) {
       O __pp_r507___result;
       __pp_r507___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 + (29.0L/4800.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (37.0L/1440.0L)*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 + (289.0L/640.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (2803.0L/2880.0L)*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 3.0L/80.0L*x_0*pow(x_1, 4) - 11.0L/90.0L*x_0*pow(x_1, 3) + (101.0L/480.0L)*x_0*pow(x_1, 2) + (21439.0L/11520.0L)*x_0*x_1 + (209569.0L/76800.0L)*x_0 - 1.0L/800.0L*pow(x_1, 6) - 79.0L/3600.0L*pow(x_1, 5) - 71.0L/480.0L*pow(x_1, 4) - 947.0L/2160.0L*pow(x_1, 3) - 1187.0L/3840.0L*pow(x_1, 2) + (128707.0L/115200.0L)*x_1 + 574799.0L/307200.0L;
       return __pp_r507___result;
    }
    static inline O __pp_r508__(const I &x_0, const I &x_1) {
       O __pp_r508___result;
       __pp_r508___result = (1.0L/4320.0L)*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 + (1.0L/144.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (25.0L/72.0L)*pow(x_0, 3)*x_1 + (125.0L/216.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (125.0L/72.0L)*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (5.0L/144.0L)*x_0*pow(x_1, 4) + (25.0L/72.0L)*x_0*pow(x_1, 3) + (125.0L/72.0L)*x_0*pow(x_1, 2) + (625.0L/144.0L)*x_0*x_1 + (625.0L/144.0L)*x_0 + (1.0L/4320.0L)*pow(x_1, 6) + (1.0L/144.0L)*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) + (125.0L/216.0L)*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) + (625.0L/144.0L)*x_1 + 3125.0L/864.0L;
       return __pp_r508___result;
    }
    static inline O __pp_r509__(const I &x_0, const I &x_1) {
       O __pp_r509___result;
       __pp_r509___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (31.0L/2400.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/60.0L*pow(x_0, 4)*x_1 + (1.0L/36.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 3)*x_1 - 3.0L/64.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/16.0L*pow(x_0, 2)*x_1 - 2539.0L/11520.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) + (11.0L/480.0L)*x_0*pow(x_1, 4) + (23.0L/288.0L)*x_0*pow(x_1, 3) + (15.0L/64.0L)*x_0*pow(x_1, 2) + (79.0L/11520.0L)*x_0*x_1 - 2477.0L/19200.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (59.0L/4800.0L)*pow(x_1, 5) + (157.0L/2304.0L)*pow(x_1, 4) + (115.0L/1152.0L)*pow(x_1, 3) + (2321.0L/46080.0L)*pow(x_1, 2) + (2179.0L/76800.0L)*x_1 + 77437.0L/552960.0L;
       return __pp_r509___result;
    }
    static inline O __pp_r510__(const I &x_0, const I &x_1) {
       O __pp_r510___result;
       __pp_r510___result = (11.0L/7200.0L)*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (1.0L/72.0L)*pow(x_0, 5) + (1.0L/120.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/144.0L*pow(x_0, 4)*x_1 + (209.0L/5760.0L)*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/80.0L)*pow(x_0, 3)*x_1 - 31.0L/4320.0L*pow(x_0, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/9.0L)*pow(x_0, 2)*pow(x_1, 3) + (193.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (253.0L/1440.0L)*pow(x_0, 2)*x_1 - 2677.0L/23040.0L*pow(x_0, 2) + (41.0L/3600.0L)*x_0*pow(x_1, 5) + (29.0L/288.0L)*x_0*pow(x_1, 4) + (169.0L/480.0L)*x_0*pow(x_1, 3) + (2047.0L/2880.0L)*x_0*pow(x_1, 2) + (1627.0L/3840.0L)*x_0*x_1 + (389.0L/23040.0L)*x_0 + (43.0L/14400.0L)*pow(x_1, 6) + (25.0L/576.0L)*pow(x_1, 5) + (2353.0L/11520.0L)*pow(x_1, 4) + (7213.0L/17280.0L)*pow(x_1, 3) + (21529.0L/46080.0L)*pow(x_1, 2) + (14753.0L/46080.0L)*x_1 + 622483.0L/2764800.0L;
       return __pp_r510___result;
    }
    static inline O __pp_r511__(const I &x_0, const I &x_1) {
       O __pp_r511___result;
       __pp_r511___result = (1.0L/1200.0L)*pow(x_0, 6) - 19.0L/3600.0L*pow(x_0, 5)*x_1 + (1.0L/720.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 5.0L/72.0L*pow(x_0, 4)*x_1 - 331.0L/5760.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/36.0L*pow(x_0, 3)*pow(x_1, 2) - 29.0L/80.0L*pow(x_0, 3)*x_1 - 1651.0L/4320.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 1367.0L/1440.0L*pow(x_0, 2)*x_1 - 22117.0L/23040.0L*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) + (11.0L/288.0L)*x_0*pow(x_1, 4) - 11.0L/480.0L*x_0*pow(x_1, 3) - 1193.0L/2880.0L*x_0*pow(x_1, 2) - 4853.0L/3840.0L*x_0*x_1 - 22939.0L/23040.0L*x_0 + (11.0L/4800.0L)*pow(x_1, 6) + (89.0L/2880.0L)*pow(x_1, 5) + (1273.0L/11520.0L)*pow(x_1, 4) + (733.0L/17280.0L)*pow(x_1, 3) - 17351.0L/46080.0L*pow(x_1, 2) - 31903.0L/46080.0L*x_1 - 777197.0L/2764800.0L;
       return __pp_r511___result;
    }
    static inline O __pp_r512__(const I &x_0, const I &x_1) {
       O __pp_r512___result;
       __pp_r512___result = (17.0L/21600.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 + (1.0L/7200.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 59.0L/720.0L*pow(x_0, 4)*x_1 - 103.0L/1440.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 7.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) - 19.0L/40.0L*pow(x_0, 3)*x_1 - 4031.0L/8640.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 41.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 131.0L/90.0L*pow(x_0, 2)*x_1 - 14339.0L/11520.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) - 89.0L/1440.0L*x_0*pow(x_1, 4) - 227.0L/480.0L*x_0*pow(x_1, 3) - 4109.0L/2880.0L*x_0*pow(x_1, 2) - 9227.0L/3840.0L*x_0*x_1 - 10859.0L/7200.0L*x_0 - 29.0L/43200.0L*pow(x_1, 6) - 131.0L/14400.0L*pow(x_1, 5) - 1319.0L/11520.0L*pow(x_1, 4) - 10931.0L/17280.0L*pow(x_1, 3) - 69839.0L/46080.0L*pow(x_1, 2) - 395711.0L/230400.0L*x_1 - 1840079.0L/2764800.0L;
       return __pp_r512___result;
    }
    static inline O __pp_r513__(const I &x_0, const I &x_1) {
       O __pp_r513___result;
       __pp_r513___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (31.0L/2400.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/60.0L*pow(x_0, 4)*x_1 + (1.0L/36.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 3)*x_1 - 3.0L/64.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/16.0L*pow(x_0, 2)*x_1 - 2539.0L/11520.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) + (11.0L/480.0L)*x_0*pow(x_1, 4) + (23.0L/288.0L)*x_0*pow(x_1, 3) + (15.0L/64.0L)*x_0*pow(x_1, 2) + (79.0L/11520.0L)*x_0*x_1 - 2477.0L/19200.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (59.0L/4800.0L)*pow(x_1, 5) + (157.0L/2304.0L)*pow(x_1, 4) + (115.0L/1152.0L)*pow(x_1, 3) + (2321.0L/46080.0L)*pow(x_1, 2) + (2179.0L/76800.0L)*x_1 + 77437.0L/552960.0L;
       return __pp_r513___result;
    }
    static inline O __pp_r514__(const I &x_0, const I &x_1) {
       O __pp_r514___result;
       __pp_r514___result = (11.0L/7200.0L)*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (1.0L/72.0L)*pow(x_0, 5) + (1.0L/120.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/144.0L*pow(x_0, 4)*x_1 + (209.0L/5760.0L)*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (1.0L/80.0L)*pow(x_0, 3)*x_1 - 31.0L/4320.0L*pow(x_0, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/9.0L)*pow(x_0, 2)*pow(x_1, 3) + (193.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (253.0L/1440.0L)*pow(x_0, 2)*x_1 - 2677.0L/23040.0L*pow(x_0, 2) + (41.0L/3600.0L)*x_0*pow(x_1, 5) + (29.0L/288.0L)*x_0*pow(x_1, 4) + (169.0L/480.0L)*x_0*pow(x_1, 3) + (2047.0L/2880.0L)*x_0*pow(x_1, 2) + (1627.0L/3840.0L)*x_0*x_1 + (389.0L/23040.0L)*x_0 + (43.0L/14400.0L)*pow(x_1, 6) + (25.0L/576.0L)*pow(x_1, 5) + (2353.0L/11520.0L)*pow(x_1, 4) + (7213.0L/17280.0L)*pow(x_1, 3) + (21529.0L/46080.0L)*pow(x_1, 2) + (14753.0L/46080.0L)*x_1 + 622483.0L/2764800.0L;
       return __pp_r514___result;
    }
    static inline O __pp_r515__(const I &x_0, const I &x_1) {
       O __pp_r515___result;
       __pp_r515___result = (1.0L/1200.0L)*pow(x_0, 6) - 19.0L/3600.0L*pow(x_0, 5)*x_1 + (1.0L/720.0L)*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 5.0L/72.0L*pow(x_0, 4)*x_1 - 331.0L/5760.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/36.0L*pow(x_0, 3)*pow(x_1, 2) - 29.0L/80.0L*pow(x_0, 3)*x_1 - 1651.0L/4320.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 1367.0L/1440.0L*pow(x_0, 2)*x_1 - 22117.0L/23040.0L*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) + (11.0L/288.0L)*x_0*pow(x_1, 4) - 11.0L/480.0L*x_0*pow(x_1, 3) - 1193.0L/2880.0L*x_0*pow(x_1, 2) - 4853.0L/3840.0L*x_0*x_1 - 22939.0L/23040.0L*x_0 + (11.0L/4800.0L)*pow(x_1, 6) + (89.0L/2880.0L)*pow(x_1, 5) + (1273.0L/11520.0L)*pow(x_1, 4) + (733.0L/17280.0L)*pow(x_1, 3) - 17351.0L/46080.0L*pow(x_1, 2) - 31903.0L/46080.0L*x_1 - 777197.0L/2764800.0L;
       return __pp_r515___result;
    }
    static inline O __pp_r516__(const I &x_0, const I &x_1) {
       O __pp_r516___result;
       __pp_r516___result = (17.0L/21600.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 + (1.0L/7200.0L)*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 59.0L/720.0L*pow(x_0, 4)*x_1 - 103.0L/1440.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 7.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) - 19.0L/40.0L*pow(x_0, 3)*x_1 - 4031.0L/8640.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 41.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 131.0L/90.0L*pow(x_0, 2)*x_1 - 14339.0L/11520.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) - 89.0L/1440.0L*x_0*pow(x_1, 4) - 227.0L/480.0L*x_0*pow(x_1, 3) - 4109.0L/2880.0L*x_0*pow(x_1, 2) - 9227.0L/3840.0L*x_0*x_1 - 10859.0L/7200.0L*x_0 - 29.0L/43200.0L*pow(x_1, 6) - 131.0L/14400.0L*pow(x_1, 5) - 1319.0L/11520.0L*pow(x_1, 4) - 10931.0L/17280.0L*pow(x_1, 3) - 69839.0L/46080.0L*pow(x_1, 2) - 395711.0L/230400.0L*x_1 - 1840079.0L/2764800.0L;
       return __pp_r516___result;
    }
    static inline O __pp_r517__(const I &x_0, const I &x_1) {
       O __pp_r517___result;
       __pp_r517___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 11.0L/2880.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 5.0L/288.0L*pow(x_0, 4)*x_1 - 287.0L/11520.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/36.0L*pow(x_0, 3)*pow(x_1, 2) - 49.0L/480.0L*pow(x_0, 3)*x_1 - 979.0L/17280.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 859.0L/2880.0L*pow(x_0, 2)*x_1 + (2641.0L/46080.0L)*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) + (11.0L/288.0L)*x_0*pow(x_1, 4) - 11.0L/480.0L*x_0*pow(x_1, 3) - 1193.0L/2880.0L*x_0*pow(x_1, 2) - 9.0L/20.0L*x_0*x_1 + (19747.0L/46080.0L)*x_0 + (11.0L/4800.0L)*pow(x_1, 6) + (89.0L/2880.0L)*pow(x_1, 5) + (1273.0L/11520.0L)*pow(x_1, 4) + (733.0L/17280.0L)*pow(x_1, 3) - 17351.0L/46080.0L*pow(x_1, 2) - 13153.0L/46080.0L*x_1 + 666089.0L/1382400.0L;
       return __pp_r517___result;
    }
    static inline O __pp_r518__(const I &x_0, const I &x_1) {
       O __pp_r518___result;
       __pp_r518___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 73.0L/14400.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 43.0L/1440.0L*pow(x_0, 4)*x_1 - 449.0L/11520.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 7.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) - 103.0L/480.0L*pow(x_0, 3)*x_1 - 2437.0L/17280.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 41.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 2317.0L/2880.0L*pow(x_0, 2)*x_1 - 10481.0L/46080.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) - 89.0L/1440.0L*x_0*pow(x_1, 4) - 227.0L/480.0L*x_0*pow(x_1, 3) - 4109.0L/2880.0L*x_0*pow(x_1, 2) - 1017.0L/640.0L*x_0*x_1 - 19363.0L/230400.0L*x_0 - 29.0L/43200.0L*pow(x_1, 6) - 131.0L/14400.0L*pow(x_1, 5) - 1319.0L/11520.0L*pow(x_1, 4) - 10931.0L/17280.0L*pow(x_1, 3) - 69839.0L/46080.0L*pow(x_1, 2) - 301961.0L/230400.0L*x_1 + 16831.0L/172800.0L;
       return __pp_r518___result;
    }
    static inline O __pp_r519__(const I &x_0, const I &x_1) {
       O __pp_r519___result;
       __pp_r519___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 + (29.0L/4800.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (37.0L/1440.0L)*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 + (289.0L/640.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (2803.0L/2880.0L)*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 1.0L/160.0L*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) + (337.0L/960.0L)*x_0*pow(x_1, 2) + (11327.0L/5760.0L)*x_0*x_1 + (211999.0L/76800.0L)*x_0 - 1.0L/4800.0L*pow(x_1, 6) + (29.0L/14400.0L)*pow(x_1, 5) - 13.0L/3840.0L*pow(x_1, 4) - 691.0L/17280.0L*pow(x_1, 3) + (4027.0L/15360.0L)*pow(x_1, 2) + (353399.0L/230400.0L)*x_1 + 38279.0L/19200.0L;
       return __pp_r519___result;
    }
    static inline O __pp_r520__(const I &x_0, const I &x_1) {
       O __pp_r520___result;
       __pp_r520___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 + (29.0L/4800.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (37.0L/1440.0L)*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 + (289.0L/640.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (2803.0L/2880.0L)*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 3.0L/80.0L*x_0*pow(x_1, 4) - 11.0L/90.0L*x_0*pow(x_1, 3) + (101.0L/480.0L)*x_0*pow(x_1, 2) + (21439.0L/11520.0L)*x_0*x_1 + (209569.0L/76800.0L)*x_0 - 1.0L/800.0L*pow(x_1, 6) - 79.0L/3600.0L*pow(x_1, 5) - 71.0L/480.0L*pow(x_1, 4) - 947.0L/2160.0L*pow(x_1, 3) - 1187.0L/3840.0L*pow(x_1, 2) + (128707.0L/115200.0L)*x_1 + 574799.0L/307200.0L;
       return __pp_r520___result;
    }
    static inline O __pp_r521__(const I &x_0, const I &x_1) {
       O __pp_r521___result;
       __pp_r521___result = (1.0L/4320.0L)*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 + (1.0L/144.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (25.0L/72.0L)*pow(x_0, 3)*x_1 + (125.0L/216.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (125.0L/72.0L)*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (5.0L/144.0L)*x_0*pow(x_1, 4) + (25.0L/72.0L)*x_0*pow(x_1, 3) + (125.0L/72.0L)*x_0*pow(x_1, 2) + (625.0L/144.0L)*x_0*x_1 + (625.0L/144.0L)*x_0 + (1.0L/4320.0L)*pow(x_1, 6) + (1.0L/144.0L)*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) + (125.0L/216.0L)*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) + (625.0L/144.0L)*x_1 + 3125.0L/864.0L;
       return __pp_r521___result;
    }
    static inline O __pp_r522__(const I &x_0, const I &x_1) {
       O __pp_r522___result;
       __pp_r522___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 + (29.0L/4800.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (37.0L/1440.0L)*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 + (289.0L/640.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (2803.0L/2880.0L)*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 1.0L/160.0L*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) + (337.0L/960.0L)*x_0*pow(x_1, 2) + (11327.0L/5760.0L)*x_0*x_1 + (211999.0L/76800.0L)*x_0 - 1.0L/4800.0L*pow(x_1, 6) + (29.0L/14400.0L)*pow(x_1, 5) - 13.0L/3840.0L*pow(x_1, 4) - 691.0L/17280.0L*pow(x_1, 3) + (4027.0L/15360.0L)*pow(x_1, 2) + (353399.0L/230400.0L)*x_1 + 38279.0L/19200.0L;
       return __pp_r522___result;
    }
    static inline O __pp_r523__(const I &x_0, const I &x_1) {
       O __pp_r523___result;
       __pp_r523___result = -1.0L/4800.0L*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 - 11.0L/2880.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) - 5.0L/288.0L*pow(x_0, 4)*x_1 - 287.0L/11520.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/36.0L*pow(x_0, 3)*pow(x_1, 2) - 49.0L/480.0L*pow(x_0, 3)*x_1 - 979.0L/17280.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) - 77.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 859.0L/2880.0L*pow(x_0, 2)*x_1 + (2641.0L/46080.0L)*pow(x_0, 2) + (13.0L/1800.0L)*x_0*pow(x_1, 5) + (11.0L/288.0L)*x_0*pow(x_1, 4) - 11.0L/480.0L*x_0*pow(x_1, 3) - 1193.0L/2880.0L*x_0*pow(x_1, 2) - 9.0L/20.0L*x_0*x_1 + (19747.0L/46080.0L)*x_0 + (11.0L/4800.0L)*pow(x_1, 6) + (89.0L/2880.0L)*pow(x_1, 5) + (1273.0L/11520.0L)*pow(x_1, 4) + (733.0L/17280.0L)*pow(x_1, 3) - 17351.0L/46080.0L*pow(x_1, 2) - 13153.0L/46080.0L*x_1 + 666089.0L/1382400.0L;
       return __pp_r523___result;
    }
    static inline O __pp_r524__(const I &x_0, const I &x_1) {
       O __pp_r524___result;
       __pp_r524___result = -11.0L/43200.0L*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 - 73.0L/14400.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) - 43.0L/1440.0L*pow(x_0, 4)*x_1 - 449.0L/11520.0L*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) - 7.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) - 103.0L/480.0L*pow(x_0, 3)*x_1 - 2437.0L/17280.0L*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 41.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) - 239.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 2317.0L/2880.0L*pow(x_0, 2)*x_1 - 10481.0L/46080.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) - 89.0L/1440.0L*x_0*pow(x_1, 4) - 227.0L/480.0L*x_0*pow(x_1, 3) - 4109.0L/2880.0L*x_0*pow(x_1, 2) - 1017.0L/640.0L*x_0*x_1 - 19363.0L/230400.0L*x_0 - 29.0L/43200.0L*pow(x_1, 6) - 131.0L/14400.0L*pow(x_1, 5) - 1319.0L/11520.0L*pow(x_1, 4) - 10931.0L/17280.0L*pow(x_1, 3) - 69839.0L/46080.0L*pow(x_1, 2) - 301961.0L/230400.0L*x_1 + 16831.0L/172800.0L;
       return __pp_r524___result;
    }
    static inline O __pp_r525__(const I &x_0, const I &x_1) {
       O __pp_r525___result;
       __pp_r525___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 + (29.0L/4800.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (37.0L/1440.0L)*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 + (289.0L/640.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (2803.0L/2880.0L)*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 1.0L/160.0L*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) + (337.0L/960.0L)*x_0*pow(x_1, 2) + (11327.0L/5760.0L)*x_0*x_1 + (211999.0L/76800.0L)*x_0 - 1.0L/4800.0L*pow(x_1, 6) + (29.0L/14400.0L)*pow(x_1, 5) - 13.0L/3840.0L*pow(x_1, 4) - 691.0L/17280.0L*pow(x_1, 3) + (4027.0L/15360.0L)*pow(x_1, 2) + (353399.0L/230400.0L)*x_1 + 38279.0L/19200.0L;
       return __pp_r525___result;
    }
    static inline O __pp_r526__(const I &x_0, const I &x_1) {
       O __pp_r526___result;
       __pp_r526___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 + (29.0L/4800.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (37.0L/1440.0L)*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 + (289.0L/640.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (2803.0L/2880.0L)*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 3.0L/80.0L*x_0*pow(x_1, 4) - 11.0L/90.0L*x_0*pow(x_1, 3) + (101.0L/480.0L)*x_0*pow(x_1, 2) + (21439.0L/11520.0L)*x_0*x_1 + (209569.0L/76800.0L)*x_0 - 1.0L/800.0L*pow(x_1, 6) - 79.0L/3600.0L*pow(x_1, 5) - 71.0L/480.0L*pow(x_1, 4) - 947.0L/2160.0L*pow(x_1, 3) - 1187.0L/3840.0L*pow(x_1, 2) + (128707.0L/115200.0L)*x_1 + 574799.0L/307200.0L;
       return __pp_r526___result;
    }
    static inline O __pp_r527__(const I &x_0, const I &x_1) {
       O __pp_r527___result;
       __pp_r527___result = (1.0L/4320.0L)*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 + (1.0L/144.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 4)*x_1 + (25.0L/288.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (5.0L/72.0L)*pow(x_0, 3)*pow(x_1, 2) + (25.0L/72.0L)*pow(x_0, 3)*x_1 + (125.0L/216.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (5.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (25.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) + (125.0L/72.0L)*pow(x_0, 2)*x_1 + (625.0L/288.0L)*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (5.0L/144.0L)*x_0*pow(x_1, 4) + (25.0L/72.0L)*x_0*pow(x_1, 3) + (125.0L/72.0L)*x_0*pow(x_1, 2) + (625.0L/144.0L)*x_0*x_1 + (625.0L/144.0L)*x_0 + (1.0L/4320.0L)*pow(x_1, 6) + (1.0L/144.0L)*pow(x_1, 5) + (25.0L/288.0L)*pow(x_1, 4) + (125.0L/216.0L)*pow(x_1, 3) + (625.0L/288.0L)*pow(x_1, 2) + (625.0L/144.0L)*x_1 + 3125.0L/864.0L;
       return __pp_r527___result;
    }
    static inline O __pp_r528__(const I &x_0, const I &x_1) {
       O __pp_r528___result;
       __pp_r528___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 + (29.0L/4800.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) + (37.0L/1440.0L)*pow(x_0, 4)*x_1 + (277.0L/3840.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/30.0L)*pow(x_0, 3)*pow(x_1, 2) + (331.0L/1440.0L)*pow(x_0, 3)*x_1 + (289.0L/640.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (27.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) + (2803.0L/2880.0L)*pow(x_0, 2)*x_1 + (23813.0L/15360.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 1.0L/160.0L*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) + (337.0L/960.0L)*x_0*pow(x_1, 2) + (11327.0L/5760.0L)*x_0*x_1 + (211999.0L/76800.0L)*x_0 - 1.0L/4800.0L*pow(x_1, 6) + (29.0L/14400.0L)*pow(x_1, 5) - 13.0L/3840.0L*pow(x_1, 4) - 691.0L/17280.0L*pow(x_1, 3) + (4027.0L/15360.0L)*pow(x_1, 2) + (353399.0L/230400.0L)*x_1 + 38279.0L/19200.0L;
       return __pp_r528___result;
    }
    static inline O __pp_r529__(const I &x_0, const I &x_1) {
       O __pp_r529___result;
       __pp_r529___result = -7.0L/4800.0L*pow(x_0, 6) + (17.0L/3600.0L)*pow(x_0, 5)*x_1 - 19.0L/14400.0L*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) + (31.0L/1440.0L)*pow(x_0, 4)*x_1 + (67.0L/3840.0L)*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 + (209.0L/17280.0L)*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) + (19.0L/2880.0L)*pow(x_0, 2)*x_1 - 1121.0L/15360.0L*pow(x_0, 2) + (2.0L/225.0L)*x_0*pow(x_1, 5) - 1.0L/720.0L*x_0*pow(x_1, 4) - 1.0L/90.0L*x_0*pow(x_1, 3) - 1.0L/1440.0L*x_0*pow(x_1, 2) + (17.0L/11520.0L)*x_0*x_1 + (401.0L/230400.0L)*x_0 + (1.0L/400.0L)*pow(x_1, 6) + (1.0L/7200.0L)*pow(x_1, 5) + (1.0L/120.0L)*pow(x_1, 4) + (1.0L/8640.0L)*pow(x_1, 3) - 51.0L/640.0L*pow(x_1, 2) + (1.0L/7200.0L)*x_1 + 166747.0L/921600.0L;
       return __pp_r529___result;
    }
    static inline O __pp_r530__(const I &x_0, const I &x_1) {
       O __pp_r530___result;
       __pp_r530___result = -1.0L/675.0L*pow(x_0, 6) + (1.0L/225.0L)*pow(x_0, 5)*x_1 - 1.0L/800.0L*pow(x_0, 5) - 1.0L/180.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/45.0L)*pow(x_0, 4)*x_1 + (5.0L/288.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 2) + (1.0L/36.0L)*pow(x_0, 3)*x_1 + (7.0L/576.0L)*pow(x_0, 3) - 1.0L/720.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/90.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/144.0L)*pow(x_0, 2)*x_1 - 841.0L/11520.0L*pow(x_0, 2) - 7.0L/1800.0L*x_0*pow(x_1, 5) + (1.0L/40.0L)*x_0*pow(x_1, 4) - 5.0L/144.0L*x_0*pow(x_1, 3) + (1.0L/96.0L)*x_0*pow(x_1, 2) - 7.0L/5760.0L*x_0*x_1 + (77.0L/38400.0L)*x_0 - 23.0L/21600.0L*pow(x_1, 6) + (1.0L/3600.0L)*pow(x_1, 5) + (23.0L/1152.0L)*pow(x_1, 4) - 13.0L/864.0L*pow(x_1, 3) - 1643.0L/23040.0L*pow(x_1, 2) - 119.0L/57600.0L*x_1 + 50087.0L/276480.0L;
       return __pp_r530___result;
    }
    static inline O __pp_r531__(const I &x_0, const I &x_1) {
       O __pp_r531___result;
       __pp_r531___result = -7.0L/4800.0L*pow(x_0, 6) + (17.0L/3600.0L)*pow(x_0, 5)*x_1 - 19.0L/14400.0L*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) + (31.0L/1440.0L)*pow(x_0, 4)*x_1 + (67.0L/3840.0L)*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (41.0L/1440.0L)*pow(x_0, 3)*x_1 + (209.0L/17280.0L)*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) + (19.0L/2880.0L)*pow(x_0, 2)*x_1 - 1121.0L/15360.0L*pow(x_0, 2) + (1.0L/1800.0L)*x_0*pow(x_1, 5) + (7.0L/360.0L)*x_0*pow(x_1, 4) - 23.0L/720.0L*x_0*pow(x_1, 3) + (7.0L/720.0L)*x_0*pow(x_1, 2) - 13.0L/11520.0L*x_0*x_1 + (461.0L/230400.0L)*x_0 + (1.0L/2400.0L)*pow(x_1, 6) - 7.0L/3600.0L*pow(x_1, 5) + (41.0L/1920.0L)*pow(x_1, 4) - 67.0L/4320.0L*pow(x_1, 3) - 547.0L/7680.0L*pow(x_1, 2) - 239.0L/115200.0L*x_1 + 166957.0L/921600.0L;
       return __pp_r531___result;
    }
    static inline O __pp_r532__(const I &x_0, const I &x_1) {
       O __pp_r532___result;
       __pp_r532___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 - 7.0L/4800.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (29.0L/1440.0L)*pow(x_0, 4)*x_1 + (199.0L/11520.0L)*pow(x_0, 4) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) + (13.0L/480.0L)*pow(x_0, 3)*x_1 + (23.0L/1920.0L)*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) + (17.0L/2880.0L)*pow(x_0, 2)*x_1 - 673.0L/9216.0L*pow(x_0, 2) - 1.0L/80.0L*x_0*pow(x_1, 4) - 1.0L/60.0L*x_0*pow(x_1, 3) - 1.0L/480.0L*x_0*pow(x_1, 2) + (1.0L/768.0L)*x_0*x_1 + (133.0L/76800.0L)*x_0 - 1.0L/2160.0L*pow(x_1, 6) - 31.0L/7200.0L*pow(x_1, 5) + (1.0L/180.0L)*pow(x_1, 4) - 7.0L/8640.0L*pow(x_1, 3) - 23.0L/288.0L*pow(x_1, 2) + (7.0L/57600.0L)*x_1 + 500239.0L/2764800.0L;
       return __pp_r532___result;
    }
    static inline O __pp_r533__(const I &x_0, const I &x_1) {
       O __pp_r533___result;
       __pp_r533___result = -11.0L/8640.0L*pow(x_0, 6) + (1.0L/180.0L)*pow(x_0, 5)*x_1 - 1.0L/14400.0L*pow(x_0, 5) - 1.0L/288.0L*pow(x_0, 4)*pow(x_1, 2) + (13.0L/480.0L)*pow(x_0, 4)*x_1 + (239.0L/11520.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (59.0L/1440.0L)*pow(x_0, 3)*x_1 + (287.0L/17280.0L)*pow(x_0, 3) - 1.0L/288.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/240.0L)*pow(x_0, 2)*pow(x_1, 2) + (19.0L/960.0L)*pow(x_0, 2)*x_1 - 641.0L/9216.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) - 1.0L/180.0L*x_0*pow(x_1, 4) - 1.0L/360.0L*x_0*pow(x_1, 3) + (17.0L/1440.0L)*x_0*pow(x_1, 2) + (19.0L/2304.0L)*x_0*x_1 + (719.0L/230400.0L)*x_0 - 1.0L/4320.0L*pow(x_1, 6) - 7.0L/2400.0L*pow(x_1, 5) + (13.0L/1440.0L)*pow(x_1, 4) + (11.0L/2880.0L)*pow(x_1, 3) - 11.0L/144.0L*pow(x_1, 2) + (29.0L/19200.0L)*x_1 + 500879.0L/2764800.0L;
       return __pp_r533___result;
    }
    static inline O __pp_r534__(const I &x_0, const I &x_1) {
       O __pp_r534___result;
       __pp_r534___result = -7.0L/5400.0L*pow(x_0, 6) + (19.0L/3600.0L)*pow(x_0, 5)*x_1 - 1.0L/3600.0L*pow(x_0, 5) - 7.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (1.0L/40.0L)*pow(x_0, 4)*x_1 + (23.0L/1152.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (5.0L/144.0L)*pow(x_0, 3)*x_1 + (13.0L/864.0L)*pow(x_0, 3) - 13.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) + (1.0L/96.0L)*pow(x_0, 2)*x_1 - 1643.0L/23040.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/45.0L*x_0*pow(x_1, 4) - 1.0L/36.0L*x_0*pow(x_1, 3) - 1.0L/144.0L*x_0*pow(x_1, 2) + (7.0L/5760.0L)*x_0*x_1 + (119.0L/57600.0L)*x_0 - 37.0L/21600.0L*pow(x_1, 6) - 23.0L/2400.0L*pow(x_1, 5) - 1.0L/288.0L*pow(x_1, 4) - 5.0L/576.0L*pow(x_1, 3) - 961.0L/11520.0L*pow(x_1, 2) - 23.0L/38400.0L*x_1 + 10003.0L/55296.0L;
       return __pp_r534___result;
    }
    static inline O __pp_r535__(const I &x_0, const I &x_1) {
       O __pp_r535___result;
       __pp_r535___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (47.0L/4800.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/160.0L*pow(x_0, 4)*x_1 + (601.0L/11520.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) - 59.0L/1440.0L*pow(x_0, 3)*x_1 + (403.0L/5760.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (31.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 77.0L/960.0L*pow(x_0, 2)*x_1 - 863.0L/46080.0L*pow(x_0, 2) + (31.0L/3600.0L)*x_0*pow(x_1, 5) + (1.0L/480.0L)*x_0*pow(x_1, 4) - 41.0L/1440.0L*x_0*pow(x_1, 3) + (41.0L/960.0L)*x_0*pow(x_1, 2) - 19.0L/360.0L*x_0*x_1 + (739.0L/25600.0L)*x_0 + (109.0L/43200.0L)*pow(x_1, 6) - 1.0L/4800.0L*pow(x_1, 5) + (121.0L/11520.0L)*pow(x_1, 4) - 41.0L/5760.0L*pow(x_1, 3) - 3047.0L/46080.0L*pow(x_1, 2) - 1031.0L/76800.0L*x_1 + 257933.0L/1382400.0L;
       return __pp_r535___result;
    }
    static inline O __pp_r536__(const I &x_0, const I &x_1) {
       O __pp_r536___result;
       __pp_r536___result = (71.0L/7200.0L)*pow(x_0, 5) - 1.0L/180.0L*pow(x_0, 4)*x_1 + (5.0L/96.0L)*pow(x_0, 4) + (7.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 1.0L/24.0L*pow(x_0, 3)*x_1 + (121.0L/1728.0L)*pow(x_0, 3) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 23.0L/288.0L*pow(x_0, 2)*x_1 - 3.0L/160.0L*pow(x_0, 2) - 1.0L/240.0L*x_0*pow(x_1, 5) + (41.0L/1440.0L)*x_0*pow(x_1, 4) - 5.0L/96.0L*x_0*pow(x_1, 3) + (31.0L/576.0L)*x_0*pow(x_1, 2) - 71.0L/1280.0L*x_0*x_1 + (839.0L/28800.0L)*x_0 - 1.0L/960.0L*pow(x_1, 6) - 1.0L/14400.0L*pow(x_1, 5) + (17.0L/768.0L)*pow(x_1, 4) - 77.0L/3456.0L*pow(x_1, 3) - 887.0L/15360.0L*pow(x_1, 2) - 3601.0L/230400.0L*x_1 + 34433.0L/184320.0L;
       return __pp_r536___result;
    }
    static inline O __pp_r537__(const I &x_0, const I &x_1) {
       O __pp_r537___result;
       __pp_r537___result = (1.0L/43200.0L)*pow(x_0, 6) + (1.0L/3600.0L)*pow(x_0, 5)*x_1 + (47.0L/4800.0L)*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/160.0L*pow(x_0, 4)*x_1 + (601.0L/11520.0L)*pow(x_0, 4) + (1.0L/270.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) - 59.0L/1440.0L*pow(x_0, 3)*x_1 + (403.0L/5760.0L)*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) + (31.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 77.0L/960.0L*pow(x_0, 2)*x_1 - 863.0L/46080.0L*pow(x_0, 2) + (1.0L/3600.0L)*x_0*pow(x_1, 5) + (11.0L/480.0L)*x_0*pow(x_1, 4) - 71.0L/1440.0L*x_0*pow(x_1, 3) + (17.0L/320.0L)*x_0*pow(x_1, 2) - 319.0L/5760.0L*x_0*x_1 + (2237.0L/76800.0L)*x_0 + (19.0L/43200.0L)*pow(x_1, 6) - 11.0L/4800.0L*pow(x_1, 5) + (271.0L/11520.0L)*pow(x_1, 4) - 131.0L/5760.0L*pow(x_1, 3) - 2657.0L/46080.0L*pow(x_1, 2) - 1201.0L/76800.0L*x_1 + 32281.0L/172800.0L;
       return __pp_r537___result;
    }
    static inline O __pp_r538__(const I &x_0, const I &x_1) {
       O __pp_r538___result;
       __pp_r538___result = -1.0L/43200.0L*pow(x_0, 6) - 1.0L/3600.0L*pow(x_0, 5)*x_1 + (139.0L/14400.0L)*pow(x_0, 5) - 1.0L/720.0L*pow(x_0, 4)*pow(x_1, 2) - 11.0L/1440.0L*pow(x_0, 4)*x_1 + (599.0L/11520.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) - 61.0L/1440.0L*pow(x_0, 3)*x_1 + (1207.0L/17280.0L)*pow(x_0, 3) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 4) - 7.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (29.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 233.0L/2880.0L*pow(x_0, 2)*x_1 - 173.0L/9216.0L*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) - 13.0L/1440.0L*x_0*pow(x_1, 4) - 49.0L/1440.0L*x_0*pow(x_1, 3) + (119.0L/2880.0L)*x_0*pow(x_1, 2) - 61.0L/1152.0L*x_0*x_1 + (6649.0L/230400.0L)*x_0 - 19.0L/43200.0L*pow(x_1, 6) - 67.0L/14400.0L*pow(x_1, 5) + (89.0L/11520.0L)*pow(x_1, 4) - 139.0L/17280.0L*pow(x_1, 3) - 611.0L/9216.0L*pow(x_1, 2) - 3097.0L/230400.0L*x_1 + 64483.0L/345600.0L;
       return __pp_r538___result;
    }
    static inline O __pp_r539__(const I &x_0, const I &x_1) {
       O __pp_r539___result;
       __pp_r539___result = (1.0L/4800.0L)*pow(x_0, 6) + (1.0L/900.0L)*pow(x_0, 5)*x_1 + (53.0L/4800.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/1440.0L*pow(x_0, 4)*x_1 + (71.0L/1280.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/40.0L)*pow(x_0, 3)*pow(x_1, 2) - 41.0L/1440.0L*pow(x_0, 3)*x_1 + (143.0L/1920.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 193.0L/2880.0L*pow(x_0, 2)*x_1 - 47.0L/3072.0L*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 1.0L/480.0L*x_0*pow(x_1, 4) - 29.0L/1440.0L*x_0*pow(x_1, 3) + (53.0L/960.0L)*x_0*pow(x_1, 2) - 53.0L/1152.0L*x_0*x_1 + (2323.0L/76800.0L)*x_0 - 1.0L/4800.0L*pow(x_1, 6) - 47.0L/14400.0L*pow(x_1, 5) + (43.0L/3840.0L)*pow(x_1, 4) - 59.0L/17280.0L*pow(x_1, 3) - 193.0L/3072.0L*pow(x_1, 2) - 2777.0L/230400.0L*x_1 + 21521.0L/115200.0L;
       return __pp_r539___result;
    }
    static inline O __pp_r540__(const I &x_0, const I &x_1) {
       O __pp_r540___result;
       __pp_r540___result = (1.0L/5400.0L)*pow(x_0, 6) + (1.0L/1200.0L)*pow(x_0, 5)*x_1 + (13.0L/1200.0L)*pow(x_0, 5) + (1.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/360.0L*pow(x_0, 4)*x_1 + (7.0L/128.0L)*pow(x_0, 4) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) - 5.0L/144.0L*pow(x_0, 3)*x_1 + (7.0L/96.0L)*pow(x_0, 3) - 11.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/45.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 11.0L/144.0L*pow(x_0, 2)*x_1 - 131.0L/7680.0L*pow(x_0, 2) - 1.0L/300.0L*x_0*pow(x_1, 5) - 3.0L/160.0L*x_0*pow(x_1, 4) - 13.0L/288.0L*x_0*pow(x_1, 3) + (7.0L/192.0L)*x_0*pow(x_1, 2) - 611.0L/11520.0L*x_0*x_1 + (1121.0L/38400.0L)*x_0 - 73.0L/43200.0L*pow(x_1, 6) - 143.0L/14400.0L*pow(x_1, 5) - 1.0L/768.0L*pow(x_1, 4) - 55.0L/3456.0L*pow(x_1, 3) - 1073.0L/15360.0L*pow(x_1, 2) - 3263.0L/230400.0L*x_1 + 6877.0L/36864.0L;
       return __pp_r540___result;
    }
    static inline O __pp_r541__(const I &x_0, const I &x_1) {
       O __pp_r541___result;
       __pp_r541___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 + (1.0L/100.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 23.0L/720.0L*pow(x_0, 4)*x_1 + (1.0L/80.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/40.0L)*pow(x_0, 3)*pow(x_1, 2) - 11.0L/90.0L*pow(x_0, 3)*x_1 - 43.0L/480.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 299.0L/1440.0L*pow(x_0, 2)*x_1 - 221.0L/768.0L*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 1.0L/480.0L*x_0*pow(x_1, 4) - 29.0L/1440.0L*x_0*pow(x_1, 3) + (53.0L/960.0L)*x_0*pow(x_1, 2) - 349.0L/2304.0L*x_0*x_1 - 7141.0L/38400.0L*x_0 - 1.0L/4800.0L*pow(x_1, 6) - 47.0L/14400.0L*pow(x_1, 5) + (43.0L/3840.0L)*pow(x_1, 4) - 59.0L/17280.0L*pow(x_1, 3) - 193.0L/3072.0L*pow(x_1, 2) - 10067.0L/230400.0L*x_1 + 110203.0L/921600.0L;
       return __pp_r541___result;
    }
    static inline O __pp_r542__(const I &x_0, const I &x_1) {
       O __pp_r542___result;
       __pp_r542___result = (53.0L/43200.0L)*pow(x_0, 6) - 1.0L/300.0L*pow(x_0, 5)*x_1 + (47.0L/4800.0L)*pow(x_0, 5) + (1.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 49.0L/1440.0L*pow(x_0, 4)*x_1 + (3.0L/256.0L)*pow(x_0, 4) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) - 37.0L/288.0L*pow(x_0, 3)*x_1 - 35.0L/384.0L*pow(x_0, 3) - 11.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/45.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 125.0L/576.0L*pow(x_0, 2)*x_1 - 4447.0L/15360.0L*pow(x_0, 2) - 1.0L/300.0L*x_0*pow(x_1, 5) - 3.0L/160.0L*x_0*pow(x_1, 4) - 13.0L/288.0L*x_0*pow(x_1, 3) + (7.0L/192.0L)*x_0*pow(x_1, 2) - 913.0L/5760.0L*x_0*x_1 - 14363.0L/76800.0L*x_0 - 73.0L/43200.0L*pow(x_1, 6) - 143.0L/14400.0L*pow(x_1, 5) - 1.0L/768.0L*pow(x_1, 4) - 55.0L/3456.0L*pow(x_1, 3) - 1073.0L/15360.0L*pow(x_1, 2) - 10553.0L/230400.0L*x_1 + 2749.0L/23040.0L;
       return __pp_r542___result;
    }
    static inline O __pp_r543__(const I &x_0, const I &x_1) {
       O __pp_r543___result;
       __pp_r543___result = (7.0L/4800.0L)*pow(x_0, 6) - 7.0L/3600.0L*pow(x_0, 5)*x_1 + (181.0L/14400.0L)*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) - 29.0L/1440.0L*pow(x_0, 4)*x_1 + (59.0L/2304.0L)*pow(x_0, 4) + (1.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) + (2.0L/45.0L)*pow(x_0, 3)*pow(x_1, 2) - 7.0L/96.0L*pow(x_0, 3)*x_1 - 187.0L/3456.0L*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 61.0L/576.0L*pow(x_0, 2)*x_1 - 10781.0L/46080.0L*pow(x_0, 2) - 7.0L/3600.0L*x_0*pow(x_1, 5) - 7.0L/1440.0L*x_0*pow(x_1, 4) + (1.0L/96.0L)*x_0*pow(x_1, 3) + (85.0L/576.0L)*x_0*pow(x_1, 2) - 91.0L/1920.0L*x_0*x_1 - 32849.0L/230400.0L*x_0 - 7.0L/4800.0L*pow(x_1, 6) - 103.0L/14400.0L*pow(x_1, 5) + (29.0L/2304.0L)*pow(x_1, 4) + (73.0L/3456.0L)*pow(x_1, 3) - 659.0L/46080.0L*pow(x_1, 2) - 313.0L/230400.0L*x_1 + 9271.0L/69120.0L;
       return __pp_r543___result;
    }
    static inline O __pp_r544__(const I &x_0, const I &x_1) {
       O __pp_r544___result;
       __pp_r544___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (31.0L/2400.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/60.0L*pow(x_0, 4)*x_1 + (1.0L/36.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 3)*x_1 - 3.0L/64.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/16.0L*pow(x_0, 2)*x_1 - 2539.0L/11520.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) + (11.0L/480.0L)*x_0*pow(x_1, 4) + (23.0L/288.0L)*x_0*pow(x_1, 3) + (15.0L/64.0L)*x_0*pow(x_1, 2) + (79.0L/11520.0L)*x_0*x_1 - 2477.0L/19200.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (19.0L/4800.0L)*pow(x_1, 5) + (109.0L/2304.0L)*pow(x_1, 4) + (91.0L/1152.0L)*pow(x_1, 3) + (1841.0L/46080.0L)*pow(x_1, 2) + (1979.0L/76800.0L)*x_1 + 77293.0L/552960.0L;
       return __pp_r544___result;
    }
    static inline O __pp_r545__(const I &x_0, const I &x_1) {
       O __pp_r545___result;
       __pp_r545___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 + (1.0L/100.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 23.0L/720.0L*pow(x_0, 4)*x_1 + (1.0L/80.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/40.0L)*pow(x_0, 3)*pow(x_1, 2) - 11.0L/90.0L*pow(x_0, 3)*x_1 - 43.0L/480.0L*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (13.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 299.0L/1440.0L*pow(x_0, 2)*x_1 - 221.0L/768.0L*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 1.0L/480.0L*x_0*pow(x_1, 4) - 29.0L/1440.0L*x_0*pow(x_1, 3) + (53.0L/960.0L)*x_0*pow(x_1, 2) - 349.0L/2304.0L*x_0*x_1 - 7141.0L/38400.0L*x_0 - 1.0L/4800.0L*pow(x_1, 6) - 47.0L/14400.0L*pow(x_1, 5) + (43.0L/3840.0L)*pow(x_1, 4) - 59.0L/17280.0L*pow(x_1, 3) - 193.0L/3072.0L*pow(x_1, 2) - 10067.0L/230400.0L*x_1 + 110203.0L/921600.0L;
       return __pp_r545___result;
    }
    static inline O __pp_r546__(const I &x_0, const I &x_1) {
       O __pp_r546___result;
       __pp_r546___result = (53.0L/43200.0L)*pow(x_0, 6) - 1.0L/300.0L*pow(x_0, 5)*x_1 + (47.0L/4800.0L)*pow(x_0, 5) + (1.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 49.0L/1440.0L*pow(x_0, 4)*x_1 + (3.0L/256.0L)*pow(x_0, 4) - 1.0L/360.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/60.0L)*pow(x_0, 3)*pow(x_1, 2) - 37.0L/288.0L*pow(x_0, 3)*x_1 - 35.0L/384.0L*pow(x_0, 3) - 11.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/45.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 125.0L/576.0L*pow(x_0, 2)*x_1 - 4447.0L/15360.0L*pow(x_0, 2) - 1.0L/300.0L*x_0*pow(x_1, 5) - 3.0L/160.0L*x_0*pow(x_1, 4) - 13.0L/288.0L*x_0*pow(x_1, 3) + (7.0L/192.0L)*x_0*pow(x_1, 2) - 913.0L/5760.0L*x_0*x_1 - 14363.0L/76800.0L*x_0 - 73.0L/43200.0L*pow(x_1, 6) - 143.0L/14400.0L*pow(x_1, 5) - 1.0L/768.0L*pow(x_1, 4) - 55.0L/3456.0L*pow(x_1, 3) - 1073.0L/15360.0L*pow(x_1, 2) - 10553.0L/230400.0L*x_1 + 2749.0L/23040.0L;
       return __pp_r546___result;
    }
    static inline O __pp_r547__(const I &x_0, const I &x_1) {
       O __pp_r547___result;
       __pp_r547___result = (7.0L/4800.0L)*pow(x_0, 6) - 7.0L/3600.0L*pow(x_0, 5)*x_1 + (181.0L/14400.0L)*pow(x_0, 5) + (1.0L/240.0L)*pow(x_0, 4)*pow(x_1, 2) - 29.0L/1440.0L*pow(x_0, 4)*x_1 + (59.0L/2304.0L)*pow(x_0, 4) + (1.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) + (2.0L/45.0L)*pow(x_0, 3)*pow(x_1, 2) - 7.0L/96.0L*pow(x_0, 3)*x_1 - 187.0L/3456.0L*pow(x_0, 3) - 1.0L/240.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 61.0L/576.0L*pow(x_0, 2)*x_1 - 10781.0L/46080.0L*pow(x_0, 2) - 7.0L/3600.0L*x_0*pow(x_1, 5) - 7.0L/1440.0L*x_0*pow(x_1, 4) + (1.0L/96.0L)*x_0*pow(x_1, 3) + (85.0L/576.0L)*x_0*pow(x_1, 2) - 91.0L/1920.0L*x_0*x_1 - 32849.0L/230400.0L*x_0 - 7.0L/4800.0L*pow(x_1, 6) - 103.0L/14400.0L*pow(x_1, 5) + (29.0L/2304.0L)*pow(x_1, 4) + (73.0L/3456.0L)*pow(x_1, 3) - 659.0L/46080.0L*pow(x_1, 2) - 313.0L/230400.0L*x_1 + 9271.0L/69120.0L;
       return __pp_r547___result;
    }
    static inline O __pp_r548__(const I &x_0, const I &x_1) {
       O __pp_r548___result;
       __pp_r548___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/600.0L*pow(x_0, 5)*x_1 + (31.0L/2400.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/60.0L*pow(x_0, 4)*x_1 + (1.0L/36.0L)*pow(x_0, 4) + (1.0L/180.0L)*pow(x_0, 3)*pow(x_1, 3) + (7.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 3)*x_1 - 3.0L/64.0L*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/30.0L)*pow(x_0, 2)*pow(x_1, 3) + (19.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 1.0L/16.0L*pow(x_0, 2)*x_1 - 2539.0L/11520.0L*pow(x_0, 2) + (1.0L/400.0L)*x_0*pow(x_1, 5) + (11.0L/480.0L)*x_0*pow(x_1, 4) + (23.0L/288.0L)*x_0*pow(x_1, 3) + (15.0L/64.0L)*x_0*pow(x_1, 2) + (79.0L/11520.0L)*x_0*x_1 - 2477.0L/19200.0L*x_0 + (1.0L/43200.0L)*pow(x_1, 6) + (19.0L/4800.0L)*pow(x_1, 5) + (109.0L/2304.0L)*pow(x_1, 4) + (91.0L/1152.0L)*pow(x_1, 3) + (1841.0L/46080.0L)*pow(x_1, 2) + (1979.0L/76800.0L)*x_1 + 77293.0L/552960.0L;
       return __pp_r548___result;
    }
    static inline O __pp_r549__(const I &x_0, const I &x_1) {
       O __pp_r549___result;
       __pp_r549___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (91.0L/7200.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 7.0L/360.0L*pow(x_0, 4)*x_1 + (19.0L/288.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (17.0L/360.0L)*pow(x_0, 3)*pow(x_1, 2) - 7.0L/72.0L*pow(x_0, 3)*x_1 + (185.0L/1728.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 11.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 55.0L/288.0L*pow(x_0, 2)*x_1 + (53.0L/1440.0L)*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) + (61.0L/1440.0L)*x_0*pow(x_1, 4) - 31.0L/288.0L*x_0*pow(x_1, 3) + (95.0L/576.0L)*x_0*pow(x_1, 2) - 1919.0L/11520.0L*x_0*x_1 + (2119.0L/28800.0L)*x_0 - 7.0L/8640.0L*pow(x_1, 6) - 41.0L/14400.0L*pow(x_1, 5) + (83.0L/2304.0L)*pow(x_1, 4) - 205.0L/3456.0L*pow(x_1, 3) - 101.0L/46080.0L*pow(x_1, 2) - 13841.0L/230400.0L*x_1 + 111491.0L/552960.0L;
       return __pp_r549___result;
    }
    static inline O __pp_r550__(const I &x_0, const I &x_1) {
       O __pp_r550___result;
       __pp_r550___result = (23.0L/7200.0L)*pow(x_0, 6) - 37.0L/3600.0L*pow(x_0, 5)*x_1 + (7.0L/160.0L)*pow(x_0, 5) + (7.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 7.0L/72.0L*pow(x_0, 4)*x_1 + (97.0L/480.0L)*pow(x_0, 4) - 13.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/8.0L)*pow(x_0, 3)*pow(x_1, 2) - 133.0L/360.0L*pow(x_0, 3)*x_1 + (1223.0L/2880.0L)*pow(x_0, 3) + (1.0L/160.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) - 961.0L/1440.0L*pow(x_0, 2)*x_1 + (871.0L/1920.0L)*pow(x_0, 2) - 11.0L/1800.0L*x_0*pow(x_1, 5) + (5.0L/96.0L)*x_0*pow(x_1, 4) - 253.0L/1440.0L*x_0*pow(x_1, 3) + (129.0L/320.0L)*x_0*pow(x_1, 2) - 6721.0L/11520.0L*x_0*x_1 + (1403.0L/3840.0L)*x_0 - 11.0L/14400.0L*pow(x_1, 6) - 11.0L/2880.0L*pow(x_1, 5) + (57.0L/1280.0L)*pow(x_1, 4) - 1711.0L/17280.0L*pow(x_1, 3) + (1567.0L/15360.0L)*pow(x_1, 2) - 9491.0L/46080.0L*x_1 + 264251.0L/921600.0L;
       return __pp_r550___result;
    }
    static inline O __pp_r551__(const I &x_0, const I &x_1) {
       O __pp_r551___result;
       __pp_r551___result = (31.0L/14400.0L)*pow(x_0, 6) - 11.0L/1800.0L*pow(x_0, 5)*x_1 + (19.0L/960.0L)*pow(x_0, 5) + (7.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 19.0L/288.0L*pow(x_0, 4)*x_1 + (221.0L/3840.0L)*pow(x_0, 4) - 13.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/8.0L)*pow(x_0, 3)*pow(x_1, 2) - 397.0L/1440.0L*pow(x_0, 3)*x_1 + (151.0L/5760.0L)*pow(x_0, 3) + (1.0L/160.0L)*pow(x_0, 2)*pow(x_1, 4) - 5.0L/72.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/20.0L)*pow(x_0, 2)*pow(x_1, 2) - 1517.0L/2880.0L*pow(x_0, 2)*x_1 - 1807.0L/15360.0L*pow(x_0, 2) - 11.0L/1800.0L*x_0*pow(x_1, 5) + (5.0L/96.0L)*x_0*pow(x_1, 4) - 253.0L/1440.0L*x_0*pow(x_1, 3) + (129.0L/320.0L)*x_0*pow(x_1, 2) - 2753.0L/5760.0L*x_0*x_1 - 787.0L/15360.0L*x_0 - 11.0L/14400.0L*pow(x_1, 6) - 11.0L/2880.0L*pow(x_1, 5) + (57.0L/1280.0L)*pow(x_1, 4) - 1711.0L/17280.0L*pow(x_1, 3) + (1567.0L/15360.0L)*pow(x_1, 2) - 8033.0L/46080.0L*x_1 + 18907.0L/115200.0L;
       return __pp_r551___result;
    }
    static inline O __pp_r552__(const I &x_0, const I &x_1) {
       O __pp_r552___result;
       __pp_r552___result = (11.0L/43200.0L)*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (181.0L/14400.0L)*pow(x_0, 5) + (7.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 29.0L/1440.0L*pow(x_0, 4)*x_1 + (761.0L/11520.0L)*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (2.0L/45.0L)*pow(x_0, 3)*pow(x_1, 2) - 139.0L/1440.0L*pow(x_0, 3)*x_1 + (1849.0L/17280.0L)*pow(x_0, 3) + (13.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 13.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (71.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 551.0L/2880.0L*pow(x_0, 2)*x_1 + (1697.0L/46080.0L)*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) + (53.0L/1440.0L)*x_0*pow(x_1, 4) - 151.0L/1440.0L*x_0*pow(x_1, 3) + (473.0L/2880.0L)*x_0*pow(x_1, 2) - 959.0L/5760.0L*x_0*x_1 + (16951.0L/230400.0L)*x_0 + (29.0L/43200.0L)*pow(x_1, 6) - 73.0L/14400.0L*pow(x_1, 5) + (431.0L/11520.0L)*pow(x_1, 4) - 1033.0L/17280.0L*pow(x_1, 3) - 97.0L/46080.0L*pow(x_1, 2) - 13843.0L/230400.0L*x_1 + 34841.0L/172800.0L;
       return __pp_r552___result;
    }
    static inline O __pp_r553__(const I &x_0, const I &x_1) {
       O __pp_r553___result;
       __pp_r553___result = (139.0L/43200.0L)*pow(x_0, 6) - 1.0L/100.0L*pow(x_0, 5)*x_1 + (629.0L/14400.0L)*pow(x_0, 5) + (23.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 47.0L/480.0L*pow(x_0, 4)*x_1 + (2329.0L/11520.0L)*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (11.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) - 59.0L/160.0L*pow(x_0, 3)*x_1 + (7337.0L/17280.0L)*pow(x_0, 3) + (17.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 3.0L/40.0L*pow(x_0, 2)*pow(x_1, 3) + (169.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 641.0L/960.0L*pow(x_0, 2)*x_1 + (4181.0L/9216.0L)*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) + (67.0L/1440.0L)*x_0*pow(x_1, 4) - 83.0L/480.0L*x_0*pow(x_1, 3) + (1159.0L/2880.0L)*x_0*pow(x_1, 2) - 7.0L/12.0L*x_0*x_1 + (84179.0L/230400.0L)*x_0 + (31.0L/43200.0L)*pow(x_1, 6) - 29.0L/4800.0L*pow(x_1, 5) + (529.0L/11520.0L)*pow(x_1, 4) - 191.0L/1920.0L*pow(x_1, 3) + (941.0L/9216.0L)*pow(x_1, 2) - 5273.0L/25600.0L*x_1 + 396377.0L/1382400.0L;
       return __pp_r553___result;
    }
    static inline O __pp_r554__(const I &x_0, const I &x_1) {
       O __pp_r554___result;
       __pp_r554___result = (47.0L/21600.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 + (71.0L/3600.0L)*pow(x_0, 5) + (23.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/15.0L*pow(x_0, 4)*x_1 + (83.0L/1440.0L)*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (11.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) - 11.0L/40.0L*pow(x_0, 3)*x_1 + (113.0L/4320.0L)*pow(x_0, 3) + (17.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 3.0L/40.0L*pow(x_0, 2)*pow(x_1, 3) + (169.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 253.0L/480.0L*pow(x_0, 2)*x_1 - 271.0L/2304.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) + (67.0L/1440.0L)*x_0*pow(x_1, 4) - 83.0L/480.0L*x_0*pow(x_1, 3) + (1159.0L/2880.0L)*x_0*pow(x_1, 2) - 367.0L/768.0L*x_0*x_1 - 5903.0L/115200.0L*x_0 + (31.0L/43200.0L)*pow(x_1, 6) - 29.0L/4800.0L*pow(x_1, 5) + (529.0L/11520.0L)*pow(x_1, 4) - 191.0L/1920.0L*pow(x_1, 3) + (941.0L/9216.0L)*pow(x_1, 2) - 4463.0L/25600.0L*x_1 + 453769.0L/2764800.0L;
       return __pp_r554___result;
    }
    static inline O __pp_r555__(const I &x_0, const I &x_1) {
       O __pp_r555___result;
       __pp_r555___result = -7.0L/8640.0L*pow(x_0, 6) + (1.0L/360.0L)*pow(x_0, 5)*x_1 - 97.0L/4800.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (49.0L/1440.0L)*pow(x_0, 4)*x_1 - 643.0L/3840.0L*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/40.0L)*pow(x_0, 3)*pow(x_1, 2) + (251.0L/1440.0L)*pow(x_0, 3)*x_1 - 3737.0L/5760.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 7.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (1.0L/80.0L)*pow(x_0, 2)*pow(x_1, 2) + (1399.0L/2880.0L)*pow(x_0, 2)*x_1 - 19303.0L/15360.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) + (19.0L/480.0L)*x_0*pow(x_1, 4) - 91.0L/1440.0L*x_0*pow(x_1, 3) - 33.0L/320.0L*x_0*pow(x_1, 2) + (119.0L/180.0L)*x_0*x_1 - 82667.0L/76800.0L*x_0 - 7.0L/8640.0L*pow(x_1, 6) - 37.0L/14400.0L*pow(x_1, 5) + (39.0L/1280.0L)*pow(x_1, 4) - 253.0L/17280.0L*pow(x_1, 3) - 2807.0L/15360.0L*pow(x_1, 2) + (77933.0L/230400.0L)*x_1 - 101519.0L/460800.0L;
       return __pp_r555___result;
    }
    static inline O __pp_r556__(const I &x_0, const I &x_1) {
       O __pp_r556___result;
       __pp_r556___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 - 157.0L/4800.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (139.0L/1440.0L)*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/10.0L*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 - 5897.0L/5760.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (19.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (4639.0L/2880.0L)*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) - 11.0L/480.0L*x_0*pow(x_1, 4) + (449.0L/1440.0L)*x_0*pow(x_1, 3) - 393.0L/320.0L*x_0*pow(x_1, 2) + (1691.0L/720.0L)*x_0*x_1 - 160427.0L/76800.0L*x_0 - 13.0L/8640.0L*pow(x_1, 6) + (143.0L/14400.0L)*pow(x_1, 5) - 81.0L/1280.0L*pow(x_1, 4) + (6227.0L/17280.0L)*pow(x_1, 3) - 15767.0L/15360.0L*pow(x_1, 2) + (311213.0L/230400.0L)*x_1 - 334799.0L/460800.0L;
       return __pp_r556___result;
    }
    static inline O __pp_r557__(const I &x_0, const I &x_1) {
       O __pp_r557___result;
       __pp_r557___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 - 157.0L/4800.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (139.0L/1440.0L)*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/10.0L*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 - 5897.0L/5760.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (19.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (4639.0L/2880.0L)*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) - 1.0L/720.0L*x_0*pow(x_1, 5) - 11.0L/480.0L*x_0*pow(x_1, 4) + (449.0L/1440.0L)*x_0*pow(x_1, 3) - 393.0L/320.0L*x_0*pow(x_1, 2) + (1691.0L/720.0L)*x_0*x_1 - 160427.0L/76800.0L*x_0 - 13.0L/8640.0L*pow(x_1, 6) + (143.0L/14400.0L)*pow(x_1, 5) - 81.0L/1280.0L*pow(x_1, 4) + (6227.0L/17280.0L)*pow(x_1, 3) - 15767.0L/15360.0L*pow(x_1, 2) + (311213.0L/230400.0L)*x_1 - 334799.0L/460800.0L;
       return __pp_r557___result;
    }
    static inline O __pp_r558__(const I &x_0, const I &x_1) {
       O __pp_r558___result;
       __pp_r558___result = -17.0L/21600.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 73.0L/3600.0L*pow(x_0, 5) + (7.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/30.0L)*pow(x_0, 4)*x_1 - 241.0L/1440.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/45.0L)*pow(x_0, 3)*pow(x_1, 2) + (7.0L/40.0L)*pow(x_0, 3)*x_1 - 2803.0L/4320.0L*pow(x_0, 3) + (13.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/40.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (233.0L/480.0L)*pow(x_0, 2)*x_1 - 14477.0L/11520.0L*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) + (49.0L/1440.0L)*x_0*pow(x_1, 4) - 29.0L/480.0L*x_0*pow(x_1, 3) - 299.0L/2880.0L*x_0*pow(x_1, 2) + (2539.0L/3840.0L)*x_0*x_1 - 124001.0L/115200.0L*x_0 + (29.0L/43200.0L)*pow(x_1, 6) - 23.0L/4800.0L*pow(x_1, 5) + (367.0L/11520.0L)*pow(x_1, 4) - 29.0L/1920.0L*pow(x_1, 3) - 8417.0L/46080.0L*pow(x_1, 2) + (8659.0L/25600.0L)*x_1 - 609113.0L/2764800.0L;
       return __pp_r558___result;
    }
    static inline O __pp_r559__(const I &x_0, const I &x_1) {
       O __pp_r559___result;
       __pp_r559___result = -1.0L/675.0L*pow(x_0, 6) + (13.0L/1800.0L)*pow(x_0, 5)*x_1 - 59.0L/1800.0L*pow(x_0, 5) - 1.0L/180.0L*pow(x_0, 4)*pow(x_1, 2) + (23.0L/240.0L)*pow(x_0, 4)*x_1 - 47.0L/180.0L*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) - 37.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (11.0L/20.0L)*pow(x_0, 3)*x_1 - 4423.0L/4320.0L*pow(x_0, 3) - 1.0L/720.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/10.0L)*pow(x_0, 2)*pow(x_1, 3) - 263.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (773.0L/480.0L)*pow(x_0, 2)*x_1 - 24197.0L/11520.0L*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) - 41.0L/1440.0L*x_0*pow(x_1, 4) + (151.0L/480.0L)*x_0*pow(x_1, 3) - 3539.0L/2880.0L*x_0*pow(x_1, 2) + (9019.0L/3840.0L)*x_0*x_1 - 240641.0L/115200.0L*x_0 - 1.0L/43200.0L*pow(x_1, 6) + (37.0L/4800.0L)*pow(x_1, 5) - 713.0L/11520.0L*pow(x_1, 4) + (691.0L/1920.0L)*pow(x_1, 3) - 47297.0L/46080.0L*pow(x_1, 2) + (34579.0L/25600.0L)*x_1 - 2008793.0L/2764800.0L;
       return __pp_r559___result;
    }
    static inline O __pp_r560__(const I &x_0, const I &x_1) {
       O __pp_r560___result;
       __pp_r560___result = -1.0L/675.0L*pow(x_0, 6) + (13.0L/1800.0L)*pow(x_0, 5)*x_1 - 59.0L/1800.0L*pow(x_0, 5) - 1.0L/180.0L*pow(x_0, 4)*pow(x_1, 2) + (23.0L/240.0L)*pow(x_0, 4)*x_1 - 47.0L/180.0L*pow(x_0, 4) + (7.0L/540.0L)*pow(x_0, 3)*pow(x_1, 3) - 37.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (11.0L/20.0L)*pow(x_0, 3)*x_1 - 4423.0L/4320.0L*pow(x_0, 3) - 1.0L/720.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/10.0L)*pow(x_0, 2)*pow(x_1, 3) - 263.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) + (773.0L/480.0L)*pow(x_0, 2)*x_1 - 24197.0L/11520.0L*pow(x_0, 2) + (11.0L/3600.0L)*x_0*pow(x_1, 5) - 41.0L/1440.0L*x_0*pow(x_1, 4) + (151.0L/480.0L)*x_0*pow(x_1, 3) - 3539.0L/2880.0L*x_0*pow(x_1, 2) + (9019.0L/3840.0L)*x_0*x_1 - 240641.0L/115200.0L*x_0 - 1.0L/43200.0L*pow(x_1, 6) + (37.0L/4800.0L)*pow(x_1, 5) - 713.0L/11520.0L*pow(x_1, 4) + (691.0L/1920.0L)*pow(x_1, 3) - 47297.0L/46080.0L*pow(x_1, 2) + (34579.0L/25600.0L)*x_1 - 2008793.0L/2764800.0L;
       return __pp_r560___result;
    }
    static inline O __pp_r561__(const I &x_0, const I &x_1) {
       O __pp_r561___result;
       __pp_r561___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 - 157.0L/4800.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (139.0L/1440.0L)*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/10.0L*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 - 5897.0L/5760.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (19.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (4639.0L/2880.0L)*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) - 13.0L/240.0L*x_0*pow(x_1, 4) + (73.0L/180.0L)*x_0*pow(x_1, 3) - 219.0L/160.0L*x_0*pow(x_1, 2) + (28271.0L/11520.0L)*x_0*x_1 - 162857.0L/76800.0L*x_0 - 1.0L/2160.0L*pow(x_1, 6) + (79.0L/7200.0L)*pow(x_1, 5) - 17.0L/160.0L*pow(x_1, 4) + (4531.0L/8640.0L)*pow(x_1, 3) - 1247.0L/960.0L*pow(x_1, 2) + (90257.0L/57600.0L)*x_1 - 731563.0L/921600.0L;
       return __pp_r561___result;
    }
    static inline O __pp_r562__(const I &x_0, const I &x_1) {
       O __pp_r562___result;
       __pp_r562___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 - 157.0L/4800.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (139.0L/1440.0L)*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/10.0L*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 - 5897.0L/5760.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (19.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (4639.0L/2880.0L)*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) - 13.0L/240.0L*x_0*pow(x_1, 4) + (73.0L/180.0L)*x_0*pow(x_1, 3) - 219.0L/160.0L*x_0*pow(x_1, 2) + (28271.0L/11520.0L)*x_0*x_1 - 162857.0L/76800.0L*x_0 - 1.0L/2160.0L*pow(x_1, 6) + (79.0L/7200.0L)*pow(x_1, 5) - 17.0L/160.0L*pow(x_1, 4) + (4531.0L/8640.0L)*pow(x_1, 3) - 1247.0L/960.0L*pow(x_1, 2) + (90257.0L/57600.0L)*x_1 - 731563.0L/921600.0L;
       return __pp_r562___result;
    }
    static inline O __pp_r563__(const I &x_0, const I &x_1) {
       O __pp_r563___result;
       __pp_r563___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 - 157.0L/4800.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (139.0L/1440.0L)*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/10.0L*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 - 5897.0L/5760.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (19.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (4639.0L/2880.0L)*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) - 13.0L/240.0L*x_0*pow(x_1, 4) + (73.0L/180.0L)*x_0*pow(x_1, 3) - 219.0L/160.0L*x_0*pow(x_1, 2) + (28271.0L/11520.0L)*x_0*x_1 - 162857.0L/76800.0L*x_0 - 1.0L/2160.0L*pow(x_1, 6) + (79.0L/7200.0L)*pow(x_1, 5) - 17.0L/160.0L*pow(x_1, 4) + (4531.0L/8640.0L)*pow(x_1, 3) - 1247.0L/960.0L*pow(x_1, 2) + (90257.0L/57600.0L)*x_1 - 731563.0L/921600.0L;
       return __pp_r563___result;
    }
    static inline O __pp_r564__(const I &x_0, const I &x_1) {
       O __pp_r564___result;
       __pp_r564___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/144.0L)*pow(x_0, 5)*x_1 - 157.0L/4800.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (139.0L/1440.0L)*pow(x_0, 4)*x_1 - 1003.0L/3840.0L*pow(x_0, 4) + (1.0L/108.0L)*pow(x_0, 3)*pow(x_1, 3) - 1.0L/10.0L*pow(x_0, 3)*pow(x_1, 2) + (791.0L/1440.0L)*pow(x_0, 3)*x_1 - 5897.0L/5760.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) + (19.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 11.0L/20.0L*pow(x_0, 2)*pow(x_1, 2) + (4639.0L/2880.0L)*pow(x_0, 2)*x_1 - 32263.0L/15360.0L*pow(x_0, 2) + (1.0L/360.0L)*x_0*pow(x_1, 5) - 13.0L/240.0L*x_0*pow(x_1, 4) + (73.0L/180.0L)*x_0*pow(x_1, 3) - 219.0L/160.0L*x_0*pow(x_1, 2) + (28271.0L/11520.0L)*x_0*x_1 - 162857.0L/76800.0L*x_0 - 1.0L/2160.0L*pow(x_1, 6) + (79.0L/7200.0L)*pow(x_1, 5) - 17.0L/160.0L*pow(x_1, 4) + (4531.0L/8640.0L)*pow(x_1, 3) - 1247.0L/960.0L*pow(x_1, 2) + (90257.0L/57600.0L)*x_1 - 731563.0L/921600.0L;
       return __pp_r564___result;
    }
    static inline O __pp_r565__(const I &x_0, const I &x_1) {
       O __pp_r565___result;
       __pp_r565___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 - 311.0L/14400.0L*pow(x_0, 5) + (59.0L/1440.0L)*pow(x_0, 4)*x_1 - 1729.0L/11520.0L*pow(x_0, 4) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (151.0L/1440.0L)*pow(x_0, 3)*x_1 - 7451.0L/17280.0L*pow(x_0, 3) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) - 481.0L/2880.0L*pow(x_0, 2)*x_1 - 14869.0L/46080.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 4) - 7.0L/180.0L*x_0*pow(x_1, 3) + (589.0L/1440.0L)*x_0*pow(x_1, 2) - 12689.0L/11520.0L*x_0*x_1 + (166789.0L/230400.0L)*x_0 - 1.0L/7200.0L*pow(x_1, 5) + (7.0L/1440.0L)*pow(x_1, 4) - 589.0L/8640.0L*pow(x_1, 3) + (1379.0L/2880.0L)*pow(x_1, 2) - 73583.0L/57600.0L*x_1 + 3048191.0L/2764800.0L;
       return __pp_r565___result;
    }
    static inline O __pp_r566__(const I &x_0, const I &x_1) {
       O __pp_r566___result;
       __pp_r566___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 - 311.0L/14400.0L*pow(x_0, 5) + (59.0L/1440.0L)*pow(x_0, 4)*x_1 - 1729.0L/11520.0L*pow(x_0, 4) + (1.0L/90.0L)*pow(x_0, 3)*pow(x_1, 2) + (151.0L/1440.0L)*pow(x_0, 3)*x_1 - 7451.0L/17280.0L*pow(x_0, 3) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/60.0L)*pow(x_0, 2)*pow(x_1, 2) - 481.0L/2880.0L*pow(x_0, 2)*x_1 - 14869.0L/46080.0L*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 4) - 7.0L/180.0L*x_0*pow(x_1, 3) + (589.0L/1440.0L)*x_0*pow(x_1, 2) - 12689.0L/11520.0L*x_0*x_1 + (166789.0L/230400.0L)*x_0 - 1.0L/7200.0L*pow(x_1, 5) + (7.0L/1440.0L)*pow(x_1, 4) - 589.0L/8640.0L*pow(x_1, 3) + (1379.0L/2880.0L)*pow(x_1, 2) - 73583.0L/57600.0L*x_1 + 3048191.0L/2764800.0L;
       return __pp_r566___result;
    }
    static inline O __pp_r567__(const I &x_0, const I &x_1) {
       O __pp_r567___result;
       __pp_r567___result = (19.0L/43200.0L)*pow(x_0, 6) - 1.0L/3600.0L*pow(x_0, 5)*x_1 + (7.0L/960.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/32.0L*pow(x_0, 4)*x_1 + (65.0L/768.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/12.0L)*pow(x_0, 3)*pow(x_1, 2) - 35.0L/96.0L*pow(x_0, 3)*x_1 + (75.0L/128.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 325.0L/192.0L*pow(x_0, 2)*x_1 + (6625.0L/3072.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/96.0L)*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) + (75.0L/64.0L)*x_0*pow(x_1, 2) - 1375.0L/384.0L*x_0*x_1 + (12125.0L/3072.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 1.0L/960.0L*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) - 25.0L/128.0L*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) - 8875.0L/3072.0L*x_1 + 4375.0L/1536.0L;
       return __pp_r567___result;
    }
    static inline O __pp_r568__(const I &x_0, const I &x_1) {
       O __pp_r568___result;
       __pp_r568___result = (19.0L/43200.0L)*pow(x_0, 6) - 1.0L/3600.0L*pow(x_0, 5)*x_1 + (7.0L/960.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/32.0L*pow(x_0, 4)*x_1 + (65.0L/768.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/12.0L)*pow(x_0, 3)*pow(x_1, 2) - 35.0L/96.0L*pow(x_0, 3)*x_1 + (75.0L/128.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 325.0L/192.0L*pow(x_0, 2)*x_1 + (6625.0L/3072.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/96.0L)*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) + (75.0L/64.0L)*x_0*pow(x_1, 2) - 1375.0L/384.0L*x_0*x_1 + (12125.0L/3072.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 1.0L/960.0L*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) - 25.0L/128.0L*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) - 8875.0L/3072.0L*x_1 + 4375.0L/1536.0L;
       return __pp_r568___result;
    }
    static inline O __pp_r569__(const I &x_0, const I &x_1) {
       O __pp_r569___result;
       __pp_r569___result = (7.0L/3600.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 + (11.0L/600.0L)*pow(x_0, 5) + (1.0L/80.0L)*pow(x_0, 4)*pow(x_1, 2) - 53.0L/720.0L*pow(x_0, 4)*x_1 + (13.0L/240.0L)*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) + (13.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) - 13.0L/45.0L*pow(x_0, 3)*x_1 + (31.0L/1440.0L)*pow(x_0, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 4) - 4.0L/45.0L*pow(x_0, 2)*pow(x_1, 3) + (53.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 779.0L/1440.0L*pow(x_0, 2)*x_1 - 31.0L/256.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (19.0L/480.0L)*x_0*pow(x_1, 4) - 269.0L/1440.0L*x_0*pow(x_1, 3) + (373.0L/960.0L)*x_0*pow(x_1, 2) - 1117.0L/2304.0L*x_0*x_1 - 2021.0L/38400.0L*x_0 + (7.0L/14400.0L)*pow(x_1, 6) - 167.0L/14400.0L*pow(x_1, 5) + (203.0L/3840.0L)*pow(x_1, 4) - 1979.0L/17280.0L*pow(x_1, 3) + (319.0L/3072.0L)*pow(x_1, 2) - 40787.0L/230400.0L*x_1 + 151163.0L/921600.0L;
       return __pp_r569___result;
    }
    static inline O __pp_r570__(const I &x_0, const I &x_1) {
       O __pp_r570___result;
       __pp_r570___result = (83.0L/43200.0L)*pow(x_0, 6) - 3.0L/400.0L*pow(x_0, 5)*x_1 + (29.0L/1600.0L)*pow(x_0, 5) + (1.0L/90.0L)*pow(x_0, 4)*pow(x_1, 2) - 109.0L/1440.0L*pow(x_0, 4)*x_1 + (41.0L/768.0L)*pow(x_0, 4) - 1.0L/60.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/10.0L)*pow(x_0, 3)*pow(x_1, 2) - 85.0L/288.0L*pow(x_0, 3)*x_1 + (23.0L/1152.0L)*pow(x_0, 3) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 4) - 19.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) + (5.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 317.0L/576.0L*pow(x_0, 2)*x_1 - 629.0L/5120.0L*pow(x_0, 2) - 3.0L/400.0L*x_0*pow(x_1, 5) + (11.0L/480.0L)*x_0*pow(x_1, 4) - 61.0L/288.0L*x_0*pow(x_1, 3) + (71.0L/192.0L)*x_0*pow(x_1, 2) - 2833.0L/5760.0L*x_0*x_1 - 4123.0L/76800.0L*x_0 - 43.0L/43200.0L*pow(x_1, 6) - 263.0L/14400.0L*pow(x_1, 5) + (31.0L/768.0L)*pow(x_1, 4) - 439.0L/3456.0L*pow(x_1, 3) + (1487.0L/15360.0L)*pow(x_1, 2) - 41273.0L/230400.0L*x_1 + 3773.0L/23040.0L;
       return __pp_r570___result;
    }
    static inline O __pp_r571__(const I &x_0, const I &x_1) {
       O __pp_r571___result;
       __pp_r571___result = (31.0L/14400.0L)*pow(x_0, 6) - 11.0L/1800.0L*pow(x_0, 5)*x_1 + (301.0L/14400.0L)*pow(x_0, 5) + (7.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 89.0L/1440.0L*pow(x_0, 4)*x_1 + (155.0L/2304.0L)*pow(x_0, 4) - 13.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (23.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) - 23.0L/96.0L*pow(x_0, 3)*x_1 + (197.0L/3456.0L)*pow(x_0, 3) + (1.0L/160.0L)*pow(x_0, 2)*pow(x_1, 4) - 7.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) + (19.0L/48.0L)*pow(x_0, 2)*pow(x_1, 2) - 253.0L/576.0L*pow(x_0, 2)*x_1 - 3101.0L/46080.0L*pow(x_0, 2) - 11.0L/1800.0L*x_0*pow(x_1, 5) + (53.0L/1440.0L)*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) + (277.0L/576.0L)*x_0*pow(x_1, 2) - 731.0L/1920.0L*x_0*x_1 - 2129.0L/230400.0L*x_0 - 11.0L/14400.0L*pow(x_1, 6) - 223.0L/14400.0L*pow(x_1, 5) + (125.0L/2304.0L)*pow(x_1, 4) - 311.0L/3456.0L*pow(x_1, 3) + (7021.0L/46080.0L)*pow(x_1, 2) - 31033.0L/230400.0L*x_1 + 12343.0L/69120.0L;
       return __pp_r571___result;
    }
    static inline O __pp_r572__(const I &x_0, const I &x_1) {
       O __pp_r572___result;
       __pp_r572___result = (47.0L/21600.0L)*pow(x_0, 6) - 7.0L/1200.0L*pow(x_0, 5)*x_1 + (17.0L/800.0L)*pow(x_0, 5) + (23.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 7.0L/120.0L*pow(x_0, 4)*x_1 + (5.0L/72.0L)*pow(x_0, 4) - 1.0L/120.0L*pow(x_0, 3)*pow(x_1, 3) + (17.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) - 2.0L/9.0L*pow(x_0, 3)*x_1 + (37.0L/576.0L)*pow(x_0, 3) + (17.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/20.0L*pow(x_0, 2)*pow(x_1, 3) + (43.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 19.0L/48.0L*pow(x_0, 2)*x_1 - 619.0L/11520.0L*pow(x_0, 2) - 1.0L/600.0L*x_0*pow(x_1, 5) + (31.0L/480.0L)*x_0*pow(x_1, 4) - 25.0L/288.0L*x_0*pow(x_1, 3) + (109.0L/192.0L)*x_0*pow(x_1, 2) - 3761.0L/11520.0L*x_0*x_1 + (83.0L/19200.0L)*x_0 + (31.0L/43200.0L)*pow(x_1, 6) - 7.0L/1600.0L*pow(x_1, 5) + (205.0L/2304.0L)*pow(x_1, 4) - 37.0L/1152.0L*pow(x_1, 3) + (9521.0L/46080.0L)*pow(x_1, 2) - 8261.0L/76800.0L*x_1 + 101869.0L/552960.0L;
       return __pp_r572___result;
    }
    static inline O __pp_r573__(const I &x_0, const I &x_1) {
       O __pp_r573___result;
       __pp_r573___result = -11.0L/10800.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 - 13.0L/600.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (19.0L/720.0L)*pow(x_0, 4)*x_1 - 41.0L/240.0L*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (29.0L/180.0L)*pow(x_0, 3)*x_1 - 941.0L/1440.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 7.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 1.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (679.0L/1440.0L)*pow(x_0, 2)*x_1 - 1613.0L/1280.0L*pow(x_0, 2) - 1.0L/400.0L*x_0*pow(x_1, 5) + (13.0L/480.0L)*x_0*pow(x_1, 4) - 107.0L/1440.0L*x_0*pow(x_1, 3) - 113.0L/960.0L*x_0*pow(x_1, 2) + (7537.0L/11520.0L)*x_0*x_1 - 41387.0L/38400.0L*x_0 + (19.0L/43200.0L)*pow(x_1, 6) - 149.0L/14400.0L*pow(x_1, 5) + (149.0L/3840.0L)*pow(x_1, 4) - 521.0L/17280.0L*pow(x_1, 3) - 2779.0L/15360.0L*pow(x_1, 2) + (77311.0L/230400.0L)*x_1 - 203131.0L/921600.0L;
       return __pp_r573___result;
    }
    static inline O __pp_r574__(const I &x_0, const I &x_1) {
       O __pp_r574___result;
       __pp_r574___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 - 7.0L/320.0L*pow(x_0, 5) + (7.0L/288.0L)*pow(x_0, 4)*x_1 - 659.0L/3840.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) + (223.0L/1440.0L)*pow(x_0, 3)*x_1 - 3773.0L/5760.0L*pow(x_0, 3) - 1.0L/18.0L*pow(x_0, 2)*pow(x_1, 3) - 1.0L/40.0L*pow(x_0, 2)*pow(x_1, 2) + (1331.0L/2880.0L)*pow(x_0, 2)*x_1 - 6461.0L/5120.0L*pow(x_0, 2) - 1.0L/144.0L*x_0*pow(x_1, 5) + (1.0L/96.0L)*x_0*pow(x_1, 4) - 143.0L/1440.0L*x_0*pow(x_1, 3) - 131.0L/960.0L*x_0*pow(x_1, 2) + (233.0L/360.0L)*x_0*x_1 - 16571.0L/15360.0L*x_0 - 1.0L/960.0L*pow(x_1, 6) - 49.0L/2880.0L*pow(x_1, 5) + (101.0L/3840.0L)*pow(x_1, 4) - 737.0L/17280.0L*pow(x_1, 3) - 2887.0L/15360.0L*pow(x_1, 2) + (3073.0L/9216.0L)*x_1 - 101687.0L/460800.0L;
       return __pp_r574___result;
    }
    static inline O __pp_r575__(const I &x_0, const I &x_1) {
       O __pp_r575___result;
       __pp_r575___result = -7.0L/8640.0L*pow(x_0, 6) + (1.0L/360.0L)*pow(x_0, 5)*x_1 - 11.0L/576.0L*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) + (11.0L/288.0L)*pow(x_0, 4)*x_1 - 1817.0L/11520.0L*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/36.0L)*pow(x_0, 3)*pow(x_1, 2) + (101.0L/480.0L)*pow(x_0, 3)*x_1 - 10679.0L/17280.0L*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/36.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) + (1651.0L/2880.0L)*pow(x_0, 2)*x_1 - 55589.0L/46080.0L*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) + (7.0L/288.0L)*x_0*pow(x_1, 4) - 7.0L/160.0L*x_0*pow(x_1, 3) - 73.0L/2880.0L*x_0*pow(x_1, 2) + (91.0L/120.0L)*x_0*x_1 - 9533.0L/9216.0L*x_0 - 7.0L/8640.0L*pow(x_1, 6) - 41.0L/2880.0L*pow(x_1, 5) + (463.0L/11520.0L)*pow(x_1, 4) - 97.0L/17280.0L*pow(x_1, 3) - 6101.0L/46080.0L*pow(x_1, 2) + (17413.0L/46080.0L)*x_1 - 284581.0L/1382400.0L;
       return __pp_r575___result;
    }
    static inline O __pp_r576__(const I &x_0, const I &x_1) {
       O __pp_r576___result;
       __pp_r576___result = -17.0L/21600.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 3.0L/160.0L*pow(x_0, 5) + (7.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) + (1.0L/24.0L)*pow(x_0, 4)*x_1 - 7.0L/45.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) + (41.0L/180.0L)*pow(x_0, 3)*x_1 - 1759.0L/2880.0L*pow(x_0, 3) + (13.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (53.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (37.0L/60.0L)*pow(x_0, 2)*x_1 - 13741.0L/11520.0L*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) + (5.0L/96.0L)*x_0*pow(x_1, 4) + (37.0L/1440.0L)*x_0*pow(x_1, 3) + (59.0L/960.0L)*x_0*pow(x_1, 2) + (9361.0L/11520.0L)*x_0*x_1 - 49.0L/48.0L*x_0 + (29.0L/43200.0L)*pow(x_1, 6) - 1.0L/320.0L*pow(x_1, 5) + (863.0L/11520.0L)*pow(x_1, 4) + (301.0L/5760.0L)*pow(x_1, 3) - 3601.0L/46080.0L*pow(x_1, 2) + (6221.0L/15360.0L)*x_1 - 553537.0L/2764800.0L;
       return __pp_r576___result;
    }
    static inline O __pp_r577__(const I &x_0, const I &x_1) {
       O __pp_r577___result;
       __pp_r577___result = (7.0L/3600.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 + (11.0L/600.0L)*pow(x_0, 5) + (1.0L/80.0L)*pow(x_0, 4)*pow(x_1, 2) - 53.0L/720.0L*pow(x_0, 4)*x_1 + (13.0L/240.0L)*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) + (13.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) - 13.0L/45.0L*pow(x_0, 3)*x_1 + (31.0L/1440.0L)*pow(x_0, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 4) - 4.0L/45.0L*pow(x_0, 2)*pow(x_1, 3) + (53.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 779.0L/1440.0L*pow(x_0, 2)*x_1 - 31.0L/256.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (19.0L/480.0L)*x_0*pow(x_1, 4) - 269.0L/1440.0L*x_0*pow(x_1, 3) + (373.0L/960.0L)*x_0*pow(x_1, 2) - 1117.0L/2304.0L*x_0*x_1 - 2021.0L/38400.0L*x_0 + (7.0L/14400.0L)*pow(x_1, 6) - 107.0L/14400.0L*pow(x_1, 5) + (163.0L/3840.0L)*pow(x_1, 4) - 1799.0L/17280.0L*pow(x_1, 3) + (101.0L/1024.0L)*pow(x_1, 2) - 40487.0L/230400.0L*x_1 + 151043.0L/921600.0L;
       return __pp_r577___result;
    }
    static inline O __pp_r578__(const I &x_0, const I &x_1) {
       O __pp_r578___result;
       __pp_r578___result = (7.0L/3600.0L)*pow(x_0, 6) - 13.0L/1800.0L*pow(x_0, 5)*x_1 + (11.0L/600.0L)*pow(x_0, 5) + (1.0L/80.0L)*pow(x_0, 4)*pow(x_1, 2) - 53.0L/720.0L*pow(x_0, 4)*x_1 + (13.0L/240.0L)*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) + (13.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) - 13.0L/45.0L*pow(x_0, 3)*x_1 + (31.0L/1440.0L)*pow(x_0, 3) + (1.0L/120.0L)*pow(x_0, 2)*pow(x_1, 4) - 4.0L/45.0L*pow(x_0, 2)*pow(x_1, 3) + (53.0L/160.0L)*pow(x_0, 2)*pow(x_1, 2) - 779.0L/1440.0L*pow(x_0, 2)*x_1 - 31.0L/256.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) + (19.0L/480.0L)*x_0*pow(x_1, 4) - 269.0L/1440.0L*x_0*pow(x_1, 3) + (373.0L/960.0L)*x_0*pow(x_1, 2) - 1117.0L/2304.0L*x_0*x_1 - 2021.0L/38400.0L*x_0 + (7.0L/14400.0L)*pow(x_1, 6) - 107.0L/14400.0L*pow(x_1, 5) + (163.0L/3840.0L)*pow(x_1, 4) - 1799.0L/17280.0L*pow(x_1, 3) + (101.0L/1024.0L)*pow(x_1, 2) - 40487.0L/230400.0L*x_1 + 151043.0L/921600.0L;
       return __pp_r578___result;
    }
    static inline O __pp_r579__(const I &x_0, const I &x_1) {
       O __pp_r579___result;
       __pp_r579___result = -11.0L/10800.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 - 13.0L/600.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (19.0L/720.0L)*pow(x_0, 4)*x_1 - 41.0L/240.0L*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (29.0L/180.0L)*pow(x_0, 3)*x_1 - 941.0L/1440.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 7.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 1.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (679.0L/1440.0L)*pow(x_0, 2)*x_1 - 1613.0L/1280.0L*pow(x_0, 2) - 1.0L/400.0L*x_0*pow(x_1, 5) + (13.0L/480.0L)*x_0*pow(x_1, 4) - 107.0L/1440.0L*x_0*pow(x_1, 3) - 113.0L/960.0L*x_0*pow(x_1, 2) + (7537.0L/11520.0L)*x_0*x_1 - 41387.0L/38400.0L*x_0 + (19.0L/43200.0L)*pow(x_1, 6) - 89.0L/14400.0L*pow(x_1, 5) + (109.0L/3840.0L)*pow(x_1, 4) - 341.0L/17280.0L*pow(x_1, 3) - 953.0L/5120.0L*pow(x_1, 2) + (77611.0L/230400.0L)*x_1 - 203251.0L/921600.0L;
       return __pp_r579___result;
    }
    static inline O __pp_r580__(const I &x_0, const I &x_1) {
       O __pp_r580___result;
       __pp_r580___result = -11.0L/10800.0L*pow(x_0, 6) + (1.0L/600.0L)*pow(x_0, 5)*x_1 - 13.0L/600.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) + (19.0L/720.0L)*pow(x_0, 4)*x_1 - 41.0L/240.0L*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 2) + (29.0L/180.0L)*pow(x_0, 3)*x_1 - 941.0L/1440.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 7.0L/180.0L*pow(x_0, 2)*pow(x_1, 3) - 1.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (679.0L/1440.0L)*pow(x_0, 2)*x_1 - 1613.0L/1280.0L*pow(x_0, 2) - 1.0L/400.0L*x_0*pow(x_1, 5) + (13.0L/480.0L)*x_0*pow(x_1, 4) - 107.0L/1440.0L*x_0*pow(x_1, 3) - 113.0L/960.0L*x_0*pow(x_1, 2) + (7537.0L/11520.0L)*x_0*x_1 - 41387.0L/38400.0L*x_0 + (19.0L/43200.0L)*pow(x_1, 6) - 89.0L/14400.0L*pow(x_1, 5) + (109.0L/3840.0L)*pow(x_1, 4) - 341.0L/17280.0L*pow(x_1, 3) - 953.0L/5120.0L*pow(x_1, 2) + (77611.0L/230400.0L)*x_1 - 203251.0L/921600.0L;
       return __pp_r580___result;
    }
    static inline O __pp_r581__(const I &x_0, const I &x_1) {
       O __pp_r581___result;
       __pp_r581___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 - 41.0L/1200.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (4.0L/45.0L)*pow(x_0, 4)*x_1 - 127.0L/480.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 - 1481.0L/1440.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (31.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (2299.0L/1440.0L)*pow(x_0, 2)*x_1 - 2693.0L/1280.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 17.0L/480.0L*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) - 1193.0L/960.0L*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 - 80267.0L/38400.0L*x_0 - 11.0L/43200.0L*pow(x_1, 6) + (91.0L/14400.0L)*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) + (6139.0L/17280.0L)*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) + (310891.0L/230400.0L)*x_1 - 669811.0L/921600.0L;
       return __pp_r581___result;
    }
    static inline O __pp_r582__(const I &x_0, const I &x_1) {
       O __pp_r582___result;
       __pp_r582___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 - 41.0L/1200.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (4.0L/45.0L)*pow(x_0, 4)*x_1 - 127.0L/480.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 - 1481.0L/1440.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (31.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (2299.0L/1440.0L)*pow(x_0, 2)*x_1 - 2693.0L/1280.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 17.0L/480.0L*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) - 1193.0L/960.0L*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 - 80267.0L/38400.0L*x_0 - 11.0L/43200.0L*pow(x_1, 6) + (91.0L/14400.0L)*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) + (6139.0L/17280.0L)*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) + (310891.0L/230400.0L)*x_1 - 669811.0L/921600.0L;
       return __pp_r582___result;
    }
    static inline O __pp_r583__(const I &x_0, const I &x_1) {
       O __pp_r583___result;
       __pp_r583___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 - 41.0L/1200.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (4.0L/45.0L)*pow(x_0, 4)*x_1 - 127.0L/480.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 - 1481.0L/1440.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (31.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (2299.0L/1440.0L)*pow(x_0, 2)*x_1 - 2693.0L/1280.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 17.0L/480.0L*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) - 1193.0L/960.0L*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 - 80267.0L/38400.0L*x_0 - 11.0L/43200.0L*pow(x_1, 6) + (91.0L/14400.0L)*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) + (6139.0L/17280.0L)*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) + (310891.0L/230400.0L)*x_1 - 669811.0L/921600.0L;
       return __pp_r583___result;
    }
    static inline O __pp_r584__(const I &x_0, const I &x_1) {
       O __pp_r584___result;
       __pp_r584___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 - 41.0L/1200.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (4.0L/45.0L)*pow(x_0, 4)*x_1 - 127.0L/480.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 - 1481.0L/1440.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (31.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (2299.0L/1440.0L)*pow(x_0, 2)*x_1 - 2693.0L/1280.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 17.0L/480.0L*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) - 1193.0L/960.0L*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 - 80267.0L/38400.0L*x_0 - 11.0L/43200.0L*pow(x_1, 6) + (91.0L/14400.0L)*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) + (6139.0L/17280.0L)*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) + (310891.0L/230400.0L)*x_1 - 669811.0L/921600.0L;
       return __pp_r584___result;
    }
    static inline O __pp_r585__(const I &x_0, const I &x_1) {
       O __pp_r585___result;
       __pp_r585___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 - 3.0L/100.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (4.0L/45.0L)*pow(x_0, 4)*x_1 - 17.0L/80.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 - 553.0L/720.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (31.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (2299.0L/1440.0L)*pow(x_0, 2)*x_1 - 5579.0L/3840.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 17.0L/480.0L*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) - 1193.0L/960.0L*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 - 16339.0L/12800.0L*x_0 - 11.0L/43200.0L*pow(x_1, 6) + (91.0L/14400.0L)*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) + (6139.0L/17280.0L)*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) + (310891.0L/230400.0L)*x_1 - 294811.0L/921600.0L;
       return __pp_r585___result;
    }
    static inline O __pp_r586__(const I &x_0, const I &x_1) {
       O __pp_r586___result;
       __pp_r586___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 + (1.0L/36.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 4)*x_1 + (37.0L/144.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/36.0L)*pow(x_0, 3)*pow(x_1, 2) - 29.0L/72.0L*pow(x_0, 3)*x_1 + (547.0L/432.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 419.0L/288.0L*pow(x_0, 2)*x_1 + (8077.0L/2304.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 5.0L/288.0L*x_0*pow(x_1, 4) + (19.0L/288.0L)*x_0*pow(x_1, 3) + (163.0L/576.0L)*x_0*pow(x_1, 2) - 6029.0L/2304.0L*x_0*x_1 + (119107.0L/23040.0L)*x_0 - 1.0L/4800.0L*pow(x_1, 6) + (13.0L/2880.0L)*pow(x_1, 5) - 83.0L/2304.0L*pow(x_1, 4) + (349.0L/3456.0L)*pow(x_1, 3) + (1933.0L/9216.0L)*pow(x_1, 2) - 86339.0L/46080.0L*x_1 + 1753837.0L/552960.0L;
       return __pp_r586___result;
    }
    static inline O __pp_r587__(const I &x_0, const I &x_1) {
       O __pp_r587___result;
       __pp_r587___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 + (1.0L/30.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/12.0L*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/12.0L)*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 + (25.0L/16.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/32.0L*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/96.0L)*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) + (75.0L/64.0L)*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 + (3375.0L/512.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 1.0L/960.0L*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) - 25.0L/128.0L*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) - 3375.0L/1024.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r587___result;
    }
    static inline O __pp_r588__(const I &x_0, const I &x_1) {
       O __pp_r588___result;
       __pp_r588___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 - 3.0L/100.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (4.0L/45.0L)*pow(x_0, 4)*x_1 - 17.0L/80.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 - 553.0L/720.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (31.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (2299.0L/1440.0L)*pow(x_0, 2)*x_1 - 5579.0L/3840.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 17.0L/480.0L*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) - 1193.0L/960.0L*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 - 16339.0L/12800.0L*x_0 - 11.0L/43200.0L*pow(x_1, 6) + (91.0L/14400.0L)*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) + (6139.0L/17280.0L)*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) + (310891.0L/230400.0L)*x_1 - 294811.0L/921600.0L;
       return __pp_r588___result;
    }
    static inline O __pp_r589__(const I &x_0, const I &x_1) {
       O __pp_r589___result;
       __pp_r589___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 + (1.0L/36.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 4)*x_1 + (37.0L/144.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/36.0L)*pow(x_0, 3)*pow(x_1, 2) - 29.0L/72.0L*pow(x_0, 3)*x_1 + (547.0L/432.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 419.0L/288.0L*pow(x_0, 2)*x_1 + (8077.0L/2304.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 5.0L/288.0L*x_0*pow(x_1, 4) + (19.0L/288.0L)*x_0*pow(x_1, 3) + (163.0L/576.0L)*x_0*pow(x_1, 2) - 6029.0L/2304.0L*x_0*x_1 + (119107.0L/23040.0L)*x_0 - 1.0L/4800.0L*pow(x_1, 6) + (13.0L/2880.0L)*pow(x_1, 5) - 83.0L/2304.0L*pow(x_1, 4) + (349.0L/3456.0L)*pow(x_1, 3) + (1933.0L/9216.0L)*pow(x_1, 2) - 86339.0L/46080.0L*x_1 + 1753837.0L/552960.0L;
       return __pp_r589___result;
    }
    static inline O __pp_r590__(const I &x_0, const I &x_1) {
       O __pp_r590___result;
       __pp_r590___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 + (1.0L/30.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/12.0L*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/12.0L)*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 + (25.0L/16.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/32.0L*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/96.0L)*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) + (75.0L/64.0L)*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 + (3375.0L/512.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 1.0L/960.0L*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) - 25.0L/128.0L*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) - 3375.0L/1024.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r590___result;
    }
    static inline O __pp_r591__(const I &x_0, const I &x_1) {
       O __pp_r591___result;
       __pp_r591___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 + (1.0L/30.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/12.0L*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/12.0L)*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 + (25.0L/16.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/32.0L*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/96.0L)*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) + (75.0L/64.0L)*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 + (3375.0L/512.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 1.0L/960.0L*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) - 25.0L/128.0L*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) - 3375.0L/1024.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r591___result;
    }
    static inline O __pp_r592__(const I &x_0, const I &x_1) {
       O __pp_r592___result;
       __pp_r592___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 + (1.0L/30.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/12.0L*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/12.0L)*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 + (25.0L/16.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/32.0L*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/96.0L)*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) + (75.0L/64.0L)*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 + (3375.0L/512.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 1.0L/960.0L*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) - 25.0L/128.0L*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) - 3375.0L/1024.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r592___result;
    }
    static inline O __pp_r593__(const I &x_0, const I &x_1) {
       O __pp_r593___result;
       __pp_r593___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 + (1.0L/30.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/12.0L*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/12.0L)*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 + (25.0L/16.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/32.0L*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/96.0L)*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) + (75.0L/64.0L)*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 + (3375.0L/512.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 1.0L/960.0L*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) - 25.0L/128.0L*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) - 3375.0L/1024.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r593___result;
    }
    static inline O __pp_r594__(const I &x_0, const I &x_1) {
       O __pp_r594___result;
       __pp_r594___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 + (1.0L/30.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/12.0L*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/12.0L)*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 + (25.0L/16.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/32.0L*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/96.0L)*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) + (75.0L/64.0L)*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 + (3375.0L/512.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 1.0L/960.0L*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) - 25.0L/128.0L*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) - 3375.0L/1024.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r594___result;
    }
    static inline O __pp_r595__(const I &x_0, const I &x_1) {
       O __pp_r595___result;
       __pp_r595___result = (1.0L/4320.0L)*pow(x_0, 6) - 1.0L/720.0L*pow(x_0, 5)*x_1 + (1.0L/360.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/72.0L*pow(x_0, 4)*x_1 + (13.0L/720.0L)*pow(x_0, 4) - 1.0L/216.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/36.0L)*pow(x_0, 3)*pow(x_1, 2) - 1.0L/20.0L*pow(x_0, 3)*x_1 + (53.0L/540.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/36.0L*pow(x_0, 2)*pow(x_1, 3) + (7.0L/120.0L)*pow(x_0, 2)*pow(x_1, 2) - 7.0L/90.0L*pow(x_0, 2)*x_1 + (137.0L/360.0L)*pow(x_0, 2) - 1.0L/180.0L*x_0*pow(x_1, 5) + (7.0L/288.0L)*x_0*pow(x_1, 4) - 7.0L/160.0L*x_0*pow(x_1, 3) - 73.0L/2880.0L*x_0*pow(x_1, 2) - 71.0L/1280.0L*x_0*x_1 + (3671.0L/4608.0L)*x_0 - 7.0L/8640.0L*pow(x_1, 6) - 41.0L/2880.0L*pow(x_1, 5) + (463.0L/11520.0L)*pow(x_1, 4) - 97.0L/17280.0L*pow(x_1, 3) - 6101.0L/46080.0L*pow(x_1, 2) - 1337.0L/46080.0L*x_1 + 1821463.0L/2764800.0L;
       return __pp_r595___result;
    }
    static inline O __pp_r596__(const I &x_0, const I &x_1) {
       O __pp_r596___result;
       __pp_r596___result = -1.0L/800.0L*pow(x_0, 6) + (11.0L/3600.0L)*pow(x_0, 5)*x_1 - 13.0L/600.0L*pow(x_0, 5) - 1.0L/480.0L*pow(x_0, 4)*pow(x_1, 2) + (17.0L/360.0L)*pow(x_0, 4)*x_1 - 3.0L/20.0L*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/30.0L*pow(x_0, 3)*pow(x_1, 2) + (103.0L/360.0L)*pow(x_0, 3)*x_1 - 373.0L/720.0L*pow(x_0, 3) + (1.0L/480.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 31.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (1219.0L/1440.0L)*pow(x_0, 2)*x_1 - 3419.0L/3840.0L*pow(x_0, 2) - 19.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/60.0L)*x_0*pow(x_1, 4) + (29.0L/720.0L)*x_0*pow(x_1, 3) - 39.0L/80.0L*x_0*pow(x_1, 2) + (7001.0L/5760.0L)*x_0*x_1 - 5773.0L/9600.0L*x_0 - 1.0L/1200.0L*pow(x_1, 6) - 97.0L/7200.0L*pow(x_1, 5) + (19.0L/640.0L)*pow(x_1, 4) + (617.0L/8640.0L)*pow(x_1, 3) - 3457.0L/7680.0L*pow(x_1, 2) + (77183.0L/115200.0L)*x_1 + 8317.0L/460800.0L;
       return __pp_r596___result;
    }
    static inline O __pp_r597__(const I &x_0, const I &x_1) {
       O __pp_r597___result;
       __pp_r597___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 - 3.0L/100.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (4.0L/45.0L)*pow(x_0, 4)*x_1 - 17.0L/80.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 - 553.0L/720.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (31.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (2299.0L/1440.0L)*pow(x_0, 2)*x_1 - 5579.0L/3840.0L*pow(x_0, 2) - 1.0L/400.0L*x_0*pow(x_1, 5) - 1.0L/40.0L*x_0*pow(x_1, 4) + (209.0L/720.0L)*x_0*pow(x_1, 3) - 99.0L/80.0L*x_0*pow(x_1, 2) + (13481.0L/5760.0L)*x_0*x_1 - 12253.0L/9600.0L*x_0 - 7.0L/5400.0L*pow(x_1, 6) - 37.0L/7200.0L*pow(x_1, 5) - 21.0L/640.0L*pow(x_1, 4) + (2777.0L/8640.0L)*pow(x_1, 3) - 7777.0L/7680.0L*pow(x_1, 2) + (154943.0L/115200.0L)*x_1 - 147203.0L/460800.0L;
       return __pp_r597___result;
    }
    static inline O __pp_r598__(const I &x_0, const I &x_1) {
       O __pp_r598___result;
       __pp_r598___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 - 3.0L/100.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (4.0L/45.0L)*pow(x_0, 4)*x_1 - 17.0L/80.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 - 553.0L/720.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (31.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (2299.0L/1440.0L)*pow(x_0, 2)*x_1 - 5579.0L/3840.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 17.0L/480.0L*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) - 1193.0L/960.0L*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 - 16339.0L/12800.0L*x_0 - 11.0L/43200.0L*pow(x_1, 6) + (91.0L/14400.0L)*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) + (6139.0L/17280.0L)*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) + (310891.0L/230400.0L)*x_1 - 294811.0L/921600.0L;
       return __pp_r598___result;
    }
    static inline O __pp_r599__(const I &x_0, const I &x_1) {
       O __pp_r599___result;
       __pp_r599___result = -37.0L/21600.0L*pow(x_0, 6) + (7.0L/1200.0L)*pow(x_0, 5)*x_1 - 3.0L/100.0L*pow(x_0, 5) - 13.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (4.0L/45.0L)*pow(x_0, 4)*x_1 - 17.0L/80.0L*pow(x_0, 4) + (1.0L/120.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/60.0L*pow(x_0, 3)*pow(x_1, 2) + (193.0L/360.0L)*pow(x_0, 3)*x_1 - 553.0L/720.0L*pow(x_0, 3) - 7.0L/1440.0L*pow(x_0, 2)*pow(x_1, 4) + (31.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 91.0L/160.0L*pow(x_0, 2)*pow(x_1, 2) + (2299.0L/1440.0L)*pow(x_0, 2)*x_1 - 5579.0L/3840.0L*pow(x_0, 2) + (1.0L/600.0L)*x_0*pow(x_1, 5) - 17.0L/480.0L*x_0*pow(x_1, 4) + (433.0L/1440.0L)*x_0*pow(x_1, 3) - 1193.0L/960.0L*x_0*pow(x_1, 2) + (26977.0L/11520.0L)*x_0*x_1 - 16339.0L/12800.0L*x_0 - 11.0L/43200.0L*pow(x_1, 6) + (91.0L/14400.0L)*pow(x_1, 5) - 251.0L/3840.0L*pow(x_1, 4) + (6139.0L/17280.0L)*pow(x_1, 3) - 5273.0L/5120.0L*pow(x_1, 2) + (310891.0L/230400.0L)*x_1 - 294811.0L/921600.0L;
       return __pp_r599___result;
    }
    static inline O __pp_r600__(const I &x_0, const I &x_1) {
       O __pp_r600___result;
       __pp_r600___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 + (1.0L/36.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 4)*x_1 + (37.0L/144.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/36.0L)*pow(x_0, 3)*pow(x_1, 2) - 29.0L/72.0L*pow(x_0, 3)*x_1 + (547.0L/432.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 419.0L/288.0L*pow(x_0, 2)*x_1 + (8077.0L/2304.0L)*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 1.0L/144.0L*x_0*pow(x_1, 4) + (1.0L/18.0L)*x_0*pow(x_1, 3) + (83.0L/288.0L)*x_0*pow(x_1, 2) - 377.0L/144.0L*x_0*x_1 + (11911.0L/2304.0L)*x_0 - 1.0L/800.0L*pow(x_1, 6) - 1.0L/144.0L*pow(x_1, 5) - 1.0L/288.0L*pow(x_1, 4) + (29.0L/432.0L)*pow(x_1, 3) + (523.0L/2304.0L)*pow(x_1, 2) - 4327.0L/2304.0L*x_1 + 10963.0L/3456.0L;
       return __pp_r600___result;
    }
    static inline O __pp_r601__(const I &x_0, const I &x_1) {
       O __pp_r601___result;
       __pp_r601___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 + (1.0L/36.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 4)*x_1 + (37.0L/144.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/36.0L)*pow(x_0, 3)*pow(x_1, 2) - 29.0L/72.0L*pow(x_0, 3)*x_1 + (547.0L/432.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 419.0L/288.0L*pow(x_0, 2)*x_1 + (8077.0L/2304.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 5.0L/288.0L*x_0*pow(x_1, 4) + (19.0L/288.0L)*x_0*pow(x_1, 3) + (163.0L/576.0L)*x_0*pow(x_1, 2) - 6029.0L/2304.0L*x_0*x_1 + (119107.0L/23040.0L)*x_0 - 1.0L/4800.0L*pow(x_1, 6) + (13.0L/2880.0L)*pow(x_1, 5) - 83.0L/2304.0L*pow(x_1, 4) + (349.0L/3456.0L)*pow(x_1, 3) + (1933.0L/9216.0L)*pow(x_1, 2) - 86339.0L/46080.0L*x_1 + 1753837.0L/552960.0L;
       return __pp_r601___result;
    }
    static inline O __pp_r602__(const I &x_0, const I &x_1) {
       O __pp_r602___result;
       __pp_r602___result = (1.0L/800.0L)*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 + (1.0L/36.0L)*pow(x_0, 5) + (1.0L/480.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/18.0L*pow(x_0, 4)*x_1 + (37.0L/144.0L)*pow(x_0, 4) + (1.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/36.0L)*pow(x_0, 3)*pow(x_1, 2) - 29.0L/72.0L*pow(x_0, 3)*x_1 + (547.0L/432.0L)*pow(x_0, 3) - 1.0L/480.0L*pow(x_0, 2)*pow(x_1, 4) + (1.0L/72.0L)*pow(x_0, 2)*pow(x_1, 3) + (13.0L/96.0L)*pow(x_0, 2)*pow(x_1, 2) - 419.0L/288.0L*pow(x_0, 2)*x_1 + (8077.0L/2304.0L)*pow(x_0, 2) + (1.0L/900.0L)*x_0*pow(x_1, 5) - 5.0L/288.0L*x_0*pow(x_1, 4) + (19.0L/288.0L)*x_0*pow(x_1, 3) + (163.0L/576.0L)*x_0*pow(x_1, 2) - 6029.0L/2304.0L*x_0*x_1 + (119107.0L/23040.0L)*x_0 - 1.0L/4800.0L*pow(x_1, 6) + (13.0L/2880.0L)*pow(x_1, 5) - 83.0L/2304.0L*pow(x_1, 4) + (349.0L/3456.0L)*pow(x_1, 3) + (1933.0L/9216.0L)*pow(x_1, 2) - 86339.0L/46080.0L*x_1 + 1753837.0L/552960.0L;
       return __pp_r602___result;
    }
    static inline O __pp_r603__(const I &x_0, const I &x_1) {
       O __pp_r603___result;
       __pp_r603___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 + (1.0L/30.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/12.0L*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/12.0L)*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 + (25.0L/16.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/32.0L*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/96.0L)*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) + (75.0L/64.0L)*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 + (3375.0L/512.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 1.0L/960.0L*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) - 25.0L/128.0L*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) - 3375.0L/1024.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r603___result;
    }
    static inline O __pp_r604__(const I &x_0, const I &x_1) {
       O __pp_r604___result;
       __pp_r604___result = (1.0L/675.0L)*pow(x_0, 6) - 1.0L/225.0L*pow(x_0, 5)*x_1 + (1.0L/30.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/12.0L*pow(x_0, 4)*x_1 + (5.0L/16.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/12.0L)*pow(x_0, 3)*pow(x_1, 2) - 5.0L/8.0L*pow(x_0, 3)*x_1 + (25.0L/16.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) + (15.0L/32.0L)*pow(x_0, 2)*pow(x_1, 2) - 75.0L/32.0L*pow(x_0, 2)*x_1 + (1125.0L/256.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (1.0L/96.0L)*x_0*pow(x_1, 4) - 5.0L/32.0L*x_0*pow(x_1, 3) + (75.0L/64.0L)*x_0*pow(x_1, 2) - 1125.0L/256.0L*x_0*x_1 + (3375.0L/512.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 1.0L/960.0L*pow(x_1, 5) + (5.0L/256.0L)*pow(x_1, 4) - 25.0L/128.0L*pow(x_1, 3) + (1125.0L/1024.0L)*pow(x_1, 2) - 3375.0L/1024.0L*x_1 + 16875.0L/4096.0L;
       return __pp_r604___result;
    }
    static inline O __pp_r605__(const I &x_0, const I &x_1) {
       O __pp_r605___result;
       __pp_r605___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 - 1.0L/960.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/32.0L*pow(x_0, 4)*x_1 - 127.0L/11520.0L*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) - 227.0L/1440.0L*pow(x_0, 3)*x_1 - 113.0L/5760.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) - 37.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 131.0L/320.0L*pow(x_0, 2)*x_1 + (5201.0L/46080.0L)*pow(x_0, 2) + (7.0L/1200.0L)*x_0*pow(x_1, 5) + (5.0L/96.0L)*x_0*pow(x_1, 4) - 113.0L/1440.0L*x_0*pow(x_1, 3) - 97.0L/320.0L*x_0*pow(x_1, 2) - 101.0L/180.0L*x_0*x_1 + (1453.0L/3072.0L)*x_0 + (109.0L/43200.0L)*pow(x_1, 6) + (9.0L/320.0L)*pow(x_1, 5) + (1433.0L/11520.0L)*pow(x_1, 4) + (31.0L/5760.0L)*pow(x_1, 3) - 14791.0L/46080.0L*pow(x_1, 2) - 1689.0L/5120.0L*x_1 + 686569.0L/1382400.0L;
       return __pp_r605___result;
    }
    static inline O __pp_r606__(const I &x_0, const I &x_1) {
       O __pp_r606___result;
       __pp_r606___result = -1.0L/43200.0L*pow(x_0, 6) - 11.0L/3600.0L*pow(x_0, 5)*x_1 - 11.0L/4800.0L*pow(x_0, 5) - 1.0L/720.0L*pow(x_0, 4)*pow(x_1, 2) - 7.0L/160.0L*pow(x_0, 4)*x_1 - 289.0L/11520.0L*pow(x_0, 4) - 7.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/20.0L*pow(x_0, 3)*pow(x_1, 2) - 389.0L/1440.0L*pow(x_0, 3)*x_1 - 599.0L/5760.0L*pow(x_0, 3) - 1.0L/180.0L*pow(x_0, 2)*pow(x_1, 4) - 17.0L/120.0L*pow(x_0, 2)*pow(x_1, 3) - 199.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 293.0L/320.0L*pow(x_0, 2)*x_1 - 7921.0L/46080.0L*pow(x_0, 2) - 11.0L/3600.0L*x_0*pow(x_1, 5) - 23.0L/480.0L*x_0*pow(x_1, 4) - 761.0L/1440.0L*x_0*pow(x_1, 3) - 421.0L/320.0L*x_0*pow(x_1, 2) - 9793.0L/5760.0L*x_0*x_1 - 3041.0L/76800.0L*x_0 - 19.0L/43200.0L*pow(x_1, 6) - 19.0L/1600.0L*pow(x_1, 5) - 1159.0L/11520.0L*pow(x_1, 4) - 3857.0L/5760.0L*pow(x_1, 3) - 67279.0L/46080.0L*pow(x_1, 2) - 34689.0L/25600.0L*x_1 + 19391.0L/172800.0L;
       return __pp_r606___result;
    }
    static inline O __pp_r607__(const I &x_0, const I &x_1) {
       O __pp_r607___result;
       __pp_r607___result = (19.0L/43200.0L)*pow(x_0, 6) - 1.0L/3600.0L*pow(x_0, 5)*x_1 + (127.0L/14400.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (17.0L/1440.0L)*pow(x_0, 4)*x_1 + (991.0L/11520.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (11.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (251.0L/1440.0L)*pow(x_0, 3)*x_1 + (8443.0L/17280.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 11.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (121.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (2483.0L/2880.0L)*pow(x_0, 2)*x_1 + (73999.0L/46080.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (11.0L/1440.0L)*x_0*pow(x_1, 4) - 121.0L/1440.0L*x_0*pow(x_1, 3) + (1331.0L/2880.0L)*x_0*pow(x_1, 2) + (10687.0L/5760.0L)*x_0*x_1 + (646237.0L/230400.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 11.0L/14400.0L*pow(x_1, 5) + (121.0L/11520.0L)*pow(x_1, 4) - 1331.0L/17280.0L*pow(x_1, 3) + (14641.0L/46080.0L)*pow(x_1, 2) + (343159.0L/230400.0L)*x_1 + 347071.0L/172800.0L;
       return __pp_r607___result;
    }
    static inline O __pp_r608__(const I &x_0, const I &x_1) {
       O __pp_r608___result;
       __pp_r608___result = (19.0L/43200.0L)*pow(x_0, 6) - 1.0L/3600.0L*pow(x_0, 5)*x_1 + (127.0L/14400.0L)*pow(x_0, 5) + (1.0L/180.0L)*pow(x_0, 4)*pow(x_1, 2) + (17.0L/1440.0L)*pow(x_0, 4)*x_1 + (991.0L/11520.0L)*pow(x_0, 4) - 1.0L/270.0L*pow(x_0, 3)*pow(x_1, 3) + (11.0L/180.0L)*pow(x_0, 3)*pow(x_1, 2) + (251.0L/1440.0L)*pow(x_0, 3)*x_1 + (8443.0L/17280.0L)*pow(x_0, 3) + (1.0L/720.0L)*pow(x_0, 2)*pow(x_1, 4) - 11.0L/360.0L*pow(x_0, 2)*pow(x_1, 3) + (121.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) + (2483.0L/2880.0L)*pow(x_0, 2)*x_1 + (73999.0L/46080.0L)*pow(x_0, 2) - 1.0L/3600.0L*x_0*pow(x_1, 5) + (11.0L/1440.0L)*x_0*pow(x_1, 4) - 121.0L/1440.0L*x_0*pow(x_1, 3) + (1331.0L/2880.0L)*x_0*pow(x_1, 2) + (10687.0L/5760.0L)*x_0*x_1 + (646237.0L/230400.0L)*x_0 + (1.0L/43200.0L)*pow(x_1, 6) - 11.0L/14400.0L*pow(x_1, 5) + (121.0L/11520.0L)*pow(x_1, 4) - 1331.0L/17280.0L*pow(x_1, 3) + (14641.0L/46080.0L)*pow(x_1, 2) + (343159.0L/230400.0L)*x_1 + 347071.0L/172800.0L;
       return __pp_r608___result;
    }
    static inline O __pp_r609__(const I &x_0, const I &x_1) {
       O __pp_r609___result;
       __pp_r609___result = -7.0L/4800.0L*pow(x_0, 6) + (7.0L/3600.0L)*pow(x_0, 5)*x_1 - 367.0L/14400.0L*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) + (43.0L/1440.0L)*pow(x_0, 4)*x_1 - 2063.0L/11520.0L*pow(x_0, 4) - 1.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) - 11.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (257.0L/1440.0L)*pow(x_0, 3)*x_1 - 10987.0L/17280.0L*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) - 79.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) + (1483.0L/2880.0L)*pow(x_0, 2)*x_1 - 53363.0L/46080.0L*pow(x_0, 2) + (11.0L/1800.0L)*x_0*pow(x_1, 5) + (2.0L/45.0L)*x_0*pow(x_1, 4) + (1.0L/180.0L)*x_0*pow(x_1, 3) - 551.0L/720.0L*x_0*pow(x_1, 2) + (8177.0L/11520.0L)*x_0*x_1 - 213127.0L/230400.0L*x_0 + (1.0L/400.0L)*pow(x_1, 6) + (13.0L/450.0L)*pow(x_1, 5) + (41.0L/360.0L)*pow(x_1, 4) + (89.0L/1080.0L)*pow(x_1, 3) - 3679.0L/5760.0L*pow(x_1, 2) + (42523.0L/115200.0L)*x_1 - 398423.0L/2764800.0L;
       return __pp_r609___result;
    }
    static inline O __pp_r610__(const I &x_0, const I &x_1) {
       O __pp_r610___result;
       __pp_r610___result = -13.0L/8640.0L*pow(x_0, 6) + (1.0L/720.0L)*pow(x_0, 5)*x_1 - 77.0L/2880.0L*pow(x_0, 5) - 1.0L/144.0L*pow(x_0, 4)*pow(x_1, 2) + (5.0L/288.0L)*pow(x_0, 4)*x_1 - 445.0L/2304.0L*pow(x_0, 4) - 1.0L/108.0L*pow(x_0, 3)*pow(x_1, 3) - 1.0L/9.0L*pow(x_0, 3)*pow(x_1, 2) + (19.0L/288.0L)*pow(x_0, 3)*x_1 - 2489.0L/3456.0L*pow(x_0, 3) - 1.0L/144.0L*pow(x_0, 2)*pow(x_1, 4) - 1.0L/9.0L*pow(x_0, 2)*pow(x_1, 3) - 2.0L/3.0L*pow(x_0, 2)*pow(x_1, 2) + (5.0L/576.0L)*pow(x_0, 2)*x_1 - 13297.0L/9216.0L*pow(x_0, 2) - 1.0L/360.0L*x_0*pow(x_1, 5) - 1.0L/18.0L*x_0*pow(x_1, 4) - 4.0L/9.0L*x_0*pow(x_1, 3) - 16.0L/9.0L*x_0*pow(x_1, 2) - 989.0L/2304.0L*x_0*x_1 - 13249.0L/9216.0L*x_0 - 1.0L/2160.0L*pow(x_1, 6) - 1.0L/90.0L*pow(x_1, 5) - 1.0L/9.0L*pow(x_1, 4) - 16.0L/27.0L*pow(x_1, 3) - 16.0L/9.0L*pow(x_1, 2) - 3023.0L/4608.0L*x_1 - 292261.0L/552960.0L;
       return __pp_r610___result;
    }
    static inline O __pp_r611__(const I &x_0, const I &x_1) {
       O __pp_r611___result;
       __pp_r611___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 - 1.0L/64.0L*pow(x_0, 5) + (7.0L/96.0L)*pow(x_0, 4)*x_1 - 21.0L/256.0L*pow(x_0, 4) + (49.0L/96.0L)*pow(x_0, 3)*x_1 - 49.0L/384.0L*pow(x_0, 3) + (343.0L/192.0L)*pow(x_0, 2)*x_1 + (343.0L/1024.0L)*pow(x_0, 2) + (2401.0L/768.0L)*x_0*x_1 + (7203.0L/5120.0L)*x_0 + (16807.0L/7680.0L)*x_1 + 16807.0L/12288.0L;
       return __pp_r611___result;
    }
    static inline O __pp_r612__(const I &x_0, const I &x_1) {
       O __pp_r612___result;
       __pp_r612___result = -1.0L/960.0L*pow(x_0, 6) + (1.0L/240.0L)*pow(x_0, 5)*x_1 - 1.0L/64.0L*pow(x_0, 5) + (7.0L/96.0L)*pow(x_0, 4)*x_1 - 21.0L/256.0L*pow(x_0, 4) + (49.0L/96.0L)*pow(x_0, 3)*x_1 - 49.0L/384.0L*pow(x_0, 3) + (343.0L/192.0L)*pow(x_0, 2)*x_1 + (343.0L/1024.0L)*pow(x_0, 2) + (2401.0L/768.0L)*x_0*x_1 + (7203.0L/5120.0L)*x_0 + (16807.0L/7680.0L)*x_1 + 16807.0L/12288.0L;
       return __pp_r612___result;
    }
    static inline O __pp_r613__(const I &x_0, const I &x_1) {
       O __pp_r613___result;
       __pp_r613___result = (11.0L/43200.0L)*pow(x_0, 6) - 1.0L/900.0L*pow(x_0, 5)*x_1 + (1.0L/320.0L)*pow(x_0, 5) + (7.0L/1440.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/96.0L*pow(x_0, 4)*x_1 + (233.0L/11520.0L)*pow(x_0, 4) - 1.0L/1080.0L*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) - 47.0L/1440.0L*pow(x_0, 3)*x_1 + (607.0L/5760.0L)*pow(x_0, 3) + (13.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (53.0L/480.0L)*pow(x_0, 2)*pow(x_1, 2) - 11.0L/320.0L*pow(x_0, 2)*x_1 + (18161.0L/46080.0L)*pow(x_0, 2) - 1.0L/900.0L*x_0*pow(x_1, 5) + (5.0L/96.0L)*x_0*pow(x_1, 4) + (37.0L/1440.0L)*x_0*pow(x_1, 3) + (59.0L/960.0L)*x_0*pow(x_1, 2) - 7.0L/5760.0L*x_0*x_1 + (2489.0L/3072.0L)*x_0 + (29.0L/43200.0L)*pow(x_1, 6) - 1.0L/320.0L*pow(x_1, 5) + (863.0L/11520.0L)*pow(x_1, 4) + (301.0L/5760.0L)*pow(x_1, 3) - 3601.0L/46080.0L*pow(x_1, 2) - 29.0L/15360.0L*x_1 + 57409.0L/86400.0L;
       return __pp_r613___result;
    }
    static inline O __pp_r614__(const I &x_0, const I &x_1) {
       O __pp_r614___result;
       __pp_r614___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 - 1.0L/960.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/32.0L*pow(x_0, 4)*x_1 - 127.0L/11520.0L*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) - 227.0L/1440.0L*pow(x_0, 3)*x_1 - 113.0L/5760.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) - 37.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 131.0L/320.0L*pow(x_0, 2)*x_1 + (5201.0L/46080.0L)*pow(x_0, 2) - 1.0L/400.0L*x_0*pow(x_1, 5) + (1.0L/32.0L)*x_0*pow(x_1, 4) - 143.0L/1440.0L*x_0*pow(x_1, 3) - 301.0L/960.0L*x_0*pow(x_1, 2) - 3247.0L/5760.0L*x_0*x_1 + (7261.0L/15360.0L)*x_0 + (19.0L/43200.0L)*pow(x_1, 6) - 7.0L/960.0L*pow(x_1, 5) + (503.0L/11520.0L)*pow(x_1, 4) - 419.0L/5760.0L*pow(x_1, 3) - 16561.0L/46080.0L*pow(x_1, 2) - 5213.0L/15360.0L*x_1 + 42829.0L/86400.0L;
       return __pp_r614___result;
    }
    static inline O __pp_r615__(const I &x_0, const I &x_1) {
       O __pp_r615___result;
       __pp_r615___result = (1.0L/43200.0L)*pow(x_0, 6) - 1.0L/400.0L*pow(x_0, 5)*x_1 - 1.0L/960.0L*pow(x_0, 5) + (1.0L/720.0L)*pow(x_0, 4)*pow(x_1, 2) - 1.0L/32.0L*pow(x_0, 4)*x_1 - 127.0L/11520.0L*pow(x_0, 4) - 1.0L/180.0L*pow(x_0, 3)*pow(x_1, 3) - 227.0L/1440.0L*pow(x_0, 3)*x_1 - 113.0L/5760.0L*pow(x_0, 3) + (1.0L/180.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/24.0L*pow(x_0, 2)*pow(x_1, 3) - 37.0L/480.0L*pow(x_0, 2)*pow(x_1, 2) - 131.0L/320.0L*pow(x_0, 2)*x_1 + (5201.0L/46080.0L)*pow(x_0, 2) - 1.0L/400.0L*x_0*pow(x_1, 5) + (1.0L/32.0L)*x_0*pow(x_1, 4) - 143.0L/1440.0L*x_0*pow(x_1, 3) - 301.0L/960.0L*x_0*pow(x_1, 2) - 3247.0L/5760.0L*x_0*x_1 + (7261.0L/15360.0L)*x_0 + (19.0L/43200.0L)*pow(x_1, 6) - 7.0L/960.0L*pow(x_1, 5) + (503.0L/11520.0L)*pow(x_1, 4) - 419.0L/5760.0L*pow(x_1, 3) - 16561.0L/46080.0L*pow(x_1, 2) - 5213.0L/15360.0L*x_1 + 42829.0L/86400.0L;
       return __pp_r615___result;
    }
    static inline O __pp_r616__(const I &x_0, const I &x_1) {
       O __pp_r616___result;
       __pp_r616___result = -53.0L/43200.0L*pow(x_0, 6) + (1.0L/300.0L)*pow(x_0, 5)*x_1 - 307.0L/14400.0L*pow(x_0, 5) - 1.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (73.0L/1440.0L)*pow(x_0, 4)*x_1 - 1703.0L/11520.0L*pow(x_0, 4) + (1.0L/360.0L)*pow(x_0, 3)*pow(x_1, 3) - 7.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (437.0L/1440.0L)*pow(x_0, 3)*x_1 - 8827.0L/17280.0L*pow(x_0, 3) + (11.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (11.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 17.0L/120.0L*pow(x_0, 2)*pow(x_1, 2) + (2563.0L/2880.0L)*pow(x_0, 2)*x_1 - 40403.0L/46080.0L*pow(x_0, 2) - 1.0L/1200.0L*x_0*pow(x_1, 5) + (2.0L/45.0L)*x_0*pow(x_1, 4) + (79.0L/720.0L)*x_0*pow(x_1, 3) - 577.0L/1440.0L*x_0*pow(x_1, 2) + (14627.0L/11520.0L)*x_0*x_1 - 135427.0L/230400.0L*x_0 + (7.0L/10800.0L)*pow(x_1, 6) - 17.0L/7200.0L*pow(x_1, 5) + (371.0L/5760.0L)*pow(x_1, 4) + (1117.0L/8640.0L)*pow(x_1, 3) - 9121.0L/23040.0L*pow(x_1, 2) + (20077.0L/28800.0L)*x_1 + 65527.0L/2764800.0L;
       return __pp_r616___result;
    }
    static inline O __pp_r617__(const I &x_0, const I &x_1) {
       O __pp_r617___result;
       __pp_r617___result = -7.0L/4800.0L*pow(x_0, 6) + (7.0L/3600.0L)*pow(x_0, 5)*x_1 - 367.0L/14400.0L*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) + (43.0L/1440.0L)*pow(x_0, 4)*x_1 - 2063.0L/11520.0L*pow(x_0, 4) - 1.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) - 11.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (257.0L/1440.0L)*pow(x_0, 3)*x_1 - 10987.0L/17280.0L*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) - 79.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) + (1483.0L/2880.0L)*pow(x_0, 2)*x_1 - 53363.0L/46080.0L*pow(x_0, 2) - 1.0L/450.0L*x_0*pow(x_1, 5) + (17.0L/720.0L)*x_0*pow(x_1, 4) - 11.0L/720.0L*x_0*pow(x_1, 3) - 1117.0L/1440.0L*x_0*pow(x_1, 2) + (8147.0L/11520.0L)*x_0*x_1 - 213187.0L/230400.0L*x_0 + (1.0L/2400.0L)*pow(x_1, 6) - 47.0L/7200.0L*pow(x_1, 5) + (191.0L/5760.0L)*pow(x_1, 4) + (37.0L/8640.0L)*pow(x_1, 3) - 15601.0L/23040.0L*pow(x_1, 2) + (10357.0L/28800.0L)*x_1 - 401033.0L/2764800.0L;
       return __pp_r617___result;
    }
    static inline O __pp_r618__(const I &x_0, const I &x_1) {
       O __pp_r618___result;
       __pp_r618___result = -7.0L/4800.0L*pow(x_0, 6) + (7.0L/3600.0L)*pow(x_0, 5)*x_1 - 367.0L/14400.0L*pow(x_0, 5) - 1.0L/240.0L*pow(x_0, 4)*pow(x_1, 2) + (43.0L/1440.0L)*pow(x_0, 4)*x_1 - 2063.0L/11520.0L*pow(x_0, 4) - 1.0L/540.0L*pow(x_0, 3)*pow(x_1, 3) - 11.0L/180.0L*pow(x_0, 3)*pow(x_1, 2) + (257.0L/1440.0L)*pow(x_0, 3)*x_1 - 10987.0L/17280.0L*pow(x_0, 3) + (1.0L/240.0L)*pow(x_0, 2)*pow(x_1, 4) - 1.0L/90.0L*pow(x_0, 2)*pow(x_1, 3) - 79.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) + (1483.0L/2880.0L)*pow(x_0, 2)*x_1 - 53363.0L/46080.0L*pow(x_0, 2) - 1.0L/450.0L*x_0*pow(x_1, 5) + (17.0L/720.0L)*x_0*pow(x_1, 4) - 11.0L/720.0L*x_0*pow(x_1, 3) - 1117.0L/1440.0L*x_0*pow(x_1, 2) + (8147.0L/11520.0L)*x_0*x_1 - 213187.0L/230400.0L*x_0 + (1.0L/2400.0L)*pow(x_1, 6) - 47.0L/7200.0L*pow(x_1, 5) + (191.0L/5760.0L)*pow(x_1, 4) + (37.0L/8640.0L)*pow(x_1, 3) - 15601.0L/23040.0L*pow(x_1, 2) + (10357.0L/28800.0L)*x_1 - 401033.0L/2764800.0L;
       return __pp_r618___result;
    }
    static inline O __pp_r619__(const I &x_0, const I &x_1) {
       O __pp_r619___result;
       __pp_r619___result = -73.0L/43200.0L*pow(x_0, 6) + (11.0L/1800.0L)*pow(x_0, 5)*x_1 - 427.0L/14400.0L*pow(x_0, 5) - 11.0L/1440.0L*pow(x_0, 4)*pow(x_1, 2) + (133.0L/1440.0L)*pow(x_0, 4)*x_1 - 2423.0L/11520.0L*pow(x_0, 4) + (13.0L/1080.0L)*pow(x_0, 3)*pow(x_1, 3) - 37.0L/360.0L*pow(x_0, 3)*pow(x_1, 2) + (797.0L/1440.0L)*pow(x_0, 3)*x_1 - 13147.0L/17280.0L*pow(x_0, 3) + (1.0L/1440.0L)*pow(x_0, 2)*pow(x_1, 4) + (41.0L/360.0L)*pow(x_0, 2)*pow(x_1, 3) - 31.0L/60.0L*pow(x_0, 2)*pow(x_1, 2) + (4723.0L/2880.0L)*pow(x_0, 2)*x_1 - 66323.0L/46080.0L*pow(x_0, 2) + (7.0L/3600.0L)*x_0*pow(x_1, 5) + (1.0L/360.0L)*x_0*pow(x_1, 4) + (259.0L/720.0L)*x_0*pow(x_1, 3) - 1657.0L/1440.0L*x_0*pow(x_1, 2) + (27587.0L/11520.0L)*x_0*x_1 - 290947.0L/230400.0L*x_0 + (1.0L/5400.0L)*pow(x_1, 6) + (43.0L/7200.0L)*pow(x_1, 5) + (11.0L/5760.0L)*pow(x_1, 4) + (3277.0L/8640.0L)*pow(x_1, 3) - 22081.0L/23040.0L*pow(x_1, 2) + (39517.0L/28800.0L)*x_1 - 867593.0L/2764800.0L;
       return __pp_r619___result;
    }
    static inline O __pp_r620__(const I &x_0, const I &x_1) {
       O __pp_r620___result;
       __pp_r620___result = -83.0L/43200.0L*pow(x_0, 6) + (17.0L/3600.0L)*pow(x_0, 5)*x_1 - 487.0L/14400.0L*pow(x_0, 5) - 1.0L/90.0L*pow(x_0, 4)*pow(x_1, 2) + (103.0L/1440.0L)*pow(x_0, 4)*x_1 - 2783.0L/11520.0L*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) - 13.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (617.0L/1440.0L)*pow(x_0, 3)*x_1 - 15307.0L/17280.0L*pow(x_0, 3) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 4) + (13.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 169.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) + (3643.0L/2880.0L)*pow(x_0, 2)*x_1 - 79283.0L/46080.0L*pow(x_0, 2) + (1.0L/1800.0L)*x_0*pow(x_1, 5) - 13.0L/720.0L*x_0*pow(x_1, 4) + (169.0L/720.0L)*x_0*pow(x_1, 3) - 2197.0L/1440.0L*x_0*pow(x_1, 2) + (21107.0L/11520.0L)*x_0*x_1 - 368707.0L/230400.0L*x_0 - 1.0L/21600.0L*pow(x_1, 6) + (13.0L/7200.0L)*pow(x_1, 5) - 169.0L/5760.0L*pow(x_1, 4) + (2197.0L/8640.0L)*pow(x_1, 3) - 28561.0L/23040.0L*pow(x_1, 2) + (29797.0L/28800.0L)*x_1 - 1334153.0L/2764800.0L;
       return __pp_r620___result;
    }
    static inline O __pp_r621__(const I &x_0, const I &x_1) {
       O __pp_r621___result;
       __pp_r621___result = -83.0L/43200.0L*pow(x_0, 6) + (17.0L/3600.0L)*pow(x_0, 5)*x_1 - 487.0L/14400.0L*pow(x_0, 5) - 1.0L/90.0L*pow(x_0, 4)*pow(x_1, 2) + (103.0L/1440.0L)*pow(x_0, 4)*x_1 - 2783.0L/11520.0L*pow(x_0, 4) + (1.0L/135.0L)*pow(x_0, 3)*pow(x_1, 3) - 13.0L/90.0L*pow(x_0, 3)*pow(x_1, 2) + (617.0L/1440.0L)*pow(x_0, 3)*x_1 - 15307.0L/17280.0L*pow(x_0, 3) - 1.0L/360.0L*pow(x_0, 2)*pow(x_1, 4) + (13.0L/180.0L)*pow(x_0, 2)*pow(x_1, 3) - 169.0L/240.0L*pow(x_0, 2)*pow(x_1, 2) + (3643.0L/2880.0L)*pow(x_0, 2)*x_1 - 79283.0L/46080.0L*pow(x_0, 2) + (1.0L/1800.0L)*x_0*pow(x_1, 5) - 13.0L/720.0L*x_0*pow(x_1, 4) + (169.0L/720.0L)*x_0*pow(x_1, 3) - 2197.0L/1440.0L*x_0*pow(x_1, 2) + (21107.0L/11520.0L)*x_0*x_1 - 368707.0L/230400.0L*x_0 - 1.0L/21600.0L*pow(x_1, 6) + (13.0L/7200.0L)*pow(x_1, 5) - 169.0L/5760.0L*pow(x_1, 4) + (2197.0L/8640.0L)*pow(x_1, 3) - 28561.0L/23040.0L*pow(x_1, 2) + (29797.0L/28800.0L)*x_1 - 1334153.0L/2764800.0L;
       return __pp_r621___result;
    }
    static inline O __pp_r622__(const I &x_0, const I &x_1) {
       O __pp_r622___result;
       __pp_r622___result = (11.0L/8640.0L)*pow(x_0, 6) - 1.0L/360.0L*pow(x_0, 5)*x_1 + (9.0L/320.0L)*pow(x_0, 5) + (1.0L/288.0L)*pow(x_0, 4)*pow(x_1, 2) - 5.0L/96.0L*pow(x_0, 4)*x_1 + (199.0L/768.0L)*pow(x_0, 4) + (1.0L/216.0L)*pow(x_0, 3)*pow(x_1, 3) + (1.0L/24.0L)*pow(x_0, 3)*pow(x_1, 2) - 37.0L/96.0L*pow(x_0, 3)*x_1 + (163.0L/128.0L)*pow(x_0, 3) + (1.0L/288.0L)*pow(x_0, 2)*pow(x_1, 4) + (1.0L/24.0L)*pow(x_0, 2)*pow(x_1, 3) + (3.0L/16.0L)*pow(x_0, 2)*pow(x_1, 2) - 271.0L/192.0L*pow(x_0, 2)*x_1 + (10811.0L/3072.0L)*pow(x_0, 2) + (1.0L/720.0L)*x_0*pow(x_1, 5) + (1.0L/48.0L)*x_0*pow(x_1, 4) + (1.0L/8.0L)*x_0*pow(x_1, 3) + (3.0L/8.0L)*x_0*pow(x_1, 2) - 1969.0L/768.0L*x_0*x_1 + (15923.0L/3072.0L)*x_0 + (1.0L/4320.0L)*pow(x_1, 6) + (1.0L/240.0L)*pow(x_1, 5) + (1.0L/32.0L)*pow(x_1, 4) + (1.0L/8.0L)*pow(x_1, 3) + (9.0L/32.0L)*pow(x_1, 2) - 2843.0L/1536.0L*x_1 + 39049.0L/12288.0L;
       return __pp_r622___result;
    }
    static inline O __pp_r623__(const I &x_0, const I &x_1) {
       O __pp_r623___result;
       __pp_r623___result = (1.0L/960.0L)*pow(x_0, 6) - 1.0L/240.0L*pow(x_0, 5)*x_1 + (23.0L/960.0L)*pow(x_0, 5) - 7.0L/96.0L*pow(x_0, 4)*x_1 + (175.0L/768.0L)*pow(x_0, 4) - 49.0L/96.0L*pow(x_0, 3)*x_1 + (147.0L/128.0L)*pow(x_0, 3) - 343.0L/192.0L*pow(x_0, 2)*x_1 + (9947.0L/3072.0L)*pow(x_0, 2) - 2401.0L/768.0L*x_0*x_1 + (74431.0L/15360.0L)*x_0 - 16807.0L/7680.0L*x_1 + 184877.0L/61440.0L;
       return __pp_r623___result;
    }
    static inline O __pp_r624__(const I &x_0, const I &x_1) {
       O __pp_r624___result;
       __pp_r624___result = (1.0L/960.0L)*pow(x_0, 6) - 1.0L/240.0L*pow(x_0, 5)*x_1 + (23.0L/960.0L)*pow(x_0, 5) - 7.0L/96.0L*pow(x_0, 4)*x_1 + (175.0L/768.0L)*pow(x_0, 4) - 49.0L/96.0L*pow(x_0, 3)*x_1 + (147.0L/128.0L)*pow(x_0, 3) - 343.0L/192.0L*pow(x_0, 2)*x_1 + (9947.0L/3072.0L)*pow(x_0, 2) - 2401.0L/768.0L*x_0*x_1 + (74431.0L/15360.0L)*x_0 - 16807.0L/7680.0L*x_1 + 184877.0L/61440.0L;
       return __pp_r624___result;
    }


    static O box_spline(const I &x_0, const I &x_1) {
            if( x_0*-0.707106781187+x_1*-0.707106781187 < 0.0) {
                if( x_0*-0.707106781187+x_1*0.707106781187 < 0.0) {
                    if( x_0*-0.707106781187+x_1*0.707106781187 < -1.41421356237) {
                        if( x_0*-1.0 < -2.5) {
                            if( x_0*-0.4472135955+x_1*-0.894427191 < -1.11803398875) {
                                if( x_1*1.0 < 0.5) {
                                    if( x_0*-0.707106781187+x_1*0.707106781187 < -2.12132034356) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -2.90688837075) {
                                            if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.12132034356) {
                                                if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                                    if(x_0*-1.0 <= -3.5) return 0; return __pp_r1__(x_0, x_1);
                                                } else {
                                                    if(x_0*-1.0 <= -3.5) return 0; return __pp_r2__(x_0, x_1);
                                                }
                                            } else {
                                             return __pp_r3__(x_0, x_1);
                                         }
                                     } else {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.12132034356) {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                             return __pp_r4__(x_0, x_1);
                                         } else {
                                             return __pp_r5__(x_0, x_1);
                                         }
                                     } else {
                                        return __pp_r6__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < -2.45967477525) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.12132034356) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                            return __pp_r7__(x_0, x_1);
                                        } else {
                                            return __pp_r8__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r9__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.12132034356) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                            return __pp_r10__(x_0, x_1);
                                        } else {
                                            return __pp_r11__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r12__(x_0, x_1);
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -2.45967477525) {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.82842712475) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                        if(x_0*-1.0 <= -3.5) return 0; return __pp_r13__(x_0, x_1);
                                    } else {
                                        if(x_0*-1.0 <= -3.5) return 0; return __pp_r14__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.01246117975) {
                                        return __pp_r15__(x_0, x_1);
                                    } else {
                                        return __pp_r16__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.82842712475) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                        return __pp_r17__(x_0, x_1);
                                    } else {
                                        return __pp_r18__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.01246117975) {
                                        return __pp_r19__(x_0, x_1);
                                    } else {
                                        return __pp_r20__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    } else {
                        if( x_0*-0.707106781187+x_1*-0.707106781187 < -1.41421356237) {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -2.90688837075) {
                                if( x_0*-0.707106781187+x_1*0.707106781187 < -2.82842712475) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                        if(x_0*-0.894427191+x_1*0.4472135955 <= -3.35410196625) return 0; return __pp_r21__(x_0, x_1);
                                    } else {
                                        if(x_0*-0.894427191+x_1*0.4472135955 <= -3.35410196625) return 0; return __pp_r22__(x_0, x_1);
                                    }
                                } else {
                                    if( x_1*1.0 < -0.5) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                            return __pp_r23__(x_0, x_1);
                                        } else {
                                            return __pp_r24__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r25__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_1*1.0 < -0.5) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                        return __pp_r26__(x_0, x_1);
                                    } else {
                                        return __pp_r27__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -2.45967477525) {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < -2.12132034356) {
                                            return __pp_r28__(x_0, x_1);
                                        } else {
                                            return __pp_r29__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r30__(x_0, x_1);
                                    }
                                }
                            }
                        } else {
                            if( x_1*1.0 < -1.5) {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -0.707106781187) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                        if(x_0*-0.894427191+x_1*0.4472135955 <= -3.35410196625) return 0; return __pp_r31__(x_0, x_1);
                                    } else {
                                        if(x_0*-0.894427191+x_1*0.4472135955 <= -3.35410196625) return 0; return __pp_r32__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                        if(x_0*-0.894427191+x_1*0.4472135955 <= -3.35410196625) return 0; return __pp_r33__(x_0, x_1);
                                    } else {
                                        if(x_0*-0.894427191+x_1*0.4472135955 <= -3.35410196625) return 0; return __pp_r34__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -2.90688837075) {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < -2.82842712475) {
                                            if(x_0*-0.894427191+x_1*0.4472135955 <= -3.35410196625) return 0; return __pp_r35__(x_0, x_1);
                                        } else {
                                            return __pp_r36__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r37__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -2.90688837075) {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < -2.82842712475) {
                                            if(x_0*-0.894427191+x_1*0.4472135955 <= -3.35410196625) return 0; return __pp_r38__(x_0, x_1);
                                        } else {
                                            return __pp_r39__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r40__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -0.707106781187) {
                        if( x_1*1.0 < -0.5) {
                            if( x_0*-0.707106781187+x_1*0.707106781187 < -2.12132034356) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < -2.45967477525) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                        return __pp_r41__(x_0, x_1);
                                    } else {
                                        return __pp_r42__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                        return __pp_r43__(x_0, x_1);
                                    } else {
                                        return __pp_r44__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < -2.01246117975) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                        return __pp_r45__(x_0, x_1);
                                    } else {
                                        return __pp_r46__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                        return __pp_r47__(x_0, x_1);
                                    } else {
                                        return __pp_r48__(x_0, x_1);
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -2.01246117975) {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -1.41421356237) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.11803398875) {
                                        return __pp_r49__(x_0, x_1);
                                    } else {
                                        return __pp_r50__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                        return __pp_r51__(x_0, x_1);
                                    } else {
                                        return __pp_r52__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -1.41421356237) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.11803398875) {
                                        return __pp_r53__(x_0, x_1);
                                    } else {
                                        return __pp_r54__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                        return __pp_r55__(x_0, x_1);
                                    } else {
                                        return __pp_r56__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    } else {
                        if( x_1*1.0 < -1.5) {
                            if( x_0*-0.707106781187+x_1*0.707106781187 < -2.82842712475) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < -2.90688837075) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                        return __pp_r57__(x_0, x_1);
                                    } else {
                                        return __pp_r58__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                        return __pp_r59__(x_0, x_1);
                                    } else {
                                        return __pp_r60__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < -2.45967477525) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                        return __pp_r61__(x_0, x_1);
                                    } else {
                                        return __pp_r62__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                        return __pp_r63__(x_0, x_1);
                                    } else {
                                        return __pp_r64__(x_0, x_1);
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -2.01246117975) {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -2.45967477525) {
                                        return __pp_r65__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < -2.12132034356) {
                                            return __pp_r66__(x_0, x_1);
                                        } else {
                                            return __pp_r67__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -2.45967477525) {
                                        return __pp_r68__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < -2.12132034356) {
                                            return __pp_r69__(x_0, x_1);
                                        } else {
                                            return __pp_r70__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -1.56524758425) {
                                        if( x_0*-1.0 < -1.5) {
                                            return __pp_r71__(x_0, x_1);
                                        } else {
                                            return __pp_r72__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r73__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -1.56524758425) {
                                        if( x_0*-1.0 < -1.5) {
                                            return __pp_r74__(x_0, x_1);
                                        } else {
                                            return __pp_r75__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r76__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if( x_0*-0.707106781187+x_1*0.707106781187 < -0.707106781187) {
                    if( x_1*1.0 < 0.5) {
                        if( x_0*-1.0 < -1.5) {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -1.56524758425) {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -1.41421356237) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.11803398875) {
                                        return __pp_r77__(x_0, x_1);
                                    } else {
                                        return __pp_r78__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                        return __pp_r79__(x_0, x_1);
                                    } else {
                                        return __pp_r80__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -1.41421356237) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.11803398875) {
                                        return __pp_r81__(x_0, x_1);
                                    } else {
                                        return __pp_r82__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                        return __pp_r83__(x_0, x_1);
                                    } else {
                                        return __pp_r84__(x_0, x_1);
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -1.11803398875) {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -0.707106781187) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                            return __pp_r85__(x_0, x_1);
                                        } else {
                                            return __pp_r86__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r87__(x_0, x_1);
                                    }
                                } else {
                                    if( x_1*1.0 < -0.5) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                            return __pp_r88__(x_0, x_1);
                                        } else {
                                            return __pp_r89__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r90__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -0.707106781187) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                            return __pp_r91__(x_0, x_1);
                                        } else {
                                            return __pp_r92__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r93__(x_0, x_1);
                                    }
                                } else {
                                    if( x_1*1.0 < -0.5) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                            return __pp_r94__(x_0, x_1);
                                        } else {
                                            return __pp_r95__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r96__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    } else {
                        if( x_0*-1.0 < -2.5) {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -2.01246117975) {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                    if( x_1*1.0 < 1.5) {
                                        return __pp_r97__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -2.90688837075) {
                                            if(x_0*-0.707106781187+x_1*-0.707106781187 <= -3.53553390593) return 0; return __pp_r98__(x_0, x_1);
                                        } else {
                                            return __pp_r99__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.82842712475) {
                                        return __pp_r100__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -2.01246117975) {
                                            return __pp_r101__(x_0, x_1);
                                        } else {
                                            return __pp_r102__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                    if( x_1*1.0 < 1.5) {
                                        return __pp_r103__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -2.90688837075) {
                                            if(x_0*-0.707106781187+x_1*-0.707106781187 <= -3.53553390593) return 0; return __pp_r104__(x_0, x_1);
                                        } else {
                                            return __pp_r105__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.82842712475) {
                                        return __pp_r106__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -2.01246117975) {
                                            return __pp_r107__(x_0, x_1);
                                        } else {
                                            return __pp_r108__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -1.56524758425) {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.12132034356) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.01246117975) {
                                        return __pp_r109__(x_0, x_1);
                                    } else {
                                        return __pp_r110__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r111__(x_0, x_1);
                                    } else {
                                        return __pp_r112__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.12132034356) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.01246117975) {
                                        return __pp_r113__(x_0, x_1);
                                    } else {
                                        return __pp_r114__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r115__(x_0, x_1);
                                    } else {
                                        return __pp_r116__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if( x_0*-1.0 < -1.5) {
                        if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                            if( x_0*-1.0 < -2.5) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < -1.56524758425) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.90688837075) {
                                        if(x_0*-0.707106781187+x_1*-0.707106781187 <= -3.53553390593) return 0; return __pp_r117__(x_0, x_1);
                                    } else {
                                        return __pp_r118__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.90688837075) {
                                        if(x_0*-0.707106781187+x_1*-0.707106781187 <= -3.53553390593) return 0; return __pp_r119__(x_0, x_1);
                                    } else {
                                        return __pp_r120__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < -1.11803398875) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.82842712475) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -2.90688837075) {
                                            return __pp_r121__(x_0, x_1);
                                        } else {
                                            return __pp_r122__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r123__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.82842712475) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -2.90688837075) {
                                            return __pp_r124__(x_0, x_1);
                                        } else {
                                            return __pp_r125__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r126__(x_0, x_1);
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -1.11803398875) {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.12132034356) {
                                    if( x_1*1.0 < 1.5) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -2.01246117975) {
                                            return __pp_r127__(x_0, x_1);
                                        } else {
                                            return __pp_r128__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r129__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r130__(x_0, x_1);
                                    } else {
                                        return __pp_r131__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.12132034356) {
                                    if( x_1*1.0 < 1.5) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -2.01246117975) {
                                            return __pp_r132__(x_0, x_1);
                                        } else {
                                            return __pp_r133__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r134__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r135__(x_0, x_1);
                                    } else {
                                        return __pp_r136__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    } else {
                        if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -0.67082039325) {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -1.41421356237) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r137__(x_0, x_1);
                                    } else {
                                        return __pp_r138__(x_0, x_1);
                                    }
                                } else {
                                    if( x_1*1.0 < 0.5) {
                                        return __pp_r139__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -1.11803398875) {
                                            return __pp_r140__(x_0, x_1);
                                        } else {
                                            return __pp_r141__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -1.41421356237) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r142__(x_0, x_1);
                                    } else {
                                        return __pp_r143__(x_0, x_1);
                                    }
                                } else {
                                    if( x_1*1.0 < 0.5) {
                                        return __pp_r144__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -1.11803398875) {
                                            return __pp_r145__(x_0, x_1);
                                        } else {
                                            return __pp_r146__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        } else {
                            if( x_0*-1.0 < -0.5) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < -0.67082039325) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -0.707106781187) {
                                        return __pp_r147__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                            return __pp_r148__(x_0, x_1);
                                        } else {
                                            return __pp_r149__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -0.707106781187) {
                                        return __pp_r150__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                            return __pp_r151__(x_0, x_1);
                                        } else {
                                            return __pp_r152__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                        return __pp_r153__(x_0, x_1);
                                    } else {
                                        return __pp_r154__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                        return __pp_r155__(x_0, x_1);
                                    } else {
                                        return __pp_r156__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if( x_0*-0.707106781187+x_1*-0.707106781187 < -1.41421356237) {
                if( x_1*1.0 < 2.5) {
                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.12132034356) {
                        if( x_0*-1.0 < -1.5) {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -0.67082039325) {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.82842712475) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.90688837075) {
                                        return __pp_r157__(x_0, x_1);
                                    } else {
                                        return __pp_r158__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                        return __pp_r159__(x_0, x_1);
                                    } else {
                                        return __pp_r160__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.82842712475) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.90688837075) {
                                        return __pp_r161__(x_0, x_1);
                                    } else {
                                        return __pp_r162__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                        return __pp_r163__(x_0, x_1);
                                    } else {
                                        return __pp_r164__(x_0, x_1);
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.707106781187+x_1*0.707106781187 < 0.707106781187) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                        return __pp_r165__(x_0, x_1);
                                    } else {
                                        return __pp_r166__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                        return __pp_r167__(x_0, x_1);
                                    } else {
                                        return __pp_r168__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                        return __pp_r169__(x_0, x_1);
                                    } else {
                                        return __pp_r170__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                        return __pp_r171__(x_0, x_1);
                                    } else {
                                        return __pp_r172__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    } else {
                        if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                            if( x_1*1.0 < 1.5) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r173__(x_0, x_1);
                                    } else {
                                        return __pp_r174__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r175__(x_0, x_1);
                                    } else {
                                        return __pp_r176__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < -2.01246117975) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                        return __pp_r177__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 0.707106781187) {
                                            return __pp_r178__(x_0, x_1);
                                        } else {
                                            return __pp_r179__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                        return __pp_r180__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 0.707106781187) {
                                            return __pp_r181__(x_0, x_1);
                                        } else {
                                            return __pp_r182__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.4472135955+x_1*-0.894427191 < -2.01246117975) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 0.67082039325) {
                                    if( x_0*-1.0 < -0.5) {
                                        return __pp_r183__(x_0, x_1);
                                    } else {
                                        return __pp_r184__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 1.11803398875) {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 1.41421356237) {
                                            return __pp_r185__(x_0, x_1);
                                        } else {
                                            return __pp_r186__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r187__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 0.67082039325) {
                                    if( x_0*-1.0 < -0.5) {
                                        return __pp_r188__(x_0, x_1);
                                    } else {
                                        return __pp_r189__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 1.11803398875) {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 1.41421356237) {
                                            return __pp_r190__(x_0, x_1);
                                        } else {
                                            return __pp_r191__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r192__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if( x_0*-0.894427191+x_1*0.4472135955 < 1.11803398875) {
                        if( x_0*-0.707106781187+x_1*0.707106781187 < 1.41421356237) {
                            if( x_0*-1.0 < -1.5) {
                                if( x_0*-0.707106781187+x_1*0.707106781187 < 0.707106781187) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.67082039325) {
                                        if(x_0*-0.4472135955+x_1*-0.894427191 <= -3.35410196625) return 0; return __pp_r193__(x_0, x_1);
                                    } else {
                                        if(x_0*-0.4472135955+x_1*-0.894427191 <= -3.35410196625) return 0; return __pp_r194__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                        if(x_0*-0.4472135955+x_1*-0.894427191 <= -3.35410196625) return 0; return __pp_r195__(x_0, x_1);
                                    } else {
                                        if(x_0*-0.4472135955+x_1*-0.894427191 <= -3.35410196625) return 0; return __pp_r196__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.82842712475) {
                                        if(x_0*-0.4472135955+x_1*-0.894427191 <= -3.35410196625) return 0; return __pp_r197__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -2.90688837075) {
                                            return __pp_r198__(x_0, x_1);
                                        } else {
                                            return __pp_r199__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.82842712475) {
                                        if(x_0*-0.4472135955+x_1*-0.894427191 <= -3.35410196625) return 0; return __pp_r200__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -2.90688837075) {
                                            return __pp_r201__(x_0, x_1);
                                        } else {
                                            return __pp_r202__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.4472135955+x_1*-0.894427191 < -2.90688837075) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 0.67082039325) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.82842712475) {
                                        if(x_0*-0.4472135955+x_1*-0.894427191 <= -3.35410196625) return 0; return __pp_r203__(x_0, x_1);
                                    } else {
                                        return __pp_r204__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-1.0 < -0.5) {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.82842712475) {
                                            if(x_0*-0.4472135955+x_1*-0.894427191 <= -3.35410196625) return 0; return __pp_r205__(x_0, x_1);
                                        } else {
                                            return __pp_r206__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r207__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-1.0 < -0.5) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 0.67082039325) {
                                        return __pp_r208__(x_0, x_1);
                                    } else {
                                        return __pp_r209__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.12132034356) {
                                        return __pp_r210__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                            return __pp_r211__(x_0, x_1);
                                        } else {
                                            return __pp_r212__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        if( x_0*-0.894427191+x_1*0.4472135955 < 1.56524758425) {
                            if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.12132034356) {
                                if( x_0*-0.707106781187+x_1*0.707106781187 < 2.12132034356) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.90688837075) {
                                        return __pp_r213__(x_0, x_1);
                                    } else {
                                        return __pp_r214__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.90688837075) {
                                        if(x_1*1.0 >= 3.5) return 0; return __pp_r215__(x_0, x_1);
                                    } else {
                                        return __pp_r216__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.707106781187+x_1*0.707106781187 < 2.12132034356) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                        return __pp_r217__(x_0, x_1);
                                    } else {
                                        return __pp_r218__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                        return __pp_r219__(x_0, x_1);
                                    } else {
                                        return __pp_r220__(x_0, x_1);
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.894427191+x_1*0.4472135955 < 2.01246117975) {
                                if( x_0*-1.0 < 0.5) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < -2.12132034356) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -2.90688837075) {
                                            if(x_1*1.0 >= 3.5) return 0; return __pp_r221__(x_0, x_1);
                                        } else {
                                            return __pp_r222__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                            return __pp_r223__(x_0, x_1);
                                        } else {
                                            return __pp_r224__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                        return __pp_r225__(x_0, x_1);
                                    } else {
                                        return __pp_r226__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < -2.45967477525) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 2.45967477525) {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 2.82842712475) {
                                            return __pp_r227__(x_0, x_1);
                                        } else {
                                            if(x_1*1.0 >= 3.5) return 0; return __pp_r228__(x_0, x_1);
                                        }
                                    } else {
                                        if(x_1*1.0 >= 3.5) return 0; return __pp_r229__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 2.45967477525) {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 2.82842712475) {
                                            return __pp_r230__(x_0, x_1);
                                        } else {
                                            return __pp_r231__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r232__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if( x_0*-0.707106781187+x_1*-0.707106781187 < -0.707106781187) {
                    if( x_0*-1.0 < 0.5) {
                        if( x_1*1.0 < 1.5) {
                            if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < -1.11803398875) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                        return __pp_r233__(x_0, x_1);
                                    } else {
                                        if( x_0*-1.0 < -0.5) {
                                            return __pp_r234__(x_0, x_1);
                                        } else {
                                            return __pp_r235__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                        return __pp_r236__(x_0, x_1);
                                    } else {
                                        if( x_0*-1.0 < -0.5) {
                                            return __pp_r237__(x_0, x_1);
                                        } else {
                                            return __pp_r238__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < -1.11803398875) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 0.67082039325) {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 0.707106781187) {
                                            return __pp_r239__(x_0, x_1);
                                        } else {
                                            return __pp_r240__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r241__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 0.67082039325) {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 0.707106781187) {
                                            return __pp_r242__(x_0, x_1);
                                        } else {
                                            return __pp_r243__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r244__(x_0, x_1);
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.707106781187+x_1*0.707106781187 < 1.41421356237) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 0.67082039325) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r245__(x_0, x_1);
                                    } else {
                                        return __pp_r246__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r247__(x_0, x_1);
                                    } else {
                                        return __pp_r248__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 1.11803398875) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r249__(x_0, x_1);
                                    } else {
                                        return __pp_r250__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r251__(x_0, x_1);
                                    } else {
                                        return __pp_r252__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    } else {
                        if( x_1*1.0 < 2.5) {
                            if( x_0*-0.707106781187+x_1*0.707106781187 < 2.12132034356) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 1.56524758425) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r253__(x_0, x_1);
                                    } else {
                                        return __pp_r254__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r255__(x_0, x_1);
                                    } else {
                                        return __pp_r256__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 2.01246117975) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r257__(x_0, x_1);
                                    } else {
                                        return __pp_r258__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r259__(x_0, x_1);
                                    } else {
                                        return __pp_r260__(x_0, x_1);
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.894427191+x_1*0.4472135955 < 2.45967477525) {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < -2.01246117975) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 2.01246117975) {
                                        return __pp_r261__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 2.82842712475) {
                                            return __pp_r262__(x_0, x_1);
                                        } else {
                                            return __pp_r263__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 2.01246117975) {
                                        return __pp_r264__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 2.82842712475) {
                                            return __pp_r265__(x_0, x_1);
                                        } else {
                                            return __pp_r266__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < -2.01246117975) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 2.90688837075) {
                                        if( x_0*-1.0 < 1.5) {
                                            return __pp_r267__(x_0, x_1);
                                        } else {
                                            return __pp_r268__(x_0, x_1);
                                        }
                                    } else {
                                        if(x_0*-0.707106781187+x_1*0.707106781187 >= 3.53553390593) return 0; return __pp_r269__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 2.90688837075) {
                                        if( x_0*-1.0 < 1.5) {
                                            return __pp_r270__(x_0, x_1);
                                        } else {
                                            return __pp_r271__(x_0, x_1);
                                        }
                                    } else {
                                        if(x_0*-0.707106781187+x_1*0.707106781187 >= 3.53553390593) return 0; return __pp_r272__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if( x_1*1.0 < 1.5) {
                        if( x_0*-0.894427191+x_1*0.4472135955 < 0.67082039325) {
                            if( x_1*1.0 < 0.5) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                        return __pp_r273__(x_0, x_1);
                                    } else {
                                        return __pp_r274__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                        return __pp_r275__(x_0, x_1);
                                    } else {
                                        return __pp_r276__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                                        return __pp_r277__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 0.707106781187) {
                                            return __pp_r278__(x_0, x_1);
                                        } else {
                                            return __pp_r279__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                                        return __pp_r280__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 0.707106781187) {
                                            return __pp_r281__(x_0, x_1);
                                        } else {
                                            return __pp_r282__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 1.11803398875) {
                                    if( x_0*-1.0 < 0.5) {
                                        return __pp_r283__(x_0, x_1);
                                    } else {
                                        return __pp_r284__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 1.56524758425) {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 1.41421356237) {
                                            return __pp_r285__(x_0, x_1);
                                        } else {
                                            return __pp_r286__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r287__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 1.11803398875) {
                                    if( x_0*-1.0 < 0.5) {
                                        return __pp_r288__(x_0, x_1);
                                    } else {
                                        return __pp_r289__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 1.56524758425) {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 1.41421356237) {
                                            return __pp_r290__(x_0, x_1);
                                        } else {
                                            return __pp_r291__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r292__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    } else {
                        if( x_0*-0.894427191+x_1*0.4472135955 < 2.45967477525) {
                            if( x_0*-0.4472135955+x_1*-0.894427191 < -1.11803398875) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 2.01246117975) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 1.56524758425) {
                                        return __pp_r293__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 2.12132034356) {
                                            return __pp_r294__(x_0, x_1);
                                        } else {
                                            return __pp_r295__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-1.0 < 1.5) {
                                        return __pp_r296__(x_0, x_1);
                                    } else {
                                        return __pp_r297__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 2.01246117975) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 1.56524758425) {
                                        return __pp_r298__(x_0, x_1);
                                    } else {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 2.12132034356) {
                                            return __pp_r299__(x_0, x_1);
                                        } else {
                                            return __pp_r300__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-1.0 < 1.5) {
                                        return __pp_r301__(x_0, x_1);
                                    } else {
                                        return __pp_r302__(x_0, x_1);
                                    }
                                }
                            }
                        } else {
                            if( x_1*1.0 < 2.5) {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < -1.11803398875) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 2.90688837075) {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 2.82842712475) {
                                            return __pp_r303__(x_0, x_1);
                                        } else {
                                            return __pp_r304__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r305__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 2.90688837075) {
                                        if( x_0*-0.707106781187+x_1*0.707106781187 < 2.82842712475) {
                                            return __pp_r306__(x_0, x_1);
                                        } else {
                                            return __pp_r307__(x_0, x_1);
                                        }
                                    } else {
                                        return __pp_r308__(x_0, x_1);
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 2.90688837075) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        return __pp_r309__(x_0, x_1);
                                    } else {
                                        return __pp_r310__(x_0, x_1);
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -1.56524758425) {
                                        if(x_0*-0.707106781187+x_1*0.707106781187 >= 3.53553390593) return 0; return __pp_r311__(x_0, x_1);
                                    } else {
                                        if(x_0*-0.707106781187+x_1*0.707106781187 >= 3.53553390593) return 0; return __pp_r312__(x_0, x_1);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        } else {
            if( x_0*-0.707106781187+x_1*0.707106781187 < 0.0) {
                if( x_0*-0.707106781187+x_1*-0.707106781187 < 1.41421356237) {
                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 0.707106781187) {
                        if( x_1*1.0 < -1.5) {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -2.45967477525) {
                                if( x_1*1.0 < -2.5) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -2.90688837075) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            if(x_0*-0.707106781187+x_1*0.707106781187 <= -3.53553390593) return 0; return __pp_r313__(x_0, x_1);
                                        } else {
                                            if(x_0*-0.707106781187+x_1*0.707106781187 <= -3.53553390593) return 0; return __pp_r314__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r315__(x_0, x_1);
                                        } else {
                                            return __pp_r316__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 1.11803398875) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -2.90688837075) {
                                            return __pp_r317__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -2.82842712475) {
                                                return __pp_r318__(x_0, x_1);
                                            } else {
                                                return __pp_r319__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -2.90688837075) {
                                            return __pp_r320__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -2.82842712475) {
                                                return __pp_r321__(x_0, x_1);
                                            } else {
                                                return __pp_r322__(x_0, x_1);
                                            }
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < 1.11803398875) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -2.01246117975) {
                                        if( x_0*-1.0 < -1.5) {
                                            return __pp_r323__(x_0, x_1);
                                        } else {
                                            return __pp_r324__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -1.56524758425) {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -2.12132034356) {
                                                return __pp_r325__(x_0, x_1);
                                            } else {
                                                return __pp_r326__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r327__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -2.01246117975) {
                                        if( x_0*-1.0 < -1.5) {
                                            return __pp_r328__(x_0, x_1);
                                        } else {
                                            return __pp_r329__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -1.56524758425) {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -2.12132034356) {
                                                return __pp_r330__(x_0, x_1);
                                            } else {
                                                return __pp_r331__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r332__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -0.67082039325) {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -1.11803398875) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -1.56524758425) {
                                            return __pp_r333__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -1.41421356237) {
                                                return __pp_r334__(x_0, x_1);
                                            } else {
                                                return __pp_r335__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-1.0 < -0.5) {
                                            return __pp_r336__(x_0, x_1);
                                        } else {
                                            return __pp_r337__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -1.11803398875) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -1.56524758425) {
                                            return __pp_r338__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -1.41421356237) {
                                                return __pp_r339__(x_0, x_1);
                                            } else {
                                                return __pp_r340__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-1.0 < -0.5) {
                                            return __pp_r341__(x_0, x_1);
                                        } else {
                                            return __pp_r342__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_1*1.0 < -0.5) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -0.707106781187) {
                                                return __pp_r343__(x_0, x_1);
                                            } else {
                                                return __pp_r344__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r345__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -0.707106781187) {
                                                return __pp_r346__(x_0, x_1);
                                            } else {
                                                return __pp_r347__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r348__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                            return __pp_r349__(x_0, x_1);
                                        } else {
                                            return __pp_r350__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                            return __pp_r351__(x_0, x_1);
                                        } else {
                                            return __pp_r352__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        if( x_0*-1.0 < -0.5) {
                            if( x_1*1.0 < -2.5) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < -2.45967477525) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 2.01246117975) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -2.90688837075) {
                                            if(x_0*-0.707106781187+x_1*0.707106781187 <= -3.53553390593) return 0; return __pp_r353__(x_0, x_1);
                                        } else {
                                            if( x_0*-1.0 < -1.5) {
                                                return __pp_r354__(x_0, x_1);
                                            } else {
                                                return __pp_r355__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -2.90688837075) {
                                            if(x_0*-0.707106781187+x_1*0.707106781187 <= -3.53553390593) return 0; return __pp_r356__(x_0, x_1);
                                        } else {
                                            if( x_0*-1.0 < -1.5) {
                                                return __pp_r357__(x_0, x_1);
                                            } else {
                                                return __pp_r358__(x_0, x_1);
                                            }
                                        }
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 2.01246117975) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -2.01246117975) {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -2.82842712475) {
                                                return __pp_r359__(x_0, x_1);
                                            } else {
                                                return __pp_r360__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r361__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -2.01246117975) {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -2.82842712475) {
                                                return __pp_r362__(x_0, x_1);
                                            } else {
                                                return __pp_r363__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r364__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.707106781187+x_1*0.707106781187 < -2.12132034356) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -2.01246117975) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r365__(x_0, x_1);
                                        } else {
                                            return __pp_r366__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r367__(x_0, x_1);
                                        } else {
                                            return __pp_r368__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -1.56524758425) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r369__(x_0, x_1);
                                        } else {
                                            return __pp_r370__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r371__(x_0, x_1);
                                        } else {
                                            return __pp_r372__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        } else {
                            if( x_1*1.0 < -1.5) {
                                if( x_0*-0.707106781187+x_1*0.707106781187 < -1.41421356237) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -1.11803398875) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r373__(x_0, x_1);
                                        } else {
                                            return __pp_r374__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r375__(x_0, x_1);
                                        } else {
                                            return __pp_r376__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.67082039325) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r377__(x_0, x_1);
                                        } else {
                                            return __pp_r378__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r379__(x_0, x_1);
                                        } else {
                                            return __pp_r380__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 1.11803398875) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -0.67082039325) {
                                            return __pp_r381__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -0.707106781187) {
                                                return __pp_r382__(x_0, x_1);
                                            } else {
                                                return __pp_r383__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -0.67082039325) {
                                            return __pp_r384__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -0.707106781187) {
                                                return __pp_r385__(x_0, x_1);
                                            } else {
                                                return __pp_r386__(x_0, x_1);
                                            }
                                        }
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 1.11803398875) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                                            if( x_0*-1.0 < 0.5) {
                                                return __pp_r387__(x_0, x_1);
                                            } else {
                                                return __pp_r388__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r389__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                                            if( x_0*-1.0 < 0.5) {
                                                return __pp_r390__(x_0, x_1);
                                            } else {
                                                return __pp_r391__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r392__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if( x_1*1.0 < -2.5) {
                        if( x_0*-0.894427191+x_1*0.4472135955 < -1.11803398875) {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -1.56524758425) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < -2.01246117975) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -2.45967477525) {
                                            return __pp_r393__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -2.82842712475) {
                                                return __pp_r394__(x_0, x_1);
                                            } else {
                                                return __pp_r395__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -2.45967477525) {
                                            if(x_1*1.0 <= -3.5) return 0; return __pp_r396__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -2.82842712475) {
                                                if(x_1*1.0 <= -3.5) return 0; return __pp_r397__(x_0, x_1);
                                            } else {
                                                return __pp_r398__(x_0, x_1);
                                            }
                                        }
                                    }
                                } else {
                                    if( x_0*-1.0 < -0.5) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                            return __pp_r399__(x_0, x_1);
                                        } else {
                                            return __pp_r400__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.12132034356) {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                                return __pp_r401__(x_0, x_1);
                                            } else {
                                                return __pp_r402__(x_0, x_1);
                                            }
                                        } else {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 2.90688837075) {
                                                return __pp_r403__(x_0, x_1);
                                            } else {
                                                if(x_1*1.0 <= -3.5) return 0; return __pp_r404__(x_0, x_1);
                                            }
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.12132034356) {
                                    if( x_0*-0.707106781187+x_1*0.707106781187 < -2.12132034356) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                            return __pp_r405__(x_0, x_1);
                                        } else {
                                            return __pp_r406__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                            return __pp_r407__(x_0, x_1);
                                        } else {
                                            return __pp_r408__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*0.707106781187 < -2.12132034356) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.90688837075) {
                                            return __pp_r409__(x_0, x_1);
                                        } else {
                                            if(x_1*1.0 <= -3.5) return 0; return __pp_r410__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.90688837075) {
                                            return __pp_r411__(x_0, x_1);
                                        } else {
                                            return __pp_r412__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.707106781187+x_1*0.707106781187 < -1.41421356237) {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < 2.90688837075) {
                                    if( x_0*-1.0 < 0.5) {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.12132034356) {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                                return __pp_r413__(x_0, x_1);
                                            } else {
                                                return __pp_r414__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r415__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -0.67082039325) {
                                            return __pp_r416__(x_0, x_1);
                                        } else {
                                            return __pp_r417__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.67082039325) {
                                        if( x_0*-1.0 < 0.5) {
                                            return __pp_r418__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.82842712475) {
                                                return __pp_r419__(x_0, x_1);
                                            } else {
                                                if(x_0*-0.4472135955+x_1*-0.894427191 >= 3.35410196625) return 0; return __pp_r420__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.82842712475) {
                                            return __pp_r421__(x_0, x_1);
                                        } else {
                                            if(x_0*-0.4472135955+x_1*-0.894427191 >= 3.35410196625) return 0; return __pp_r422__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-1.0 < 1.5) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.82842712475) {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 2.90688837075) {
                                                return __pp_r423__(x_0, x_1);
                                            } else {
                                                return __pp_r424__(x_0, x_1);
                                            }
                                        } else {
                                            if(x_0*-0.4472135955+x_1*-0.894427191 >= 3.35410196625) return 0; return __pp_r425__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.82842712475) {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 2.90688837075) {
                                                return __pp_r426__(x_0, x_1);
                                            } else {
                                                return __pp_r427__(x_0, x_1);
                                            }
                                        } else {
                                            if(x_0*-0.4472135955+x_1*-0.894427191 >= 3.35410196625) return 0; return __pp_r428__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*0.707106781187 < -0.707106781187) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                                            if(x_0*-0.4472135955+x_1*-0.894427191 >= 3.35410196625) return 0; return __pp_r429__(x_0, x_1);
                                        } else {
                                            if(x_0*-0.4472135955+x_1*-0.894427191 >= 3.35410196625) return 0; return __pp_r430__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < 0.67082039325) {
                                            if(x_0*-0.4472135955+x_1*-0.894427191 >= 3.35410196625) return 0; return __pp_r431__(x_0, x_1);
                                        } else {
                                            if(x_0*-0.4472135955+x_1*-0.894427191 >= 3.35410196625) return 0; return __pp_r432__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.12132034356) {
                            if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                if( x_0*-0.4472135955+x_1*-0.894427191 < 2.01246117975) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.67082039325) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -1.11803398875) {
                                            return __pp_r433__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -1.41421356237) {
                                                return __pp_r434__(x_0, x_1);
                                            } else {
                                                return __pp_r435__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-1.0 < 0.5) {
                                            return __pp_r436__(x_0, x_1);
                                        } else {
                                            return __pp_r437__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.67082039325) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < -1.11803398875) {
                                            return __pp_r438__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -1.41421356237) {
                                                return __pp_r439__(x_0, x_1);
                                            } else {
                                                return __pp_r440__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-1.0 < 0.5) {
                                            return __pp_r441__(x_0, x_1);
                                        } else {
                                            return __pp_r442__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_1*1.0 < -1.5) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 2.01246117975) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -0.707106781187) {
                                                return __pp_r443__(x_0, x_1);
                                            } else {
                                                return __pp_r444__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r445__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < -0.707106781187) {
                                                return __pp_r446__(x_0, x_1);
                                            } else {
                                                return __pp_r447__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r448__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r449__(x_0, x_1);
                                        } else {
                                            return __pp_r450__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r451__(x_0, x_1);
                                        } else {
                                            return __pp_r452__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        } else {
                            if( x_0*-1.0 < 1.5) {
                                if( x_0*-0.707106781187+x_1*0.707106781187 < -0.707106781187) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.22360679775) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                            return __pp_r453__(x_0, x_1);
                                        } else {
                                            return __pp_r454__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                            return __pp_r455__(x_0, x_1);
                                        } else {
                                            return __pp_r456__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                            return __pp_r457__(x_0, x_1);
                                        } else {
                                            return __pp_r458__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                            return __pp_r459__(x_0, x_1);
                                        } else {
                                            return __pp_r460__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 0.67082039325) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.82842712475) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                            return __pp_r461__(x_0, x_1);
                                        } else {
                                            return __pp_r462__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.90688837075) {
                                            return __pp_r463__(x_0, x_1);
                                        } else {
                                            return __pp_r464__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.82842712475) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                            return __pp_r465__(x_0, x_1);
                                        } else {
                                            return __pp_r466__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.90688837075) {
                                            return __pp_r467__(x_0, x_1);
                                        } else {
                                            return __pp_r468__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if( x_0*-0.707106781187+x_1*0.707106781187 < 1.41421356237) {
                    if( x_0*-0.707106781187+x_1*0.707106781187 < 0.707106781187) {
                        if( x_0*-1.0 < 1.5) {
                            if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                if( x_0*-1.0 < 0.5) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 0.22360679775) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                            return __pp_r469__(x_0, x_1);
                                        } else {
                                            return __pp_r470__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                            return __pp_r471__(x_0, x_1);
                                        } else {
                                            return __pp_r472__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 0.67082039325) {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 0.707106781187) {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                                return __pp_r473__(x_0, x_1);
                                            } else {
                                                return __pp_r474__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r475__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 0.707106781187) {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                                return __pp_r476__(x_0, x_1);
                                            } else {
                                                return __pp_r477__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r478__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 0.67082039325) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 1.41421356237) {
                                        if( x_1*1.0 < -0.5) {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 1.11803398875) {
                                                return __pp_r479__(x_0, x_1);
                                            } else {
                                                return __pp_r480__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r481__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r482__(x_0, x_1);
                                        } else {
                                            return __pp_r483__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 1.41421356237) {
                                        if( x_1*1.0 < -0.5) {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 1.11803398875) {
                                                return __pp_r484__(x_0, x_1);
                                            } else {
                                                return __pp_r485__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r486__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r487__(x_0, x_1);
                                        } else {
                                            return __pp_r488__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        } else {
                            if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 1.11803398875) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.12132034356) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r489__(x_0, x_1);
                                        } else {
                                            return __pp_r490__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_1*1.0 < -1.5) {
                                            return __pp_r491__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 2.01246117975) {
                                                return __pp_r492__(x_0, x_1);
                                            } else {
                                                return __pp_r493__(x_0, x_1);
                                            }
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.12132034356) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r494__(x_0, x_1);
                                        } else {
                                            return __pp_r495__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_1*1.0 < -1.5) {
                                            return __pp_r496__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 2.01246117975) {
                                                return __pp_r497__(x_0, x_1);
                                            } else {
                                                return __pp_r498__(x_0, x_1);
                                            }
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-1.0 < 2.5) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 1.11803398875) {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.82842712475) {
                                            return __pp_r499__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 2.90688837075) {
                                                return __pp_r500__(x_0, x_1);
                                            } else {
                                                return __pp_r501__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.82842712475) {
                                            return __pp_r502__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 2.90688837075) {
                                                return __pp_r503__(x_0, x_1);
                                            } else {
                                                return __pp_r504__(x_0, x_1);
                                            }
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 1.56524758425) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.90688837075) {
                                            return __pp_r505__(x_0, x_1);
                                        } else {
                                            if(x_0*-0.707106781187+x_1*-0.707106781187 >= 3.53553390593) return 0; return __pp_r506__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.90688837075) {
                                            return __pp_r507__(x_0, x_1);
                                        } else {
                                            if(x_0*-0.707106781187+x_1*-0.707106781187 >= 3.53553390593) return 0; return __pp_r508__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        if( x_1*1.0 < -0.5) {
                            if( x_0*-1.0 < 2.5) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 1.56524758425) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.12132034356) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r509__(x_0, x_1);
                                        } else {
                                            return __pp_r510__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.01246117975) {
                                            return __pp_r511__(x_0, x_1);
                                        } else {
                                            return __pp_r512__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.12132034356) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                            return __pp_r513__(x_0, x_1);
                                        } else {
                                            return __pp_r514__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.01246117975) {
                                            return __pp_r515__(x_0, x_1);
                                        } else {
                                            return __pp_r516__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 2.01246117975) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.82842712475) {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 2.01246117975) {
                                                return __pp_r517__(x_0, x_1);
                                            } else {
                                                return __pp_r518__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r519__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_1*1.0 < -1.5) {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 2.90688837075) {
                                                return __pp_r520__(x_0, x_1);
                                            } else {
                                                if(x_0*-0.707106781187+x_1*-0.707106781187 >= 3.53553390593) return 0; return __pp_r521__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r522__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.82842712475) {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 2.01246117975) {
                                                return __pp_r523__(x_0, x_1);
                                            } else {
                                                return __pp_r524__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r525__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_1*1.0 < -1.5) {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 2.90688837075) {
                                                return __pp_r526__(x_0, x_1);
                                            } else {
                                                if(x_0*-0.707106781187+x_1*-0.707106781187 >= 3.53553390593) return 0; return __pp_r527__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r528__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        } else {
                            if( x_0*-1.0 < 1.5) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 1.11803398875) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                        if( x_1*1.0 < 0.5) {
                                            return __pp_r529__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                                return __pp_r530__(x_0, x_1);
                                            } else {
                                                return __pp_r531__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 0.707106781187) {
                                            return __pp_r532__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                                return __pp_r533__(x_0, x_1);
                                            } else {
                                                return __pp_r534__(x_0, x_1);
                                            }
                                        }
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                        if( x_1*1.0 < 0.5) {
                                            return __pp_r535__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                                return __pp_r536__(x_0, x_1);
                                            } else {
                                                return __pp_r537__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 0.707106781187) {
                                            return __pp_r538__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                                return __pp_r539__(x_0, x_1);
                                            } else {
                                                return __pp_r540__(x_0, x_1);
                                            }
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 1.56524758425) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 1.41421356237) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                            return __pp_r541__(x_0, x_1);
                                        } else {
                                            return __pp_r542__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.11803398875) {
                                            return __pp_r543__(x_0, x_1);
                                        } else {
                                            return __pp_r544__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 1.41421356237) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                            return __pp_r545__(x_0, x_1);
                                        } else {
                                            return __pp_r546__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.11803398875) {
                                            return __pp_r547__(x_0, x_1);
                                        } else {
                                            return __pp_r548__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if( x_0*-1.0 < 2.5) {
                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 0.707106781187) {
                            if( x_1*1.0 < 1.5) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 2.01246117975) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < 1.56524758425) {
                                            return __pp_r549__(x_0, x_1);
                                        } else {
                                            if( x_0*-1.0 < 1.5) {
                                                return __pp_r550__(x_0, x_1);
                                            } else {
                                                return __pp_r551__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < 1.56524758425) {
                                            return __pp_r552__(x_0, x_1);
                                        } else {
                                            if( x_0*-1.0 < 1.5) {
                                                return __pp_r553__(x_0, x_1);
                                            } else {
                                                return __pp_r554__(x_0, x_1);
                                            }
                                        }
                                    }
                                } else {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < 2.45967477525) {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < 2.12132034356) {
                                                return __pp_r555__(x_0, x_1);
                                            } else {
                                                return __pp_r556__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r557__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < 2.45967477525) {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < 2.12132034356) {
                                                return __pp_r558__(x_0, x_1);
                                            } else {
                                                return __pp_r559__(x_0, x_1);
                                            }
                                        } else {
                                            return __pp_r560__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.707106781187+x_1*0.707106781187 < 2.82842712475) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 2.45967477525) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                            return __pp_r561__(x_0, x_1);
                                        } else {
                                            return __pp_r562__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                            return __pp_r563__(x_0, x_1);
                                        } else {
                                            return __pp_r564__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 2.90688837075) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                            return __pp_r565__(x_0, x_1);
                                        } else {
                                            return __pp_r566__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                            return __pp_r567__(x_0, x_1);
                                        } else {
                                            return __pp_r568__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        } else {
                            if( x_1*1.0 < 0.5) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 2.01246117975) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 1.41421356237) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                            return __pp_r569__(x_0, x_1);
                                        } else {
                                            return __pp_r570__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.11803398875) {
                                            return __pp_r571__(x_0, x_1);
                                        } else {
                                            return __pp_r572__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 1.41421356237) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                            return __pp_r573__(x_0, x_1);
                                        } else {
                                            return __pp_r574__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.11803398875) {
                                            return __pp_r575__(x_0, x_1);
                                        } else {
                                            return __pp_r576__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.707106781187+x_1*0.707106781187 < 2.12132034356) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 2.01246117975) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                            return __pp_r577__(x_0, x_1);
                                        } else {
                                            return __pp_r578__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                            return __pp_r579__(x_0, x_1);
                                        } else {
                                            return __pp_r580__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 2.45967477525) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                            return __pp_r581__(x_0, x_1);
                                        } else {
                                            return __pp_r582__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                            return __pp_r583__(x_0, x_1);
                                        } else {
                                            return __pp_r584__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        if( x_0*-0.4472135955+x_1*-0.894427191 < 1.11803398875) {
                            if( x_0*-0.707106781187+x_1*-0.707106781187 < 1.41421356237) {
                                if( x_1*1.0 < 1.5) {
                                    if( x_0*-0.4472135955+x_1*-0.894427191 < 0.22360679775) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < 2.90688837075) {
                                            return __pp_r585__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < 2.82842712475) {
                                                return __pp_r586__(x_0, x_1);
                                            } else {
                                                if(x_0*-0.894427191+x_1*0.4472135955 >= 3.35410196625) return 0; return __pp_r587__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < 2.90688837075) {
                                            return __pp_r588__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < 2.82842712475) {
                                                return __pp_r589__(x_0, x_1);
                                            } else {
                                                if(x_0*-0.894427191+x_1*0.4472135955 >= 3.35410196625) return 0; return __pp_r590__(x_0, x_1);
                                            }
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 0.707106781187) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -0.67082039325) {
                                            if(x_0*-0.894427191+x_1*0.4472135955 >= 3.35410196625) return 0; return __pp_r591__(x_0, x_1);
                                        } else {
                                            if(x_0*-0.894427191+x_1*0.4472135955 >= 3.35410196625) return 0; return __pp_r592__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < -0.22360679775) {
                                            if(x_0*-0.894427191+x_1*0.4472135955 >= 3.35410196625) return 0; return __pp_r593__(x_0, x_1);
                                        } else {
                                            if(x_0*-0.894427191+x_1*0.4472135955 >= 3.35410196625) return 0; return __pp_r594__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 2.90688837075) {
                                    if( x_1*1.0 < 0.5) {
                                        if( x_0*-0.894427191+x_1*0.4472135955 < 2.45967477525) {
                                            return __pp_r595__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.707106781187+x_1*0.707106781187 < 2.12132034356) {
                                                return __pp_r596__(x_0, x_1);
                                            } else {
                                                return __pp_r597__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                            return __pp_r598__(x_0, x_1);
                                        } else {
                                            return __pp_r599__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*0.707106781187 < 2.82842712475) {
                                        if( x_1*1.0 < 0.5) {
                                            return __pp_r600__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                                return __pp_r601__(x_0, x_1);
                                            } else {
                                                return __pp_r602__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 0.67082039325) {
                                            if(x_0*-0.894427191+x_1*0.4472135955 >= 3.35410196625) return 0; return __pp_r603__(x_0, x_1);
                                        } else {
                                            if(x_0*-0.894427191+x_1*0.4472135955 >= 3.35410196625) return 0; return __pp_r604__(x_0, x_1);
                                        }
                                    }
                                }
                            }
                        } else {
                            if( x_1*1.0 < -0.5) {
                                if( x_0*-0.894427191+x_1*0.4472135955 < 2.45967477525) {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.82842712475) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.01246117975) {
                                            return __pp_r605__(x_0, x_1);
                                        } else {
                                            return __pp_r606__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                            return __pp_r607__(x_0, x_1);
                                        } else {
                                            return __pp_r608__(x_0, x_1);
                                        }
                                    }
                                } else {
                                    if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.82842712475) {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.01246117975) {
                                            return __pp_r609__(x_0, x_1);
                                        } else {
                                            return __pp_r610__(x_0, x_1);
                                        }
                                    } else {
                                        if( x_0*-0.4472135955+x_1*-0.894427191 < 2.45967477525) {
                                            if(x_0*-1.0 >= 3.5) return 0; return __pp_r611__(x_0, x_1);
                                        } else {
                                            if(x_0*-1.0 >= 3.5) return 0; return __pp_r612__(x_0, x_1);
                                        }
                                    }
                                }
                            } else {
                                if( x_0*-0.707106781187+x_1*0.707106781187 < 2.12132034356) {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 2.45967477525) {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.12132034356) {
                                            return __pp_r613__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                                return __pp_r614__(x_0, x_1);
                                            } else {
                                                return __pp_r615__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.12132034356) {
                                            return __pp_r616__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                                return __pp_r617__(x_0, x_1);
                                            } else {
                                                return __pp_r618__(x_0, x_1);
                                            }
                                        }
                                    }
                                } else {
                                    if( x_0*-0.894427191+x_1*0.4472135955 < 2.90688837075) {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.12132034356) {
                                            return __pp_r619__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                                return __pp_r620__(x_0, x_1);
                                            } else {
                                                return __pp_r621__(x_0, x_1);
                                            }
                                        }
                                    } else {
                                        if( x_0*-0.707106781187+x_1*-0.707106781187 < 2.12132034356) {
                                            return __pp_r622__(x_0, x_1);
                                        } else {
                                            if( x_0*-0.4472135955+x_1*-0.894427191 < 1.56524758425) {
                                                if(x_0*-1.0 >= 3.5) return 0; return __pp_r623__(x_0, x_1);
                                            } else {
                                                if(x_0*-1.0 >= 3.5) return 0; return __pp_r624__(x_0, x_1);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        return 0;
    }

};
}

#endif 
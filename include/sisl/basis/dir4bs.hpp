
#include <math.h>
#include <sisl/basis_function.hpp>
#include <sisl/primitives.hpp>
#include <sisl/lattice.hpp>
#include <vector>

#ifndef _TP_2D_4DBS_H_
#define _TP_2D_4DBS_H_

namespace sisl{
template <class O, class I>
class dir4bs : public basis_function <O,I> {


public:
    static std::string getBasisName(){
        return std::string("dir4bs");
    }
    static const O M(const I &x, const I &y) { return (O) box_spline(x,y); };
    static const O M(const vector2<I> &p) { return (O) box_spline(p.i, p.j); };

    // This function should return the intersection of the closure of the support of
    // the generator and the lattice.
    static std::vector<std::tuple<int,int>> getSupport() {
        using namespace std;
        std::vector<std::tuple<int,int>> support;
        for(int i = -2; i <= 2; i++)
            for(int j = -2; j <= 2; j++) {
                support.push_back(make_tuple(i,j));
            }
        return support;
    };

    static std::vector<std::tuple<int,int>> getEffectiveSupport(const vector2<I> &p) {
        using namespace std;
        std::vector<std::tuple<int,int>> support;
        for(int i = -2; i <= 2; i++)
            for(int j = -2; j <= 2; j++) {
                    support.push_back(make_tuple(i,j));
                }
        return support;
    }

    static std::vector<std::tuple<int,int,O>> getBeppoLevi2Norm() {throw "basis_function()::getBeppoLevi2Norm() - Not Implemented!";};
    static std::vector<std::tuple<int,int,O>> getBeppoLevi1Norm(){throw "basis_function()::getBeppoLevi1Norm() - Not Implemented!";};
    static std::vector<std::tuple<int,int,O>> autoCorrelation(){throw "basis_function()::autoCorrelation() - Not Implemented!";};

    inline static const O convolutionSum(const vector2<I> &p, const shift_invariant_space2< dir4bs, O, I> *lattice) {
        O sum = 0;
        I dh = lattice->getScale();
        vector2<I> vox(p.i/dh, p.j/dh);     

        int vx = (int)floor(vox.i), 
            vy = (int)floor(vox.j);

        for(int i = -2; i <= 2; i++)
            for(int j = -2; j <= 2; j++) {
                sum += lattice->GV(vx + i, vy + j) *
                        dir4bs::M(vox.i - I(vx + i), vox.j - I(vy + j));
            }
        return sum;     
    }

private:
    static inline O __pp_r1__(const I &x_0, const I &x_1) {
        O __pp_r1___result;
        __pp_r1___result = (1.0L/2.0L)*pow(x_0, 2) - x_0*x_1 - 3.0L/2.0L*x_0 + (1.0L/2.0L)*pow(x_1, 2) + (3.0L/2.0L)*x_1 + 9.0L/8.0L;
        return __pp_r1___result;
    }
    static inline O __pp_r2__(const I &x_0, const I &x_1) {
        O __pp_r2___result;
        __pp_r2___result = -1.0L/2.0L*pow(x_0, 2) - 1.0L/2.0L*x_0 + (1.0L/4.0L)*pow(x_1, 2) + x_1 + 7.0L/8.0L;
        return __pp_r2___result;
    }
    static inline O __pp_r3__(const I &x_1) {
        O __pp_r3___result;
        __pp_r3___result = (1.0L/4.0L)*pow(x_1, 2) + x_1 + 1;
        return __pp_r3___result;
    }
    static inline O __pp_r4__(const I &x_0, const I &x_1) {
        O __pp_r4___result;
        __pp_r4___result = pow(x_0, 2) - x_0*x_1 - 2*x_0 + (1.0L/4.0L)*pow(x_1, 2) + x_1 + 1;
        return __pp_r4___result;
    }
    static inline O __pp_r5__(const I &x_0, const I &x_1) {
        O __pp_r5___result;
        __pp_r5___result = (1.0L/2.0L)*pow(x_0, 2) - x_0*x_1 - 3.0L/2.0L*x_0 + (1.0L/4.0L)*pow(x_1, 2) + x_1 + 7.0L/8.0L;
        return __pp_r5___result;
    }
    static inline O __pp_r6__(const I &x_0, const I &x_1) {
        O __pp_r6___result;
        __pp_r6___result = -1.0L/2.0L*pow(x_0, 2) - 1.0L/2.0L*x_0 + (1.0L/2.0L)*x_1 + 5.0L/8.0L;
        return __pp_r6___result;
    }
    static inline O __pp_r7__(const I &x_0, const I &x_1) {
        O __pp_r7___result;
        __pp_r7___result = -pow(x_0, 2) + x_0*x_1 - 1.0L/2.0L*pow(x_1, 2) + 1.0L/2.0L;
        return __pp_r7___result;
    }
    static inline O __pp_r8__(const I &x_1) {
        O __pp_r8___result;
        __pp_r8___result = (1.0L/4.0L)*pow(x_1, 2) + x_1 + 1;
        return __pp_r8___result;
    }
    static inline O __pp_r9__(const I &x_0, const I &x_1) {
        O __pp_r9___result;
        __pp_r9___result = -1.0L/2.0L*pow(x_0, 2) + x_0*x_1 + (1.0L/2.0L)*x_0 - 1.0L/4.0L*pow(x_1, 2) + (1.0L/2.0L)*x_1 + 7.0L/8.0L;
        return __pp_r9___result;
    }
    static inline O __pp_r10__(const I &x_0) {
        O __pp_r10___result;
        __pp_r10___result = (1.0L/2.0L)*pow(x_0, 2) + (3.0L/2.0L)*x_0 + 9.0L/8.0L;
        return __pp_r10___result;
    }
    static inline O __pp_r11__(const I &x_0, const I &x_1) {
        O __pp_r11___result;
        __pp_r11___result = -pow(x_0, 2) + x_0*x_1 - 1.0L/2.0L*pow(x_1, 2) + 1.0L/2.0L;
        return __pp_r11___result;
    }
    static inline O __pp_r12__(const I &x_0, const I &x_1) {
        O __pp_r12___result;
        __pp_r12___result = -1.0L/2.0L*pow(x_0, 2) + x_0*x_1 + (1.0L/2.0L)*x_0 - 1.0L/2.0L*pow(x_1, 2) + 5.0L/8.0L;
        return __pp_r12___result;
    }
    static inline O __pp_r13__(const I &x_0, const I &x_1) {
        O __pp_r13___result;
        __pp_r13___result = (1.0L/2.0L)*pow(x_0, 2) + (3.0L/2.0L)*x_0 - 1.0L/4.0L*pow(x_1, 2) - 1.0L/2.0L*x_1 + 7.0L/8.0L;
        return __pp_r13___result;
    }
    static inline O __pp_r14__(const I &x_0, const I &x_1) {
        O __pp_r14___result;
        __pp_r14___result = pow(x_0, 2) - x_0*x_1 + 2*x_0 + (1.0L/4.0L)*pow(x_1, 2) - x_1 + 1;
        return __pp_r14___result;
    }
    static inline O __pp_r15__(const I &x_0, const I &x_1) {
        O __pp_r15___result;
        __pp_r15___result = pow(x_0, 2) - x_0*x_1 - 2*x_0 + (1.0L/4.0L)*pow(x_1, 2) + x_1 + 1;
        return __pp_r15___result;
    }
    static inline O __pp_r16__(const I &x_0, const I &x_1) {
        O __pp_r16___result;
        __pp_r16___result = (1.0L/2.0L)*pow(x_0, 2) - 3.0L/2.0L*x_0 - 1.0L/4.0L*pow(x_1, 2) + (1.0L/2.0L)*x_1 + 7.0L/8.0L;
        return __pp_r16___result;
    }
    static inline O __pp_r17__(const I &x_0, const I &x_1) {
        O __pp_r17___result;
        __pp_r17___result = -1.0L/2.0L*pow(x_0, 2) + x_0*x_1 - 1.0L/2.0L*x_0 - 1.0L/2.0L*pow(x_1, 2) + 5.0L/8.0L;
        return __pp_r17___result;
    }
    static inline O __pp_r18__(const I &x_0, const I &x_1) {
        O __pp_r18___result;
        __pp_r18___result = -pow(x_0, 2) + x_0*x_1 - 1.0L/2.0L*pow(x_1, 2) + 1.0L/2.0L;
        return __pp_r18___result;
    }
    static inline O __pp_r19__(const I &x_0) {
        O __pp_r19___result;
        __pp_r19___result = (1.0L/2.0L)*pow(x_0, 2) - 3.0L/2.0L*x_0 + 9.0L/8.0L;
        return __pp_r19___result;
    }
    static inline O __pp_r20__(const I &x_0, const I &x_1) {
        O __pp_r20___result;
        __pp_r20___result = -1.0L/2.0L*pow(x_0, 2) + x_0*x_1 - 1.0L/2.0L*x_0 - 1.0L/4.0L*pow(x_1, 2) - 1.0L/2.0L*x_1 + 7.0L/8.0L;
        return __pp_r20___result;
    }
    static inline O __pp_r21__(const I &x_1) {
        O __pp_r21___result;
        __pp_r21___result = (1.0L/4.0L)*pow(x_1, 2) - x_1 + 1;
        return __pp_r21___result;
    }
    static inline O __pp_r22__(const I &x_0, const I &x_1) {
        O __pp_r22___result;
        __pp_r22___result = -pow(x_0, 2) + x_0*x_1 - 1.0L/2.0L*pow(x_1, 2) + 1.0L/2.0L;
        return __pp_r22___result;
    }
    static inline O __pp_r23__(const I &x_0, const I &x_1) {
        O __pp_r23___result;
        __pp_r23___result = -1.0L/2.0L*pow(x_0, 2) + (1.0L/2.0L)*x_0 - 1.0L/2.0L*x_1 + 5.0L/8.0L;
        return __pp_r23___result;
    }
    static inline O __pp_r24__(const I &x_0, const I &x_1) {
        O __pp_r24___result;
        __pp_r24___result = (1.0L/2.0L)*pow(x_0, 2) - x_0*x_1 + (3.0L/2.0L)*x_0 + (1.0L/4.0L)*pow(x_1, 2) - x_1 + 7.0L/8.0L;
        return __pp_r24___result;
    }
    static inline O __pp_r25__(const I &x_0, const I &x_1) {
        O __pp_r25___result;
        __pp_r25___result = pow(x_0, 2) - x_0*x_1 + 2*x_0 + (1.0L/4.0L)*pow(x_1, 2) - x_1 + 1;
        return __pp_r25___result;
    }
    static inline O __pp_r26__(const I &x_1) {
        O __pp_r26___result;
        __pp_r26___result = (1.0L/4.0L)*pow(x_1, 2) - x_1 + 1;
        return __pp_r26___result;
    }
    static inline O __pp_r27__(const I &x_0, const I &x_1) {
        O __pp_r27___result;
        __pp_r27___result = -1.0L/2.0L*pow(x_0, 2) + (1.0L/2.0L)*x_0 + (1.0L/4.0L)*pow(x_1, 2) - x_1 + 7.0L/8.0L;
        return __pp_r27___result;
    }
    static inline O __pp_r28__(const I &x_0, const I &x_1) {
        O __pp_r28___result;
        __pp_r28___result = (1.0L/2.0L)*pow(x_0, 2) - x_0*x_1 + (3.0L/2.0L)*x_0 + (1.0L/2.0L)*pow(x_1, 2) - 3.0L/2.0L*x_1 + 9.0L/8.0L;
        return __pp_r28___result;
    }

    static O box_spline(const I &x_0, const I &x_1) {
        if( x_1*1.0 < 0.0) {
            if( x_0*-0.894427191+x_1*0.4472135955 < 0.0) {
                if( x_1*1.0 < -1.0) {
                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.4472135955) {
                        if(x_0*-0.707106781187+x_1*0.707106781187 <= -1.06066017178) return 0; return __pp_r1__(x_0, x_1);
                    } else {
                        if( x_0*-1.0 < 0.5) {
                            return __pp_r2__(x_0, x_1);
                        } else {
                            if(x_1*1.0 <= -2.0) return 0; return __pp_r3__(x_1);
                        }
                    }
                } else {
                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.4472135955) {
                        if( x_0*-1.0 < -0.5) {
                            if(x_0*-0.894427191+x_1*0.4472135955 <= -0.894427191) return 0; return __pp_r4__(x_0, x_1);
                        } else {
                            return __pp_r5__(x_0, x_1);
                        }
                    } else {
                        if( x_0*-0.707106781187+x_1*0.707106781187 < -0.353553390593) {
                            return __pp_r6__(x_0, x_1);
                        } else {
                            return __pp_r7__(x_0, x_1);
                        }
                    }
                }
            } else {
                if( x_1*1.0 < -1.0) {
                    if( x_0*-0.894427191+x_1*0.4472135955 < 0.4472135955) {
                        if( x_0*-0.707106781187+x_1*0.707106781187 < -0.353553390593) {
                            if(x_1*1.0 <= -2.0) return 0; return __pp_r8__(x_1);
                        } else {
                            return __pp_r9__(x_0, x_1);
                        }
                    } else {
                        if(x_0*-1.0 >= 1.5) return 0; return __pp_r10__(x_0);
                    }
                } else {
                    if( x_0*-0.894427191+x_1*0.4472135955 < 0.4472135955) {
                        if( x_0*-1.0 < 0.5) {
                            return __pp_r11__(x_0, x_1);
                        } else {
                            return __pp_r12__(x_0, x_1);
                        }
                    } else {
                        if( x_0*-0.707106781187+x_1*0.707106781187 < 0.353553390593) {
                            return __pp_r13__(x_0, x_1);
                        } else {
                            if(x_0*-0.894427191+x_1*0.4472135955 >= 0.894427191) return 0; return __pp_r14__(x_0, x_1);
                        }
                    }
                }
            }
        } else {
            if( x_0*-0.894427191+x_1*0.4472135955 < 0.0) {
                if( x_1*1.0 < 1.0) {
                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.4472135955) {
                        if( x_0*-0.707106781187+x_1*0.707106781187 < -0.353553390593) {
                            if(x_0*-0.894427191+x_1*0.4472135955 <= -0.894427191) return 0; return __pp_r15__(x_0, x_1);
                        } else {
                            return __pp_r16__(x_0, x_1);
                        }
                    } else {
                        if( x_0*-1.0 < -0.5) {
                            return __pp_r17__(x_0, x_1);
                        } else {
                            return __pp_r18__(x_0, x_1);
                        }
                    }
                } else {
                    if( x_0*-0.894427191+x_1*0.4472135955 < -0.4472135955) {
                        if(x_0*-1.0 <= -1.5) return 0; return __pp_r19__(x_0);
                    } else {
                        if( x_0*-0.707106781187+x_1*0.707106781187 < 0.353553390593) {
                            return __pp_r20__(x_0, x_1);
                        } else {
                            if(x_1*1.0 >= 2.0) return 0; return __pp_r21__(x_1);
                        }
                    }
                }
            } else {
                if( x_1*1.0 < 1.0) {
                    if( x_0*-0.894427191+x_1*0.4472135955 < 0.4472135955) {
                        if( x_0*-0.707106781187+x_1*0.707106781187 < 0.353553390593) {
                            return __pp_r22__(x_0, x_1);
                        } else {
                            return __pp_r23__(x_0, x_1);
                        }
                    } else {
                        if( x_0*-1.0 < 0.5) {
                            return __pp_r24__(x_0, x_1);
                        } else {
                            if(x_0*-0.894427191+x_1*0.4472135955 >= 0.894427191) return 0; return __pp_r25__(x_0, x_1);
                        }
                    }
                } else {
                    if( x_0*-0.894427191+x_1*0.4472135955 < 0.4472135955) {
                        if( x_0*-1.0 < -0.5) {
                            if(x_1*1.0 >= 2.0) return 0; return __pp_r26__(x_1);
                        } else {
                            return __pp_r27__(x_0, x_1);
                        }
                    } else {
                        if(x_0*-0.707106781187+x_1*0.707106781187 >= 1.06066017178) return 0; return __pp_r28__(x_0, x_1);
                    }
                }
            }
        }
        return 0;
    }
};
}

#endif // _TP_2D_4DBS_H_
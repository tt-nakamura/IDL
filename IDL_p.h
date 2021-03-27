#ifndef __IDL_p_h__
#define __IDL_p_h__

#include<NTL/ZZ_pX.h>
#include<NTL/mat_lzz_p.h>
#include "IDL.h"
using namespace NTL;

// ideal modulo prime number p
struct IDL_p : mat_zz_p {// matrix representation
    static Mat<vec_zz_p> M;
    static void init();
};

typedef vec_zz_p IDL1_p;// generator representation

inline void set(IDL_p& X) { ident(X, IDL::degree()); }// X=1

void conv(IDL2&, const ZZ_pX&);
void conv(IDL2&, const IDL_p&, long);
void conv(IDL1&, const IDL1_p&);
void conv(IDL1_p&, const IDL_p&, const ZZ&);

void add(IDL_p&, const IDL_p&, const IDL_p&);
void add(IDL_p&, const IDL_p&, const IDL1_p&);
void mul(IDL_p&, const IDL_p&, const IDL_p&);
void mul(IDL_p&, const IDL_p&, const IDL1_p&);
void sqr(IDL_p&, const IDL_p&);
void div(IDL_p&, const IDL_p&, const IDL_p&);

inline void operator+=(IDL_p& X, const IDL_p& A) { add(X,X,A); }
inline void operator*=(IDL_p& X, const IDL_p& A) { mul(X,X,A); }
inline void operator+=(IDL_p& X, const IDL1_p& a) { add(X,X,a); }
inline void operator*=(IDL_p& X, const IDL1_p& a) { mul(X,X,a); }

inline void conv(IDL_p& X, const IDL1_p& a) { set(X), X*=a; }

void radical(IDL_p&);
long valuation(const IDL&, const IDL_p&);

void PowerMod(ZZX&, const ZZX&, long, const ZZX&);
void SuplBase(mat_zz_p&, const mat_zz_p&);
void SuplBase(mat_zz_p&, const mat_zz_p&, const mat_zz_p&);
void InvImage(mat_zz_p&, const mat_zz_p&, const mat_zz_p&);

#endif // __IDL_p__h
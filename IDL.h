#ifndef __IDL_h__
#define __IDL_h__

#include<NTL/ZZX.h>
#include<NTL/HNF.h>
using namespace NTL;

// integral ideal

struct IDL : mat_ZZ {// HNF matrix representation
    static mat_ZZ W,W1;
    static ZZ D,d,d1;
    static ZZX T;
    static Vec<ZZX> w;
    static Mat<vec_ZZ> M;
    static void init(const ZZX&);
    static inline long degree() { return deg(T); }
    static inline const mat_ZZ& IntBases() { return W; }
    static inline const ZZ& denominator() { return d; }
    static inline const ZZ& FieldDiscr() { return D; }
    static void index(ZZ&);
    static void divisors(Vec<long>&);
    ZZ N;// norm
    void SetNorm();
    inline const ZZ& norm() const { return N; }
};

struct IDL1 : vec_ZZ {// principal ideal
    ZZ N;// norm
    void SetNorm();
    inline const ZZ& norm() const { return N; }
};

struct IDL2 : IDL1 {// two element representation
    ZZ a;
};

inline void set(IDL& X) { ident(X, IDL::degree()), set(X.N); }
inline long IsOne(const IDL& A) { return IsIdent(A, IDL::degree()); }
inline long IsZero(const IDL& A) { return IsZero(A.N); }
inline long IsZero(const IDL1& a) { return IsZero(a.N); }

static ZZ ZZ_one(1);

void set(IDL&, const ZZ&);
void conv(IDL1&, const ZZX&, const ZZ& = ZZ_one);
void conv(IDL2&, const ZZ&, const ZZX&, const ZZ& = ZZ_one);
void conv(IDL&, const ZZX&, const ZZ& = ZZ_one);
void conv(IDL&, const ZZ&, const ZZX&, const ZZ& = ZZ_one);

void add(IDL&, const IDL&, const IDL&);
void add(IDL&, const IDL&, const IDL1&);
void mul(IDL&, const IDL&, const IDL&);
void mul(IDL&, const IDL&, const IDL1&);
void mul(IDL&, const IDL&, const IDL2&);
void sqr(IDL&, const IDL&);
void power(IDL&, const IDL&, const ZZ&);
void power(IDL&, const IDL&, long);

inline void operator+=(IDL& X, const IDL& A) { add(X,X,A); }
inline void operator*=(IDL& X, const IDL& A) { mul(X,X,A); }
inline void operator+=(IDL& X, const IDL1& a) { add(X,X,a); }
inline void operator*=(IDL& X, const IDL1& a) { mul(X,X,a); }
inline void operator*=(IDL& X, const IDL2& A) { mul(X,X,A); }

inline void conv(IDL& X, const IDL1& A) { set(X), X*=A; }
inline void conv(IDL& X, const IDL2& A) { set(X), X*=A; }

std::ostream& operator<<(std::ostream&, const IDL2&);

#endif // __IDL_h__
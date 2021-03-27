#ifndef __SepAlg_h__
#define __SepAlg_h__

#include "IDL_p.h"

struct SepAlg : IDL_p {// Separable Algebra
    static Mat<vec_zz_p> M;
    static mat_zz_p W;
    static mat_zz_p V;
    static void init(const IDL_p&);
    inline static long dim() { return W.NumRows(); }
    inline static long IsField() { return V.NumRows()==1; }
};

long idempot(SepAlg&);

#endif
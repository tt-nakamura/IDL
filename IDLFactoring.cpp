// uses NTL
//   http://www.shoup.net/ntl
// reference:
//   H. Cohen, "A Course in Computational Algebraic Number theory"
//     section 6.2

#include "SepAlg.h"
#include<NTL/ZZ_pXFactoring.h>
using namespace NTL;

void factor(Vec<Pair<IDL2, long> >& x, const ZZ& p)
// prime ideal factorization of prime number p
// x = (prime ideal, exponent) pair of factors of p
// reference: Cohen, algorithm 6.2.9
{
    long i,j,k,l,n(IDL::degree());
    ZZ_pX g;
    IDL2 A;
    SepAlg U;
    Vec<IDL_p> H,J,K;
    Vec<Pair<ZZ_pX, long> > y;

    ZZ_p::init(p);
    if(!divide(IDL::d, p)) {
        // section 4.8.2 of Cohen
        conv(g, IDL::T);
        CanZass(y,g);
        x.SetLength(y.length());
        for(i=0; i<y.length(); i++) {
            conv(x[i].a, y[i].a);
            x[i].b = y[i].b;
        }
        return;
    }
    conv(k,p);
    zz_p::init(k);
    IDL_p::init();

    K.SetLength(1);
    radical(K[0]);
    // section 6.2.2 of Cohen
    for(k=1; !IsZero(K[k-1]); k++) {
        K.SetLength(k+1);
        mul(K[k], K[k-1], K[0]);
    }
    J.SetLength(k);
    H.SetLength(k);
    J[0] = K[0];
    for(i=1; i<k; i++) {
        div(J[i], K[i], K[i-1]);
        div(H[i-1], J[i-1], J[i]);
    }
    H[k-1] = J[k-1];
    // section 6.2.4 of Cohen
    for(j=k=0; j<H.length(); j++) {
        if(H[j].NumRows()==n) continue;
        J.SetLength(1);
        J[0] = H[j];
        for(l=0; l>=0; l--) {
            SepAlg::init(J[l]);
            if(idempot(U)) {
                J.SetLength(l + U.NumRows());
                for(i=U.NumRows()-1; i>=0; i--)
                    add(J[l+i], J[l], U[i]);
                l = J.length();
            }
            else {
                x.SetLength(k+1);
                conv(x[k].a, J[l], SepAlg::dim());
                x[k++].b = j+1;
            }
        }
    }
}

void factor2(Vec<Pair<IDL2, long> >&x, const ZZ& p)
// another method of prime ideal factorization
// reference: Cohen, exercise 11 of chapter 6
{
    long i,j,k,n(IDL::degree());
    ZZ_pX g;
    IDL P;
    SepAlg U;
    Vec<IDL_p> H;
    Vec<Pair<ZZ_pX, long> > y;

    ZZ_p::init(p);
    if(!divide(IDL::d, p)) {
        conv(g, IDL::T);
        CanZass(y,g);
        x.SetLength(y.length());
        for(i=0; i<y.length(); i++) {
            conv(x[i].a, y[i].a);
            x[i].b = y[i].b;
        }
        return;
    }
    conv(k,p);
    zz_p::init(k);
    IDL_p::init();

    H.SetLength(1);
    radical(H[0]);
    set(P,p);

    for(j=k=0; j>=0; j--) {
        SepAlg::init(H[j]);
        if(idempot(U)) {
            H.SetLength(j + U.NumRows());
            for(i=U.NumRows()-1; i>=0; i--)
                add(H[i+j], H[j], U[i]);
            j = H.length();
        }
        else {
            x.SetLength(k+1);
            conv(x[k].a, H[j], SepAlg::dim());
            x[k++].b = valuation(P, H[j]);
        }
    }
}

void mul(IDL& X, const Vec<Pair<IDL2, long> >& v)
// X = product of (prime ideal)^(exponent)
// (inverse of factorization)
{
    long i,j;
    set(X);
    for(i=0; i<v.length(); i++)
        for(j=0; j<v[i].b; j++) X *= v[i].a;
}

void mul(IDL& X, const Vec<Pair<IDL, long> >& v)
// X = product of (prime ideal)^(exponent)
// (inverse of factorization)
{
    long i,j;
    IDL A;
    set(X);
    for(i=0; i<v.length(); i++) {
        power(A, v[i].a, v[i].b);
        X *= A;
    }
}
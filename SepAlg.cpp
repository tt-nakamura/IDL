// uses NTL
//   http://www.shoup.net/ntl
// reference:
//   H. Cohen, "A Course in Computational Algebraic Number theory"
//     section 6.2.4

#include "SepAlg.h"
#include<NTL/lzz_pXFactoring.h>

Mat<vec_zz_p> SepAlg::M;// multiplication table
mat_zz_p SepAlg::W;// Bases
mat_zz_p SepAlg::V;// Fermat map kernel

void SepAlg::init(const IDL_p& H)
// initialize separable algebra to split H
// H = square free product of prime ideals above p
// reference: Cohen, algorithm 6.2.9 steps 8--10
{
    long i,j,k,l,n(IDL::degree()),r(H.NumRows()),m(n-r);
    long p(zz_p::modulus());
    zz_p a;
    vec_zz_p v;
    mat_zz_p A(H),B;
    M.kill();
    M.SetDims(m,m);
    W.SetDims(m,n);
    A.SetDims(r+1,n);
    B.SetDims(n,m);
    clear(A[r]);
    set(A[r][0]);
    SuplBase(A,A);
    for(i=0; i<m; i++) W[i] = A[r+i];
    inv(A,A);
    for(i=0; i<n; i++)
        for(j=0; j<m; j++) B[i][j] = A[i][r+j];
    // make multiplication table
    for(i=0; i<m; i++) {
        for(j=0; j<=i; j++) {
            M[i][j].SetLength(n);
            for(k=0; k<n; k++) {
                for(l=0; l<n; l++) {
                    mul(a, W[i][k], W[j][l]);
                    mul(v, a, IDL_p::M[k][l]);
                    M[i][j] += v;
                }
            }
            M[i][j] *= B;
        }
        for(j=0; j<i; j++) M[j][i] = M[i][j];
    }
    // Fermat map
    A.SetDims(m,m);
    B.SetDims(m,m);
    for(i=0; i<m; i++) {
        for(j=0; j<m; j++) A[j] = M[i][j];
        power(A,A,p);
        B[i] = A[0];
        B[i][i]--;
    }
    kernel(V,B);
}

static void mul(vec_zz_p& x, const vec_zz_p& a, const vec_zz_p& b)
// x = a*b; assume &x!=&a, &x!=&b, and x=(0,0,...,0)
{
    long i,j;
    zz_p u;
    vec_zz_p v;
    for(i=0; i<a.length(); i++) {
        for(j=0; j<b.length(); j++) {
            mul(u, a[i], b[j]);
            mul(v, u, SepAlg::M[i][j]);
            x += v;
        }
    }
}

static long deg(const vec_zz_p& a) {// degree of polynomial
    long i;
    for(i=a.length()-1; i>0; i--)
        if(!IsZero(a[i])) break;
    return i;
}

long idempot(SepAlg& E)
// find idempotent E in separable algebra
//   (i.e., E^2=E and E!=0 and E!=1)
// return:
//   0 if not found because the algebra is a field
//   1 if found (and E is set as idempotent)
// reference:
//   Cohen, algorithm 6.2.9 steps 12--14
{
    if(SepAlg::V.NumRows()==1) return 0;
    long i,j,k,l,m(SepAlg::dim());
    zz_pX f,g,u,v;
    mat_zz_p A,B;
    Vec<Pair<zz_pX, long> > t;

    for(j=0; j<SepAlg::V.NumRows(); j++)
        if(deg(SepAlg::V[j])) break;
    
    A.SetDims(m+1,m);
    set(A[0][0]);
    for(i=0; i<m; i++)
        mul(A[i+1], A[i], SepAlg::V[j]);
    kernel(B,A);
    for(i=0, k=m+1; i<B.NumRows(); i++)
        if((j=deg(B[i])) < k) { k=j; l=i; }

    conv(f, B[l]);
    MakeMonic(f);// MinPoly

    CanZass(t,f);
    E.SetDims(t.length(), m+1);
    for(i=0; i<t.length(); i++) {
        div(g,f,t[i].a);
        XGCD(g,u,v,g,t[i].a);
        v *= t[i].a;
        VectorCopy(E[i], v.rep, m+1);
    }
    E *= A;
    E *= SepAlg::W;
    return 1;
}
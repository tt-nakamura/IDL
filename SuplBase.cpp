// uses NTL
//   http://www.shoup.net/ntl
// reference:
//   H. Cohen, "A Course in Computational Algebraic Number Theory"

#include<NTL/mat_lzz_p.h>
using namespace NTL;

void SuplBase(mat_zz_p& B, const mat_zz_p& M)
// input:
//   M = m by n matrix of rank m (m<n)
// output:
//   B = n by n invertible matrix
//       first m rows of B are those of M
// reference: Cohen, algorithm 2.3.6
{
    long i,j,k,l,m(M.NumRows()), n(M.NumCols());
    zz_p a,b;
    mat_zz_p A(M),I;
    ident(I,n);
    for(k=0; k<m; k++) {
        for(l=k; l<n; l++)
            if(!IsZero(A[k][l])) break;
        if(l==n) Error("rows of M are dependent");
        if(l>k) {
            swap(I[l], I[k]);
            swap(A[k][l], A[k][k]);
        }
        inv(a, A[k][k]);
        for(i=k+1; i<m; i++) {
            if(l>k) swap(A[i][k], A[i][l]);
            if(IsZero(A[i][k])) continue;
            A[i][k] *= a;
            for(j=0; j<n; j++) {
                if(j==k || IsZero(A[k][j])) continue;
                mul(b, A[i][k], A[k][j]);
                A[i][j] -= b;
            }
        }
    }
    if(&B!=&M) B=M;
    B.SetDims(n,n);
    for(i=m; i<n; i++) B[i] = I[i];
}

void InvImage(mat_zz_p& X, const mat_zz_p& M, const mat_zz_p& V)
// input:
//   M = m by n matrix of rank m (m<=n)
//   V = r by n matrix whose rows are in row space of M
// output:
//   X = r by m matrix as solutions of XM = V
// reference: Cohen, algorithm 2.3.5
{
    long i,j,k,l,m(M.NumRows()),n(M.NumCols()),r(V.NumRows());
    mat_zz_p A(M),B(V);
    zz_p a,b,c;
    for(k=0; k<m; k++) {
        for(l=k; l<n; l++)
            if(!IsZero(A[k][l])) break;
        if(l==n) Error("rows of M are dependent");
        if(l>k) {
            for(i=k; i<m; i++) swap(A[i][k], A[i][l]);
            for(i=0; i<r; i++) swap(B[i][k], B[i][l]);
        }
        inv(a,A[k][k]);
        for(j=k+1; j<n; j++) {
            mul(b,a,A[k][j]);
            for(i=k+1; i<m; i++) {
                mul(c,b,A[i][k]);
                A[i][j] -= c;
            }
            for(i=0; i<r; i++) {
                mul(c,b,B[i][k]);
                B[i][j] -= c;
            }
        }
    }
    for(i=0; i<r; i++)
        for(j=m; j<n; j++)
            if(!IsZero(B[i][j]))
                Error("InvImage failed");
    X.SetDims(r,m);
    for(i=0; i<r; i++) {
        for(j=m-1; j>=0; j--) {
            X[i][j] = B[i][j];
            for(k=j+1; k<m; k++) {
                mul(a,X[i][k],A[k][j]);
                X[i][j] -= a;
            }
            X[i][j] /= A[j][j];
        }
    }
}

void SuplBase(mat_zz_p& B, const mat_zz_p& V, const mat_zz_p& M)
// input:
//   M = m by n matrix of rank m (m<n)
//   V = r by n matrix (r<m) whose rows are in row space of M
// output:
//   B = m-r by n matrix whose rows are in row space of M and
//         are independent of rows of V
// reference: Cohen, algorithm 2.3.7
{
    long i,j,k,n(M.NumCols()), m(M.NumRows()), r(V.NumRows());
    zz_p a;
    mat_zz_p X;
    InvImage(X,M,V);
    SuplBase(X,X);
    for(i=r; i<m; i++) X[i-r] = X[i];
    X.SetDims(m-r,m);
    mul(B,X,M);
}
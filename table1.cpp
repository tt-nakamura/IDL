#include "IDLFactoring.h"

main() {
    long i,n,b(10);
    ZZX T;
    ZZ p;
    Vec<Pair<IDL2,long> > f;
    Vec<long> d;
    IDL A,P;
    for(n=3; n<=14; n++) {
        do {
            BuildIrred(T,n,b);
            IDL::init(T);
            IDL::divisors(d);
        } while(d.length()==0);
        std::cout << T << '\t';
        p = d[d.length()-1];
        factor(f,p);
        std::cout << f << '\n';
        set(P,p);
        mul(A,f);
        if(A!=P) break;
    }
}
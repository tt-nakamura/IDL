#include "IDL.h"
#include<NTL/pair.h>

void BuildIrred(ZZX&, long, long);
void factor(Vec<Pair<IDL2, long> >&, const ZZ&);
void factor2(Vec<Pair<IDL2, long> >&, const ZZ&);
void mul(IDL&, const Vec<Pair<IDL, long> >&);
void mul(IDL&, const Vec<Pair<IDL2, long> >&);
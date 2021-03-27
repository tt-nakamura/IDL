NTL = -lntl -lgmp -L/usr/local/lib

OBJ = IDLFactoring.o round2.o ZZlib.o SuplBase.o ZZXlib.o ZZFactoring.o rho.o mpqs.o IDL.o IDL_p.o SepAlg.o

table1: table1.o $(OBJ)
	g++ $(OBJ) table1.o $(NTL)

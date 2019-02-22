output: Source.o nrerror2new.o vectornew.o  dvectornew.o  dmatrix2new.o fill02new.o free_vectornew.o free_dvectornew.o  free_dmatrix2new.o trapzdnew.o integratnew.o mgfas4new.o mgfas5new.o nrfunc3new.o nrfunc4new.o Stat2new.o Stat3new.o
	g++ Source.o nrerror2new.o vectornew.o dvectornew.o dmatrix2new.o fill02new.o free_vectornew.o free_dvectornew.o free_dmatrix2new.o trapzdnew.o integratnew.o mgfas4new.o mgfas5new.o nrfunc3new.o nrfunc4new.o Stat2new.o Stat3new.o -o output
Source.o: Source.cpp
	g++ -c Source.cpp
nrerror2new.o: nrerror2new.cpp Header1new.h
	g++ -c nrerror2new.cpp
vectornew.o: vectornew.cpp Header2new.h
	g++ -c vectornew.cpp
dvectornew.o: dvectornew.cpp Header2new.h
	g++ -c dvectornew.cpp
free_vectornew.o: free_vectornew.cpp Header2new.h
	g++ -c free_vectornew.cpp
free_dvectornew.o: free_dvectornew.cpp Header2new.h
	g++ -c free_dvectornew.cpp
dmatrix2new.o: dmatrix2new.cpp Header2new.h
	g++ -c dmatrix2new.cpp
free_dmatrix2new.o: free_dmatrix2new.cpp Header2new.h
	g++ -c free_dmatrix2new.cpp
fill02new.o: fill02new.cpp Header2new.h
	g++ -c fill02new.cpp
trapzdnew.o: trapzdnew.cpp Header5new.h
	g++ -c trapzdnew.cpp
integratnew.o: integratnew.cpp Header6new.h
	g++ -c integratnew.cpp
mgfas4new.o: mgfas4new.cpp Header7new.h
	g++ -c mgfas4new.cpp
mgfas5new.o: mgfas5new.cpp Header7new.h
	g++ -c mgfas5new.cpp
nrfunc3new.o: nrfunc3new.cpp Header7new.h
	g++ -c nrfunc3new.cpp
nrfunc4new.o: nrfunc4new.cpp Header7new.h
	g++ -c nrfunc4new.cpp
Stat2new.o: Stat2new.cpp Header7new.h
	g++ -c Stat2new.cpp
Stat3new.o: Stat3new.cpp Header7new.h
	g++ -c Stat3new.cpp
clean:
	rm *.o output

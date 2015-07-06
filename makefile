//COMPILER = g++
COMPILER = /usr/local/bin/g++-4.9 #Because XCode is Lame
OFLAG = -O3
FLAGS = -m64 -g -fopenmp $(OFLAG)
LINKS = -lfftw3l_threads -lfftw3l #-L/opt/local/lib -lfftw3l_omp

do: fclean dirmake oclean compile 

xcmake : do
	
compile: g2header.h g2parameters.h g2model.o g2functions.o g2spectra.o g2power.o g2output.o g2main.o g2init.o 
	$(COMPILER) $(FLAGS) g2model.o g2functions.o g2spectra.o g2power.o g2output.o g2init.o g2main.o $(LINKS) -o gabe

clean : fclean oclean

g2main.o: g2header.h g2parameters.h g2main.cpp
	$(COMPILER) -c $(FLAGS)  g2main.cpp

g2spectra.o: g2header.h g2parameters.h g2spectra.cpp
	$(COMPILER) -c $(FLAGS)  g2spectra.cpp

g2power.o: g2header.h g2parameters.h g2power.cpp
	$(COMPILER) -c $(FLAGS)  g2power.cpp

g2model.o: g2header.h g2parameters.h g2model.cpp
	$(COMPILER) -c $(FLAGS)  g2model.cpp

g2functions.o: g2header.h g2parameters.h g2functions.cpp
	$(COMPILER) -c $(FLAGS)  g2functions.cpp

g2output.o: g2header.h g2parameters.h g2output.cpp
	$(COMPILER) -c $(FLAGS)  g2output.cpp

g2init.o: g2header.h g2parameters.h g2init.cpp
	$(COMPILER) -c $(FLAGS)  g2init.cpp

oclean: 
	rm -f *.o

fclean:
	rm -f gabe *.dat *.txt
	rm -f ./slices/*.dat
	rm -f ./energy/*.dat
	rm -f ./gravrhoslices/*.dat
	rm -f ./source/*.dat
		
dirmake:
	@ if test -d slices; then echo directory "'slices'" exists; else mkdir slices; echo made directory "'slices'"; fi

	@ if test -d energy; then echo directory “’energy’” exists; else mkdir energy; echo made directory “’energy’”; fi

	@ if test -d gravrhoslices; then echo directory “’gravrhoslices’” exists; else mkdir gravrhoslices; echo made directory “’gravrhoslices’”;fi

	@ if test -d source; then echo directory “source” exists; else mkdir source; echo made directory “source”;fi	


run: clean compile
	./gabe

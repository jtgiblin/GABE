COMPILER = em++
OFLAG = -O3 -ffast-math -flto
FNS = -s EXPORTED_FUNCTIONS="['_copyout_fld','_sim_init','_sim_step']"
FLAGS = --std=c++11 -s NO_EXIT_RUNTIME=1 $(OFLAG) $(FNS)
LINKS = -L/opt/local/lib
OUTFILE = gabem.js

do: oclean compile

xcmake : do
	
compile: g2header.h g2parameters.h g2model.o g2functions.o g2main.o g2init.o 
	$(COMPILER) $(FLAGS) g2model.o g2functions.o g2init.o g2main.o $(LINKS) -o $(OUTFILE)

clean : oclean

g2main.o: g2header.h g2parameters.h g2main.cpp
	$(COMPILER) -c $(FLAGS)  g2main.cpp

g2model.o: g2header.h g2parameters.h g2model.cpp
	$(COMPILER) -c $(FLAGS)  g2model.cpp

g2functions.o: g2header.h g2parameters.h g2functions.cpp
	$(COMPILER) -c $(FLAGS)  g2functions.cpp

g2init.o: g2header.h g2parameters.h g2init.cpp
	$(COMPILER) -c $(FLAGS)  g2init.cpp

oclean: 
	rm -f *.o *.js.mem

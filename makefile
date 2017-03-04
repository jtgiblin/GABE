COMPILER = em++
OFLAG = -ffast-math -O3 -flto
FLAGS = --std=c++11 --bind $(OFLAG)
LINKS = -L/opt/local/lib

do: oclean compile 

xcmake : do
	
compile: g2header.h g2parameters.h g2model.o g2functions.o g2main.o g2init.o 
	$(COMPILER) $(FLAGS) g2model.o g2functions.o g2init.o g2main.o $(LINKS) -o gabe.js

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
	rm -f *.o *.js *.js.mem

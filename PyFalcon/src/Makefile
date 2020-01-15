OPT += 
#OPT += -DDEBUG

CXX = gcc
INCL += -I/home/diptajym/physics2/dependencies/hdf5_new_install/include -DH5_USE_16_API
INCL += -L/home/diptajym/physics2/dependencies/hdf5_new_install/lib -lhdf5 -lz

all:
	@make falcon

falcon: start.c
	$(CXX) -fPIC -shared integrator.c io.c force.c initialize.c $(INCL) $(OPT) $? -o libfalcon.so

clean:
	$(RM) ./*.o ./falcon

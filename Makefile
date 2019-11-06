#OPT += -DDEBUG

CXX = gcc
INCL += -I/home/diptajym/physics2/dependencies/hdf5_new_install/include -DH5_USE_16_API
INCL += -L/home/diptajym/physics2/dependencies/hdf5_new_install/lib -lhdf5 -lz

all:
	@make falcon

falcon: main.c
	$(CXX) integrator.c io.c force.c initialize.c $(INCL) $? -o $@

clean:
	$(RM) ./*.o ./falcon

OBJECTS1 = Multi_dielectric.o 
OBJECTS2 = Operate_four.o edgenod.o
OBJECTS3 = labelinv.o hpsort.o comi.o com.o
OBJECTS = $(OBJECTS1) $(OBJECTS2) $(OBJECTS3)

FFLAGS = -O2

ansys_heyj2.x: $(OBJECTS)
	ifort $(FFLAGS) $(OBJECTS) -o ./readcdb.x

.f.o:
	ifort $(FFLAGS) -c $< $(VERBOSE) 

.c.o:
	icc $(FFLAGS) -c $< $(VERBOSE)
.PHONY:clean
clean:
	rm *.o *.x

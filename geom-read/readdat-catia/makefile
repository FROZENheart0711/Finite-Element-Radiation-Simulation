ALL: readdat-CATIA-new.x #readdat-CATIA-2017.x readdat-CATIA-2018.x
OBJECTS1 = Multi_dielectric.o
OBJECTS2 = Operate_four.o edgenod.o
OBJECTS3 = labelinv.o hpsort.o comi.o com.o
OBJECTS4 = Multi_dielectric-2017.o
OBJECTS5 = Multi_dielectric-ipat.o
OBJECTS = $(OBJECTS1) $(OBJECTS2) $(OBJECTS3) $(OBJECTS4) $(OBJECTS5)
OBJECTS2015 = $(OBJECTS1) $(OBJECTS2) $(OBJECTS3)
OBJECTS2017 = $(OBJECTS4) $(OBJECTS2) $(OBJECTS3)
OBJECTS2018 = $(OBJECTS5) $(OBJECTS2) $(OBJECTS3)

FFLAGS = -O1

#readdat-CATIA-2017.x: $(OBJECTS)
#	ifort $(FFLAGS) $(OBJECTS2017) -o ./readdat-CATIA-2017.x
readdat-CATIA-new.x: $(OBJECTS)
	ifort $(FFLAGS) $(OBJECTS2015) -o ./readdat-catia.x
#readdat-CATIA-2018.x: $(OBJECTS)
#	ifort $(FFLAGS) $(OBJECTS2018) -o ./readdat-CATIA-2018.x
.f.o:
	ifort -c -O1 $< $(VERBOSE) 

.c.o:
	icc -c -O1 $< $(VERBOSE)
.PHONY:clean
clean:
	rm *.o *.x

OBJECTS1 = ansys_heyj.o com.o comi.o dis2o.o dis.o edgenod.o fcnd.o fdposi.o fstrcat.o hpsort.o  idsame.o  
OBJECTS2 = labelfemm.o labelfemi.o labelfemo.o labeltran.o labelinv.o nodeface.o num2char.o iof.o findpecedge.o
OBJECTS3 = readansys.o sur_edge_label.o sur_vol_coor.o sur_vol_edge.o vol_edge_label.o write_inform.o
OBJECTS = $(OBJECTS1) $(OBJECTS2) $(OBJECTS3)

FFLAGS = -O1

ansys_heyj.x: $(OBJECTS)
	ifort $(FFLAGS) $(OBJECTS) -o ./ansys_heyj.x

.f.o:
	ifort $(FFLAGS) -c $< $(VERBOSE) 

.c.o:
	icc -c $< $(VERBOSE)

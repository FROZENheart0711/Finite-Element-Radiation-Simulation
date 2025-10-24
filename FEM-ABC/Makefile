CC      = icc
FC      = ifort
CL      = icc

OPTF    = -O2 -qopenmp
OPTC	= -O2 -qopenmp -Wall 
ARCH    = ar r
RANLIB  = ranlib

CDEFS   =-DADD_

INC_DIR = -I./ -I/opt/intel/composerxe-2011.0.084/mkl/include

FFLAGS  = ${OPTF}
CXXFLAGS= ${OPTC} $(CDEFS)

#LIBBLAS=-L/opt/intel/composerxe-2011.0.084/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm -openmp
LIBBLAS=-L/opt/intel/composerxe-2011.0.084/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -qopenmp

#LIBBLAS=-L/opt/intel/composerxe-2011.0.084/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lpthread #-lm -openmp
INTELLIB= -lifcore -limf


MUMPSLIB= /home/bjxia/software/MUMPS_4.10.021/lib/libcmumps.a /home/bjxia/software/MUMPS_4.10.021/libseq/libmpiseq.a \
         /home/bjxia/software/MUMPS_4.10.021/lib/libzmumps.a  /home/bjxia/software/MUMPS_4.10.021/lib/libmumps_common.a

LPORDDIR = /home/bjxia/software/MUMPS_4.10.021/PORD/lib/
IPORD    = -I/home/bjxia/software/MUMPS_4.10.021/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord

LORDERINGS =  $(LPORD)
IORDERINGS =  $(IPORD)

LIB     = $(LIBBLAS) $(MUMPSLIB) $(LORDERINGS) $(INTELLIB) #-lstdc++

INCMUMPS       = -I/home/bjxia/software/MUMPS_4.10.021/include
INCMUMPSSEQ    = -I/home/bjxia/software/MUMPS_4.10.021/libseq

INC     = $(INC_DIR) $(IORDERINGS) $(INCMUMPS) $(INCMUMPSSEQ)

FOBJ= addsource.o cdotc.o detm.o elesur_md.o farfld.o getfaceinfo.o labeltran.o signcheck.o \
      assemfem.o cdotmul.o dhpsort.o elevol1.o fem.o getsubmatrix_s.o matset.o read_fem_inp.o sign.o \
      assemsur_bd.o cmatin.o dotmul.o elevol_N.o find_face.o getsubmatrix_v.o opfile.o read_geom_fem.o sparse2csr.o \
	  rhside.o elesursrcm.o cal_et.o cal_mode.o elesrc.o assemsrc.o devaCurlBsur.o devagradNsur.o devaBsur.o devaNsur.o \
	  invJmatrixf.o areapatch.o gausspoints.o numint.o s11_zinvol.o

COBJ= main.o basic.o fstrcat.o

Solve.X : $(FOBJ) $(COBJ) $(LAPACKOBJ)
	$(CL) $(CXXFLAGS) $(FOBJ) $(COBJ) $(LIB) -o $@

.c.o:
	$(CC) $(CXXFLAGS) $(INC) $(INCSEQ) -c $< $(VERBOSE) 
.cpp.o:
	$(CC) $(CXXFLAGS) $(INC) $(INCSEQ) -c $< $(VERBOSE) 
.f.o:
	$(FC) $(FFLAGS) $(INC) $(INCSEQ) -c $< $(VERBOSE) 

.PHONY: clean
clean :
	rm -f *.o *.X *.a

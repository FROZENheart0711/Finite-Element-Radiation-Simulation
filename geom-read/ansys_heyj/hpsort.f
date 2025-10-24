 
	SUBROUTINE HPSORT(N,RA,iqof,iqfo)
	INTEGER N
	real RA(N)
	INTEGER I,IR,J,L,iqof(n),iqfo(n)
c	double precision RRA
	do 5 i=1,n
	iqof(i)=i
5	iqfo(i)=i
	IF (N.LT.2) RETURN
	L=N/2+1
	IR=N
10	CONTINUE
	IF (L.GT.1) THEN
	L=L-1
	RRA=RA(L)
	ELSE
	RRA=RA(IR)
	RA(IR)=RA(1)
	it1=iqfo(1)
	it2=iqfo(ir)
	iqfo(1)=it2
	iqfo(ir)=it1
	iqof(it2)=1
	iqof(it1)=ir
	IR=IR-1
	IF (IR.EQ.1) THEN
	RA(1)=RRA
	RETURN
	ENDIF
	ENDIF
	I=L
	J=L+L
	it1=iqfo(i)
20	IF (J.LE.IR) THEN
	   IF (J.LT.IR) THEN
	     IF (RA(J).LT.RA(J+1)) J=J+1
	     ENDIF
	     IF (RRA.LT.RA(J)) THEN
	     RA(I)=RA(J)
	     it2=iqfo(j)
	     iqof(it2)=i
	     iqfo(i)=it2
	     I=J
	     J=J+J
	     ELSE
	     J=IR+1
	     ENDIF
	GOTO 20
	ENDIF
	RA(I)=RRA
	     iqfo(i)=it1
	     iqof(it1)=i
	GOTO 10
	END	      

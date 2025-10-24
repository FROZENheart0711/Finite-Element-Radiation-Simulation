        SUBROUTINE CMATIN(NM,A,N,B,M,DET,IPVOT,INDEX,PIVOT)

C.......THIS PROGRAM IS FOR COMPLX MATRIX INVERSION AND SIMULT LINER EQS.
        DIMENSION IPVOT(NM),INDEX(NM,2)
        COMPLEX A(NM,N),B(NM,N),DET,T,PIVOT(NM)
        EQUIVALENCE(IROW,JROW),(ICOL,JCOL)
        IROW=0
        JROW=0
        ICOL=0
        JCOL=0
57      DET=(1.0,0.0)
        DO 17 J=1, N
17      IPVOT(J)=0
        DO 135 I=1, N
        T=(0.0,0.0)
        DO 9 J=1, N
        IF(IPVOT(J)-1) 13,9,13
13      DO 23 K=1, N
        IF(IPVOT(K)-1) 43,23,81
43      G1=CABS(T)
        G2=CABS(A(J,K))
        IF(G1-G2) 83,23,23
83      IROW=J
        ICOL=K
        T=A(J,K)
23      CONTINUE
9       CONTINUE
        IPVOT(ICOL)=IPVOT(ICOL)+1
        IF(IROW-ICOL) 73,109,73
73      DET=-DET
        DO 12 L=1, N
        T=A(IROW,L)
        A(IROW,L)=A(ICOL,L)
12      A(ICOL,L)=T
        IF(M) 109,109,33
33      DO 2 L=1, M
        T=B(IROW,L)
        B(IROW,L)=B(ICOL,L)
2       B(ICOL,L)=T
109     INDEX(I,1)=IROW
        INDEX(I,2)=ICOL
        PIVOT(I)=A(ICOL,ICOL)
        IF(CABS(DET).GT.1.0E+10) DET=DET/CABS(DET)
        IF(CABS(PIVOT(I)).GT.1.0E+10) DET=DET/CABS(PIVOT(I))
C        DET=DET*PIVOT(I)
        A(ICOL,ICOL)=(1.0,0.0)
        DO 505 L=1, N
505     A(ICOL,L)=A(ICOL,L)/PIVOT(I)
        IF(M) 347,347,66
66      DO 52 L=1, M
52      B(ICOL,L)=B(ICOL,L)/PIVOT(I)
347     DO 135 LI=1, N
        IF(LI-ICOL) 21,135,21
21      T=A(LI,ICOL)
        A(LI,ICOL)=(0.0,0.0)
        DO 89 L=1, N
89      A(LI,L)=A(LI,L)-A(ICOL,L)*T
        IF(M) 135,135,18
18      DO 68 L=1, M
68      B(LI,L)=B(LI,L)-B(ICOL,L)*T
135     CONTINUE
222     DO 3 I=1, N
        L=N-I+1
        IF(INDEX(L,1)-INDEX(L,2)) 19,3,19
19      JROW=INDEX(L,1)
        JCOL=INDEX(L,2)
        DO 549 K=1, N
        T=A(K,JROW)
        A(K,JROW)=A(K,JCOL)
        A(K,JCOL)=T
549     CONTINUE
3       CONTINUE
81      RETURN
        END

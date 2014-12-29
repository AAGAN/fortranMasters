      SUBROUTINE SORTRX(N,DATA1,INDEX1)

      INTEGER   N,INDEX1(N)
      DOUBLEPRECISION      DATA1(N)
 
      INTEGER   LSTK(31),RSTK(31),ISTK
      INTEGER   L,R,I,J,P,INDEXP,INDEXT
      DOUBLEPRECISION      DATAP
 
      INTEGER   M
      PARAMETER (M=9)
 
      DO 50 I=1,N
         INDEX1(I)=I
   50    CONTINUE
 
      IF (N.LE.M) GOTO 900
 
      ISTK=0
      L=1
      R=N
 
  200 CONTINUE
 
      I=L
      J=R
 
      P=(L+R)/2
      INDEXP=INDEX1(P)
      DATAP=DATA1(INDEXP)
 
      IF (DATA1(INDEX1(L)) .GT. DATAP) THEN
         INDEX1(P)=INDEX1(L)
         INDEX1(L)=INDEXP
         INDEXP=INDEX1(P)
         DATAP=DATA1(INDEXP)
      ENDIF
 
      IF (DATAP .GT. DATA1(INDEX1(R))) THEN
         IF (DATA1(INDEX1(L)) .GT. DATA1(INDEX1(R))) THEN
            INDEX1(P)=INDEX1(L)
            INDEX1(L)=INDEX1(R)
         ELSE
            INDEX1(P)=INDEX1(R)
         ENDIF
         INDEX1(R)=INDEXP
         INDEXP=INDEX1(P)
         DATAP=DATA1(INDEXP)
      ENDIF
 
  300 CONTINUE
 
         I=I+1
         IF (DATA1(INDEX1(I)).LT.DATAP) GOTO 300
 
  400 CONTINUE
 
         J=J-1
         IF (DATA1(INDEX1(J)).GT.DATAP) GOTO 400
 
      IF (I.LT.J) THEN
 
         INDEXT=INDEX1(I)
         INDEX1(I)=INDEX1(J)
         INDEX1(J)=INDEXT
         GOTO 300
      ELSE
 
         IF (R-J .GE. I-L .AND. I-L .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=J+1
            RSTK(ISTK)=R
            R=I-1
         ELSE IF (I-L .GT. R-J .AND. R-J .GT. M) THEN
            ISTK=ISTK+1
            LSTK(ISTK)=L
            RSTK(ISTK)=I-1
            L=J+1
         ELSE IF (R-J .GT. M) THEN
            L=J+1
         ELSE IF (I-L .GT. M) THEN
            R=I-1
         ELSE
            IF (ISTK.LT.1) GOTO 900
            L=LSTK(ISTK)
            R=RSTK(ISTK)
            ISTK=ISTK-1
         ENDIF
         GOTO 200
      ENDIF
 
  900 CONTINUE
 
      DO 950 I=2,N
         IF (DATA1(INDEX1(I-1)) .GT. DATA1(INDEX1(I))) THEN
            INDEXP=INDEX1(I)
            DATAP=DATA1(INDEXP)
            P=I-1
  920       CONTINUE
               INDEX1(P+1) = INDEX1(P)
               P=P-1
               IF (P.GT.0) THEN
                  IF (DATA1(INDEX1(P)).GT.DATAP) GOTO 920
               ENDIF
            INDEX1(P+1) = INDEXP
         ENDIF
  950    CONTINUE
 
      END
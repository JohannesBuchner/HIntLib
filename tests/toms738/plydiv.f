      SUBROUTINE PLYDIV (PA, PB, PQ, PR)
C
C   This version :  12 December 1991
C
C
C   ------------------------------------------------------------
C
C   The following COMMON block, used by many subroutines,
C   gives the order Q of a field, its characteristic P, and its
C   addition, multiplication, and subtraction tables.
C   The parameter MAXQ gives the order of the largest field to
C   be handled.
C
      INTEGER MAXQ
      PARAMETER (MAXQ=50)
 
      INTEGER P, Q, ADD(0:MAXQ-1,0:MAXQ-1)
      INTEGER MUL(0:MAXQ-1, 0:MAXQ-1), SUB(0:MAXQ-1, 0:MAXQ-1)
      COMMON /FIELD/ P, Q, ADD, MUL, SUB
      SAVE /FIELD/
C
C   The following definitions concern the representation of
C   polynomials.
C
      INTEGER MAXDEG, DEG
      PARAMETER (MAXDEG=50, DEG=-1)
C
C   The parameter MAXDEG gives the highest degree of polynomial
C   to be handled.  Polynomials stored as arrays have the
C   coefficient of degree n in POLY(N), and the degree of the
C   polynomial in POLY(-1).  The parameter DEG is just to remind
C   us of this last fact.  A polynomial which is identically 0
C   is given degree -1.
C
C   A polynomial can also be stored in an integer I, with
C        I = AN*Q**N + ... + A0.
C   Routines ITOP and PTOI convert between these two formats.
C
C   -----------------------------------------------------------
C
C
C
      INTEGER I, J, D, M, BINV, DEGB, DEGR, DEGQ
      INTEGER PA(-1:MAXDEG), PB(-1:MAXDEG), PQ(-1:MAXDEG)
      INTEGER PR(-1:MAXDEG)
      INTEGER PTQ(-1:MAXDEG), PTR(-1:MAXDEG)
C
C   Divides polynomial PA by PB, putting the quotient in PQ
C   and the remainder in PR.
C   Coefficients are elements of the field of order Q.
C
      IF (PB(DEG) .EQ. -1) THEN
        WRITE (*,*) ' PLYDIV :  Divide by zero'
        STOP
      ENDIF
C
      DO 10 I = -1, MAXDEG
        PTQ(I) = 0
   10   PTR(I) = PA(I)
      DEGR = PA(DEG)
      DEGB = PB(DEG)
      DEGQ = DEGR - DEGB
      IF (DEGQ .LT. 0) DEGQ = -1
C
C   Find the inverse of the leading coefficient of PB.
C
      J = PB(DEGB)
      DO 15 I = 1, P-1
        IF (MUL(I,J) .EQ. 1) BINV = I
   15 CONTINUE
C
      DO 30 D = DEGQ, 0, -1
        M = MUL (PTR(DEGR), BINV)
        DO 20 I = DEGB, 0, -1
   20     PTR(DEGR+I-DEGB) = SUB (PTR(DEGR+I-DEGB), MUL (M, PB(I)))
        DEGR = DEGR - 1
   30   PTQ(D) = M
C
      DO 40 I = 0, MAXDEG
        PQ(I) = PTQ(I)
   40   PR(I) = PTR(I)
      PQ(DEG) = DEGQ
   50 IF (PR(DEGR) .EQ. 0 .AND. DEGR .GE. 0) THEN
        DEGR = DEGR - 1
        GOTO 50
      ENDIF
      PR(DEG) = DEGR
C
      RETURN
      END
C
C     *****  end of SUBROUTINE PLYDIV

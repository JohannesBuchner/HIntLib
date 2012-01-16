      SUBROUTINE ITOP (IN, POLY)
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
      INTEGER IN, I, J, POLY(-1:MAXDEG)
C
C   Converts an integer to a polynomial with coefficients in the
C   field of order Q.
C
      DO 10 J = -1, MAXDEG
   10   POLY(J) = 0
C
      I = IN
      J = -1
   20 IF (I .GT. 0) THEN
        J = J + 1
        IF (J .GT. MAXDEG) THEN
          WRITE (*,*) ' ITOP :  Polynomial exceeds MAXDEG'
          STOP
        ENDIF
        POLY(J) = MOD (I, Q)
        I = I / Q
        GOTO 20
      ENDIF
      POLY(DEG) = J
      RETURN
      END
C
C     *****  end of SUBROUTINE ITOP

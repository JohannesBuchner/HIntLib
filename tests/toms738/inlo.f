      SUBROUTINE INLO (DIM, BASE, SKIP)
C
C   This version :  12 February 1992
C
C   See the general comments on implementing Niederreiter's
C   low-discrepancy sequences.
C
C   This subroutine calculates the values of Niederreiter's
C   C(I,J,R) and performs other initialization necessary
C   before calling GOLO.
C
C INPUT :
C   DIMEN - The dimension of the sequence to be generated.
C      {In the argument of INLO, DIMEN is called DIM because
C      DIMEN is subsequently passed via COMMON and is called
C      DIMEN there.}
C
C   BASE  - The prime or prime-power base to be used.
C   SKIP  - The number of values to throw away at the beginning
C           of the sequence.
C
C OUTPUT :
C   To GOLO, labelled common /COMM/.
C
C USES :
C   Calls CALCC to calculate the values of CJ.
C   Calls SETFLD to set up field arithmetic tables.
C   Calls CHARAC to check that base is a prime or prime-power
C     in the range we can handle.
C
C   -------------------------------------------------------------
C
C
C   This segment defines the common block /COMM/ and some
C   associated parameters.  These are for use in the general base
C   version of the generator.
C
      INTEGER MAXDIM, MAXFIG, NBITS
      PARAMETER (MAXDIM=12, MAXFIG=20, NBITS=31)
C
C   The parameter MAXDIM is the maximum dimension that will be used.
C   MAXFIG is the maximum number of base-Q digits we can handle.
C   MAXINT is the largest fixed point integer we can represent.
C   NBITS is the number of bits in a fixed-point integer, not
C     counting the sign.
C     ***** NBITS is machine dependent *****
C
      INTEGER C(MAXDIM, MAXFIG, 0:MAXFIG-1)
      INTEGER COUNT(0:MAXFIG-1), D(MAXDIM, MAXFIG)
      INTEGER NEXTQ(MAXDIM), QPOW(MAXFIG)
      INTEGER DIMEN, NFIGS
      REAL    RECIP
      COMMON  /COMM/ C, COUNT, D, NEXTQ, QPOW, DIMEN, NFIGS, RECIP
      SAVE    /COMM/
C
C   The common block /COMM/ :
C     C     - Contains the values of Niederreiter's C(I,J,R)
C     COUNT - The index of the current item in the sequence,
C             expressed as an array of base-Q digits.  COUNT(R)
C             is the same as Niederreiter's AR(N) (page 54)
C             except that N is implicit.
C     D     - The values of D(I,J).
C     NEXTQ - The numerators of the next item in the series.  These
C             are like Niederreiter's XI(N) (page 54) except that
C             N is implicit, and the NEXTQ are integers.  To obtain
C             the values of XI(N), multiply by RECIP.
C     QPOW  - To speed things up a bit. QPOW(I) = Q ** (NFIGS-I).
C     DIMEN   - The dimension of the sequence to be generated
C     NFIGS - The number of base-Q digits we are using.
C     RECIP - 1.0 / (Q ** NFIGS)
C
C   Array C of the common block is set up by subroutine CALCC.
C   The other items in the common block are set up by INLO.
C
C   ------------------------------------------------------------
C
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
      INTEGER I, J, NQ, R, DIM, SKIP , BASE, CHARAC
      REAL    TEMP
C
      DIMEN = DIM
C
C     This assignment just relabels the variable
C     for subsequent use.
C
      IF (DIMEN.LE.0 .OR. DIMEN.GT.MAXDIM) THEN
        WRITE (*,*) ' INLO :  Bad dimension'
        STOP
      ENDIF
C
C
      IF (CHARAC(BASE) .EQ. 0) THEN
        WRITE (*,*) ' INLO : Base not prime power or out of range'
        STOP
      ENDIF
C
      CALL SETFLD (BASE)
C
C   Calculate how many figures to use in base Q = BASE
C
      TEMP = LOG(2.0 ** NBITS - 1)/LOG(REAL(Q))
      NFIGS = MIN(MAXFIG, INT(TEMP))
C
      CALL CALCC
C
C   Set up RECIP and QPOW(I) = Q ** (NFIGS-I)
C
      RECIP = 1.0 / (Q ** NFIGS)
      QPOW(NFIGS) = 1
      DO 5 I = NFIGS-1, 1, -1
    5   QPOW(I) = Q * QPOW(I+1)
C
C   Initialize COUNT
C
      I = SKIP
      DO 10 R = 0, NFIGS-1
        COUNT(R) = MOD(I, Q)
        I = I / Q
   10 CONTINUE
      IF (I .NE. 0) THEN
        WRITE (*,*) ' INLO :  SKIP too long'
        STOP
      ENDIF
C
C   Initialize D
C
      DO 20 I = 1, DIMEN
        DO 20 J = 1, NFIGS
   20     D(I,J) = 0
C
      DO 50 R = 0, NFIGS-1
        IF (COUNT(R) .NE. 0) THEN
          DO 40 I = 1, DIMEN
            DO 30 J = 1, NFIGS
              D(I,J) = ADD (D(I,J), MUL (C(I,J,R), COUNT(R)))
   30       CONTINUE
   40     CONTINUE
        ENDIF
   50 CONTINUE
C
C   Initialize NEXTQ
C
      DO 70 I = 1, DIMEN
        NQ = 0
        DO 60 J = 1, NFIGS
          NQ = NQ + D(I,J) * QPOW(J)
   60   CONTINUE
        NEXTQ(I) = NQ
   70 CONTINUE
C
      RETURN
      END
C
C     *****  end of SUBROUTINE INLO

      SUBROUTINE GOLO (QUASI)
      REAL QUASI(*)
C
C This version : 21 February 1992
C
C See the general comments on implementing Niederreiter's
C low-discrepancy sequences.
C
C Call subroutine INLO before calling GOLO.  Thereafter
C GOLO generates a new quasi-random vector on each call,
C
C INPUT
C   From INLO, labelled common /COMM/ and labelled common
C     /FIELD/, properly initialized.  (INLO calls SETFLD
C     to initialize /FIELD/).
C
C OUTPUT
C   To the user's program, the next vector in the sequence in
C     array QUASI.
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
C   addition, multiplication and subtraction tables.
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
      INTEGER I, J, R, OLDCNT, DIFF, NQ
C
C Multiply the numerators in NEXTQ by RECIP to get the next
C quasi-random vector
C
      DO 5 I = 1, DIMEN
        QUASI(I) = NEXTQ(I) * RECIP
    5 CONTINUE
C
C Update COUNT, treated as a base-Q integer.  Instead of
C recalculating the values of D from scratch, we update
C them for each digit of COUNT which changes.  In terms of
C Niederreiter page 54, NEXTQ(I) corresponds to XI(N), with
C N being implicit, and D(I,J) corresponds to XI(N,J), again
C with N implicit.  Finally COUNT(R) corresponds to AR(N).
C
      R = 0
   10 IF (R .EQ. NFIGS) THEN
        WRITE (*,*) ' Too many calls on subroutine GOLO'
        STOP
      ENDIF
      OLDCNT = COUNT(R)
      IF (COUNT(R) .LT. Q-1) THEN
        COUNT(R) = COUNT(R) + 1
      ELSE
        COUNT(R) = 0
      ENDIF
      DIFF = SUB(COUNT(R), OLDCNT)
C
C Digit R has just changed.  DIFF says how much it changed
C by.  We use this to update the values of array D.
C
      DO 40 I = 1, DIMEN
        DO 30 J = 1, NFIGS
          D(I,J) = ADD(D(I,J), MUL(C(I,J,R), DIFF))
   30   CONTINUE
   40 CONTINUE
C
C If COUNT(R) is now zero, we must propagate the carry
C
      IF (COUNT(R).EQ.0) THEN
        R = R + 1
        GOTO 10
      ENDIF
C
C Now use the updated values of D to calculate NEXTQ.
C Array QPOW helps to speed things up a little :
C   QPOW(J) is Q ** (NFIGS-J).
C
      DO 60 I = 1, DIMEN
        NQ = 0
        DO 50 J = 1, NFIGS
          NQ = NQ + D(I,J) * QPOW(J)
   50   CONTINUE
        NEXTQ(I) = NQ
   60 CONTINUE
C
      RETURN
      END
C
C     ***** end of SUBROUTINE GOLO

      SUBROUTINE CALCC2
C
C   This version :  12 February 1992
C
C      *****  For base-2 only.
C
C   See the general comments on implementing Niederreiter's
C   low-discrepancy sequences.
C
C   This program calculates the values of the constants C(I,J,R).
C   As far as possible, we use Niederreiter's notation.
C   For each value of I, we first calculate all the corresponding
C   values of C :  these are held in the array CI.  All these
C   values are either 0 or 1.  Next we pack the values into the
C   array CJ, in such a way that CJ(I,R) holds the values of C
C   for the indicated values of I and R and for every value of
C   J from 1 to NBITS.  The most significant bit of CJ(I,R)
C   (not counting the sign bit) is C(I,1,R) and the least
C   significant bit is C(I,NBITS,R).
C     When all the values of CJ have been calculated, we return
C   this array to the calling program.
C
C  --------------------------------------------------------------
C
C   We define the common block /COMM2/ and some
C   associated parameters.  These are for use in the base 2
C   version of the generator.
C
      INTEGER MAXDIM, NBITS
      PARAMETER (MAXDIM=12, NBITS=31)
C
C   The parameter MAXDIM is the maximum dimension that will be used.
C   NBITS is the number of bits (not counting the sign) in a
C   fixed-point integer.
C
      INTEGER CJ(MAXDIM, 0:NBITS-1), DIMEN, COUNT, NEXTQ(MAXDIM)
      COMMON /COMM2/ CJ, DIMEN, COUNT, NEXTQ
      SAVE   /COMM2/
C
C   The common block /COMM2/ :
C     CJ    - Contains the packed values of Niederreiter's C(I,J,R)
C     DIMEN   - The dimension of the sequence to be generated
C     COUNT - The index of the current item in the sequence,
C             expressed as an array of bits.  COUNT(R) is the same
C             as Niederreiter's AR(N) (page 54) except that N is
C             implicit.
C     NEXTQ - The numerators of the next item in the series.  These
C             are like Niederreiter's XI(N) (page 54) except that
C             N is implicit, and the NEXTQ are integers.  To obtain
C             the values of XI(N), multiply by RECIP (see GOLO2).
C
C   Array CJ of the common block is set up by subroutine CALCC2.
C   The other items in the common block are set up by INLO2.
C
C   --------------------------------------------------------------
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
C   ---------------------------------------------------------------
C
C
C
C   MAXE   - We need MAXDIM irreducible polynomials over Z2.
C            MAXE is the highest degree among these.
C   MAXV   - The maximum possible index used in V.
C
      INTEGER MAXE, MAXV
      PARAMETER (MAXE=5, MAXV=NBITS+MAXE)
C
C INPUT :
C   The array CJ to be initialised, and DIMEN the number of
C   dimensions we are using, are transmitted through /COMM2/.
C
C OUTPUT :
C   The array CJ is returned to the calling program.
C
C USES :
C   Subroutine SETFLD is used to set up field arithmetic tables.
C   (Although this is a little heavy-handed for the field of
C   order 2, it makes for uniformity with the general program.)
C   Subroutine CALCV is used to for the auxiliary calculation
C   of the values V needed to get the Cs.
C
      INTEGER PX(-1:MAXDEG), B(-1:MAXDEG)
      INTEGER V(0:MAXV), CI(NBITS, 0:NBITS-1)
      INTEGER E, I, J, R, U, TERM
C
      INTEGER IRRED(MAXDIM, -1:MAXE)
      SAVE IRRED
C
C   This DATA statement supplies the coefficients and the
C   degrees of the first 12 irreducible polynomials over Z2.
C   They are taken from Lidl and Niederreiter, FINITE FIELDS,
C   Cambridge University Press (1984), page 553.
C   The order of these polynomials is the same as the order in
C   file 'irrtabs.dat' used by the general program.
C
C   In this block PX(I, -1) is the degree of the Ith polynomial,
C   and PX(I, N) is the coefficient of x**n in the Ith polynomial.
C
      DATA (IRRED(1,I), I=-1,1)  / 1,0,1 /
      DATA (IRRED(2,I), I=-1,1)  / 1,1,1 /
      DATA (IRRED(3,I), I=-1,2)  / 2,1,1,1 /
      DATA (IRRED(4,I), I=-1,3)  / 3,1,1,0,1 /
      DATA (IRRED(5,I), I=-1,3)  / 3,1,0,1,1 /
      DATA (IRRED(6,I), I=-1,4)  / 4,1,1,0,0,1 /
      DATA (IRRED(7,I), I=-1,4)  / 4,1,0,0,1,1 /
      DATA (IRRED(8,I), I=-1,4)  / 4,1,1,1,1,1 /
      DATA (IRRED(9,I), I=-1,5)  / 5,1,0,1,0,0,1 /
      DATA (IRRED(10,I), I=-1,5) / 5,1,0,0,1,0,1 /
      DATA (IRRED(11,I), I=-1,5) / 5,1,1,1,1,0,1 /
      DATA (IRRED(12,I), I=-1,5) / 5,1,1,1,0,1,1 /
C
C   Prepare to work in Z2
C
      CALL SETFLD (2)
C
      DO 1000 I = 1, DIMEN
C
C   For each dimension, we need to calculate powers of an
C   appropriate irreducible polynomial :  see Niederreiter
C   page 65, just below equation (19).
C     Copy the appropriate irreducible polynomial into PX,
C   and its degree into E.  Set polynomial B = PX ** 0 = 1.
C   M is the degree of B.  Subsequently B will hold higher
C   powers of PX.
C
        E = IRRED(I, DEG)
        DO 10 J = -1, E
          PX(J) = IRRED(I,J)
   10   CONTINUE
        B(DEG) = 0
        B(0) = 1
C
C   Niederreiter (page 56, after equation (7), defines two
C   variables Q and U.  We do not need Q explicitly, but we
C   do need U.
C
        U = 0
C
        DO 90 J = 1, NBITS
C
C   If U = 0, we need to set B to the next power of PX
C   and recalculate V.  This is done by subroutine CALCV.
C
          IF (U .EQ. 0) CALL CALCV (PX, B, V, MAXV)
C
C Now C is obtained from V.  Niederreiter
C obtains A from V (page 65, near the bottom), and then gets
C C from A (page 56, equation (7)).  However this can be done
C in one step.  Here CI(J,R) corresponds to
C Niederreiter's C(I,J,R).
C
          DO 50 R = 0, NBITS-1
            CI(J,R) = V(R+U)
   50     CONTINUE
C
C Increment U.  If U = E, then U = 0 and in Niederreiter's
C paper Q = Q + 1.  Here, however, Q is not used explicitly.
C
          U = U + 1
          IF (U .EQ. E) U = 0
  90    CONTINUE
C
C  The array CI now holds the values of C(I,J,R) for this value
C  of I.  We pack them into array CJ so that CJ(I,R) holds all
C  the values of C(I,J,R) for J from 1 to NBITS.
C
        DO 120 R = 0, NBITS-1
          TERM = 0
          DO 110 J = 1, NBITS
            TERM = 2 * TERM + CI(J,R)
  110     CONTINUE
          CJ(I,R) = TERM
  120   CONTINUE
C
 1000 CONTINUE
      END
C
C      *****  end of CALCC2

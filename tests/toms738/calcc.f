      SUBROUTINE CALCC
C
C    This version : 12 February 1992
C
C    See the general comments on implementing Niederreiter's
C    low-discrepancy sequences.
C
C    This routine calculates the values of the constants C(I,J,R).
C    As far as possible, we use Niederreiter's notation.
C    We calculate the values of C for each I in turn.
C    When all the values of C have been calculated, we return
C    this array to the calling program.
C
C    Irreducible polynomials are read from file 'irrtabs.dat'
C    through Fortran unit 1.  These polys must have been put on
C    the file beforehand by GFPOLYS.  Unit 1 is closed before
C    entry to CALCC and after returning from the subroutine.
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
C
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
C
C
C   MAXE   - We need MAXDIM irreducible polynomials over GF(Q).
C            MAXE is the highest degree among these.
C   MAXV   - The maximum index used in V.
C   GFUNIT - The unit number to read field tables.
C   NPOLS  - The number of precalculated irreducible polys.
C
      INTEGER MAXE, MAXV, GFUNIT, NPOLS
      PARAMETER (MAXE=5,  GFUNIT=1, NPOLS=25)
      PARAMETER (MAXV = MAXFIG + MAXE)
C
C INPUT :
C   DIMEN, the number of dimensions in use, and NFIGS, the number
C   of base-Q figures to be used, are passed in through common COMM.
C   Necessary field arithmetic tables are passed through common
C   FIELD.
C
C OUTPUT
C   The array C is returned to the calling program.
C
C USES
C   Subroutine CALCV is used for the auxiliary calculation
C   of the values V needed to get the Cs.
C
      INTEGER PX(-1:MAXE), B(-1:MAXDEG)
      INTEGER V(0:MAXV)
      INTEGER E, I, J, R, U
C
C   Prepare to read irreducible polynomials on unit 1.
C
      OPEN (UNIT=GFUNIT, FILE='irrtabs.dat', STATUS='old')
C
C      *****  OPEN statements are system dependent
C
   10 READ (GFUNIT, 900, END=500) I
  900 FORMAT (20I3)
      IF (I .NE. Q) THEN
        DO 20 J = 1, NPOLS
   20     READ (GFUNIT, 900)
        GOTO 10
      ENDIF
C
      DO 1000 I = 1, DIMEN
C
C For each dimension, we need to calculate powers of an
C appropriate irreducible polynomial :  see Niederreiter
C page 65, just below equation (19).
C   Read the appropriate irreducible polynomial into PX,
C and its degree into E.  Set polynomial B = PX ** 0 = 1.
C M is the degree of B.  Subsequently B will hold higher
C powers of PX.
C   The polynomial PX is stored in file 'irrtabs.dat' in the
C format
C   n  a0  a1  a2  ... an
C where n is the degree of the polynomial and the ai are
C its coefficients.
C
        READ (GFUNIT, 900) E, (PX(J), J = 0,E)
        PX(DEG) = E
        B(DEG) = 0
        B(0) = 1
C
C Niederreiter (page 56, after equation (7), defines two
C variables Q and U.  We do not need Q explicitly, but we
C do need U.
C
        U = 0
C
        DO 90 J = 1, NFIGS
C
C  If U = 0, we need to set B to the next power of PX
C  and recalculate V.  This is done by subroutine CALCV.
C
          IF (U .EQ. 0) CALL CALCV (PX, B, V, MAXV)
C
C Now C is obtained from V.  Neiderreiter
C obtains A from V (page 65, near the bottom), and
C then gets C from A (page 56, equation (7)).
C However this can be done in one step.
C
          DO 50 R = 0, NFIGS-1
            C(I,J,R) = V(R+U)
   50     CONTINUE
C
C Increment U.  If U = E, then U = 0 and in Niederreiter's
C paper Q = Q + 1.  Here, however, Q is not used explicitly.
C
          U = U + 1
          IF (U .EQ. E) U = 0
  90    CONTINUE
C
 1000 CONTINUE
C
      CLOSE (GFUNIT)
      RETURN
C
  500 WRITE (*,*) ' CALCC :  Tables for q =', Q, ' not found'
      STOP
      END
C
C     ***** end of SUBROUTINE CALCC

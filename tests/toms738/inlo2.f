      SUBROUTINE INLO2 (DIM, SKIP)
C
C   This version :  12 February 1992
C
C   See the general comments on implementing Niederreiter's
C   low-discrepancy sequences.
C
C   This subroutine calculates the values of Niederreiter's
C   C(I,J,R) and performs other initialisation necessary
C   before calling GOLO2.
C
C INPUT :
C   DIMEN - The dimension of the sequence to be generated.
C        {DIMEN is called DIM in the argument of INLO2,
C        because DIMEN is subsequently passed via COMMON
C        and is called DIMEN there.}
C
C   SKIP  - The number of values to throw away at the beginning
C           of the sequence.
C
C OUTPUT :
C   To GOLO2, labelled common /COMM2/.
C
C USES :
C   Calls CALCC2 to calculate the values of CJ.
C   ***** A non-standard function is used to compute *****
C   ***** bitwise exclusive-ors.                     *****
C
C
C   ------------------------------------------------------------
C
C
C   This file defines the common block /COMM2/ and some
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
      INTEGER CJ(MAXDIM, 0:NBITS - 1), DIMEN, COUNT
      INTEGER NEXTQ(MAXDIM)
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
C   ------------------------------------------------------------
C
C
C
      INTEGER I, R, DIM, SKIP, GRAY
C
      DIMEN = DIM
C
C       This assignment just relabels the variable for
C       subsequent use.
C
      IF (DIMEN.LE.0 .OR. DIMEN.GT.MAXDIM) THEN
        WRITE (*,*) ' INLO2 :  Bad dimension'
        STOP
      ENDIF
C
      CALL CALCC2
C
C   Translate SKIP into Gray code
C
C   ***** The bitwise exclusive-or is not standard in Fortran
C   ***** This is the Vax version :
      GRAY = IEOR (SKIP, SKIP/2)
C   ***** THIS is the Unix version
C     GRAY = XOR (SKIP, SKIP/2)
C
C   Now set up NEXTQ appropriately for this value of the Gray code
C
      DO 5 I = 1, DIMEN
    5   NEXTQ(I) = 0
C
      R = 0
   10 IF (GRAY .NE. 0) THEN
        IF (MOD(GRAY,2) .NE. 0) THEN
          DO 20 I = 1, DIMEN
C           ***** This is the non-standard exclusive-or again
C           ***** Vax version :
            NEXTQ(I) = IEOR(NEXTQ(I), CJ(I,R))
C           ***** Unix version :
C           NEXTQ(I) = XOR(NEXTQ(I), CJ(I,R))
   20     CONTINUE
        ENDIF
        GRAY = GRAY / 2
        R = R + 1
        GOTO 10
      ENDIF
C
      COUNT = SKIP
      RETURN
      END
C
C     *****  end of SUBROUTINE INLO2

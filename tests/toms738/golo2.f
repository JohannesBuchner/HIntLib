      SUBROUTINE GOLO2 (QUASI)
C
C   This version :  21 February 1992
C
C        This routine is for base 2 only.  The driver, GENIN2,
C        calls it after proper set-up.
C
C See the general comments on implementing Niederreiter's
C low-discrepancy sequences.
C
C This subroutine generates a new quasi-random vector
C on each call.
C
C INPUT
C   From INLO2, labelled common /COMM2/, properly initialized.
C
C OUTPUT
C   To the caller, the next vector in the sequence in the
C   array QUASI.
C
C USES
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
C   ------------------------------------------------------------
C
C
C
      REAL RECIP
      PARAMETER (RECIP=2.0**(-NBITS))
C
C   The parameter RECIP is the multiplier which changes the
C   integers in NEXTQ into the required real values in QUASI.
C
      INTEGER I, R
      REAL    QUASI(*)
C
C Multiply the numerators in NEXTQ by RECIP to get the next
C quasi-random vector
C
      DO 5 I = 1, DIMEN
        QUASI(I) = NEXTQ(I) * RECIP
    5 CONTINUE
C
C Find the position of the right-hand zero in COUNT.  This
C is the bit that changes in the Gray-code representation as
C we go from COUNT to COUNT+1.
C
      R = 0
      I = COUNT
   10 IF (MOD(I,2).NE.0) THEN
        R = R + 1
        I = I/2
        GOTO 10
      ENDIF
C
C Check that we have not passed 2**NBITS calls on GOLO2
C
      IF (R .GE. NBITS) THEN
        WRITE (*,*) ' GOLO2 :  Too many calls'
        STOP
      ENDIF
C
C Compute the new numerators in vector NEXTQ
C
      DO 20 I = 1, DIMEN
C       ***** Bitwise exclusive-or is not standard in Fortran
C       ***** This is the Vax version :
       NEXTQ(I) = IEOR(NEXTQ(I), CJ(I,R))
C       ***** This is the Unix version
C      NEXTQ(I) = XOR(NEXTQ(I), CJ(I,R))
   20 CONTINUE
C
      COUNT = COUNT + 1
      RETURN
      END
C
C     *****   end of PROCEDURE GOLO2

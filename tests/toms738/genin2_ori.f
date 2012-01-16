      PROGRAM GENIN2
C
C        *****  This is the driver for the base-2 programs.
C        *****  More accurately, it is a driver skeleton.
C        *****  There is a default set of test integrals
C        *****  in TESTF.  The user can replace TESTF with
C        *****  another subroutine called TESTF containing
C        *****  integrals of interest to him.
C
C   This version :  2 Feb 1992
C
C   This program tests the accuracy of numerical integration
C   using the low-discrepancy binary sequences of
C   H. Niederreiter (1988) as implemented in INLO2, GOLO2, and
C   related programs.  Various possible test integrals are
C   provided by the function TESTF.  GENIN2 generates only
C   sequences with base 2.
C
C   Interactive input and output go through the Fortran units
C   denoted by * ;  at most installations, * corresponds to
C   the user's terminal.  Fortran unit 2 is used to save the output.
C
C   These programs assume that your computer's word length
C   is 31 bits, excluding sign.  If this is not the case,
C   modify the parameter NBITS throughout accordingly.
C
C
C
C    GENIN2 and its associated subroutines are listed below.
C
C       GENIN2
C          INLO2
C          GOLO2
C          CALCC2
C          CALCV
C          CHARAC
C          SETFLD
C          PLYMUL
C          TESTF
C
C    The suffix 2 indicates  routines for use only by
C    the set of programs tailored for base 2.  
C
C    The other routines are also used by the general-base programs.
C
      INTEGER MAXDIM, OUTUNT, POWER
      PARAMETER (MAXDIM=12, OUTUNT=2, POWER  = 12)
C
C   The parameter MAXDIM gives the maximum dimension that will
C   be used.  OUTUNT is the unit to save the output.
C   POWER is used in a possible warm-up formula.
C
      INTEGER I, NUM, DIMEN, SEQLEN, SKIP, STEP, PBASE
      REAL QUASI(MAXDIM), TESTF, EXACTF 
      DOUBLE PRECISION SUM
C
      WRITE (*,*) ' This is program GENIN2'
      WRITE (*,*) ' Output file name is OUTFIL.DAT '
      OPEN (unit = OUTUNT, file = 'OUTFIL.DAT', status = 'UNKNOWN')
C
C      *****   OPEN statements may well have to be modified
C      *****   to conform to local computer-center requirements.
C
C
    5 WRITE (*,*) ' Choose a test integral (1 to 4) or 0 to quit :'
      READ (*,*) NUM
      IF (NUM.LE.0) THEN
        WRITE (*,*) ' End of program GENIN2'
        CLOSE (OUTUNT)
        STOP
      ENDIF
      IF (NUM.GT.4) THEN
        WRITE (*,*) ' No such test integral'
        GOTO 5
      ENDIF
C
C       *****  Each test integral is parameterized by
C       *****  its dimension.
C
   10 WRITE (*,*) ' Enter dimension :'
      READ (*,*) DIMEN
      IF (DIMEN.GT.MAXDIM) THEN
        WRITE (*,*) ' Dimension may not exceed', MAXDIM
        GOTO 10
      ENDIF
C
C        *****  The sequence length is the number of
C        *****  quasi-random points used to estimate the
C        *****  integral, excluding warm-up.
C
C        *****  Some users may wish to rewrite
C        *****  the driver to test a [heuristic] "convergence"
C        *****  criterion, stopping the generation of points
C        *****  when it is passed or when a specified number of
C        *****  points have been generated  -- whichever occurs
C        *****  first.
C
   15 WRITE (*,*) ' Choose sequence length :'
C
C         *****  Except when comparing results with other
C         *****  bases, we suggest taking SEQLEN to be a power
C         *****  of 2.  Examples:
C
      WRITE (*,*) ' 2 ** 10 =  ',  2 ** 10
      WRITE (*,*) ' 2 ** 15 =  ',  2 ** 15
      WRITE (*,*) ' 2 ** 20 =  ',  2 ** 20
      WRITE (*,*) ' Enter SEQLEN  (possibly a power of 2 above) '
      READ (*,*) SEQLEN
      IF (SEQLEN.LE.0) THEN
        WRITE (*,*) ' Length must be strictly positive'
        GOTO 15
      ENDIF
C
   20 WRITE (*,*) ' Choose the number of values to skip :'
      WRITE (*,*) ' There is reason to believe that BASE * e, '
      WRITE (*,*) ' where e is defined for example in '
      WRITE (*,*) ' Bratley, Fox, and Niederreiter [1991], '
      WRITE (*,*) ' would be a good choice.  Our formula has '
      WRITE (*,*) ' has the form  SKIP = 2 ** POWER, where '
      WRITE (*,*) ' POWER is  chosen so that SKIP comes nowhere '
      WRITE (*,*) ' near the largest possible machine-representable'
      WRITE (*,*) ' integer.  It does not hurt to choose '
      WRITE (*,*) ' POWER larger than e, because skipping is '
      WRITE (*,*) ' done implicitly in O(1) time. '
      WRITE (*,*) ' The maximum value of e for a fixed dimension '
      WRITE (*,*) ' s grows like log s.  We allow some "fat" for '
      WRITE (*,*) ' the implicit constant. '
      WRITE (*,*) ' Numerically, 2 ** POWER = ', 2 ** POWER
      WRITE (*,*) ' Enter SKIP (not necessarily the value above)'
      READ (*,*) SKIP
      IF (SKIP.LT.0) THEN
        WRITE (*,*) ' Number must be nonnegative'
        GOTO 20
      ENDIF
C
C
      CALL INLO2 (DIMEN, SKIP)
      WRITE (*,*) ' GENIN2 :  Initialization complete'
C
C Write titles and the exact value of the integral
C
      WRITE (OUTUNT,27) NUM
   27 FORMAT (/,' Test integral ',I2)
      WRITE (OUTUNT,28) DIMEN, SEQLEN, SKIP
   28 FORMAT (/,' Dimension ',I6,',    Base 2 (GENIN2)',
     1 /,' Sequence ',I7,',    Skipped ',I4)
      WRITE (OUTUNT,29) EXACTF(NUM, DIMEN)
   29 FORMAT (/,' Correct value is ',G16.7)
      WRITE (OUTUNT,30)
   30 FORMAT(/,'      Iteration     Estimate of integral',/)
C
C Now estimate the integral
C
      SUM = 0
      PBASE = 2 ** 6
      WRITE (*,*) ' Odd-looking iteration numbers are powers of 2 '
      STEP = 500
      DO 100 I = 1, SEQLEN
        CALL GOLO2 (QUASI)
        SUM = SUM + TESTF(NUM, DIMEN, QUASI)
        IF (MOD(I,STEP).EQ.0) THEN
          WRITE (OUTUNT,99) I, SUM/I
        ENDIF
        IF (MOD(I,PBASE) .EQ. 0) THEN
          WRITE (OUTUNT,99) I, SUM/I
          PBASE = 2 * PBASE
C
C            There is reason to believe that the subsequence
C            of estimates along powers of the base [here 2]
C            converges faster than the original sequence or
C            the subsequence corresponding to STEP.
C
        ENDIF
C
   99     FORMAT (I12,G24.7)
          IF (I .EQ. 5000) STEP = 1000
          IF (I .EQ. 10000) STEP = 5000
  100 CONTINUE
C
      WRITE (*,*) ' GENIN2 :  iteration ends'
      GOTO 5
C
      END
C
C     ***** end of PROGRAM GENIN2

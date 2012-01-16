      PROGRAM GENIN
C
C        *****  This is the driver for the general base programs.
C        *****  More accurately, it is a driver skeleton.
C        *****  There is a default set of test integrals
C        *****  in TESTF.  The user can replace TESTF with
C        *****  another subroutine called TESTF containing
C        *****  integrals of interest to him.
C
C   This version :  14 Aug 1992
C
C   This program tests the accuracy of numerical integration
C   using the low-discrepancy binary sequences of
C   H. Niederreiter (1988) as implemented in INLO, GOLO, and
C   related programs.  Various possible test integrals are
C   provided by the function TESTF.
C
C   Interactive input and output go through the Fortran units
C   denoted by * ;  at most installations, * corresponds to
C   the user's terminal.  For a prime-power base, Fortran unit 1
C   is used to read field arithmetic tables ;  for any base, it
C   is used to read irreducible polynomials.
C   Fortran unit 2 is used to save the output.
C
C   Our VAX implementation does not measure elapsed time.
C   It can be modified, in a system-dependent way, to do so.
C
C   GENIN and its associated subroutines are listed below.
C   An asterisk denotes a subroutine also used by GFARIT
C   (to set up the field-arithmetic tables gftabs.dat) and
C   by GFPLYS (to set up the irreducible polynomials in
C   irrtabs.dat).  
C
C       GENIN
C          INLO
C          CALCC
C          CALCV          %
C          CHARAC *       %
C          SETFLD *       %
C          PLYMUL *       %
C          GOLO
C          TESTF          %
C
C      A percent sign above denotes a routine also used
C      in the set of programs tailored for base 2.
C
C
C     Both the general-base and base-2 programs assume that
C     your computer's word length is 31 bits, excluding sign.
C     If this is not the case, modify the parameter NBITS
C     throughout the PARAMETER statements in this set of
C     programs accordingly.
C
      INTEGER MAXDIM, OUTUNT, MAXBAS, READY
      PARAMETER (MAXDIM=12, OUTUNT=2, MAXBAS = 13)
C
C   The parameter MAXDIM gives the maximum dimension that will
C   be used.  OUTUNT is the unit to save the output.
C   MAXBAS gives the maximum asymptotically-optimal base
C   used up to MAXDIM.  The user enters an appropriate value
C   of READY depending on whether or not the required files
C   indicated above have been set up.
C
      INTEGER I, NUM, DIMEN, SEQLEN, SKIP, STEP, BASE, CHARAC
      INTEGER OPTBAS(2:MAXDIM), PBASE, POWER(2:MAXBAS)
      REAL QUASI(MAXDIM), TESTF, EXACTF
      DOUBLE PRECISION SUM
C
C     The DATA statement below gives the asymptotically-optimal
C     base for the respective dimension.
C
      DATA  (OPTBAS(I), I = 2,MAXDIM) / 2,3,3,5,7,7,9,9,11,11,13 /
C
C
C     The data statement below gives values used in a possible
C     warm-up calculation.
C
      DATA  (POWER(I), I = 2,MAXBAS) / 12,8,8,6,6,6,4,4,4,4,4,4 /
C
C     There are theoretical reasons to believe that BASE ** e,
C     where e is defined for example in Bratley, Fox, and
C     Niederreiter (1991), would be a good choice.  However,
C     we don't want to come anywhere near the largest possible
C     machine-representable integer; hence, the compromise
C     exponents above.  Note: subject to this conditon,
C     it can't hurt to take an exponent greater than e, because
C     warm-up skipping of initial values is done implicitly
C     in O(1) time.  The maximum value of e for a fixed dimension
C     s grows like log s.  We allow some "fat" for the implicit
C     constant in our choice of POWER.
C
C     WRITE (*,*) ' This is program GENIN'
C
C     WRITE(*,*) ' If you wish to use the base 2, the '
C     WRITE(*,*) ' alternative set of programs tailored for '
C     WRITE(*,*) ' base 2 will run much faster. '
C
C     WRITE (*,*) ' If the files gftabs.dat and irrtab.dat '
C     WRITE (*,*) ' have not already been set up, then '
C     WRITE (*,*) ' consult the guide to the programs for '
C     WRITE (*,*) ' how to do so, set up those files per '
C     WRITE (*,*) ' the guide, and [temporarily] quit this '
C     WRITE (*,*) ' set of programs by entering the integer 0 '
C     WRITE (*,*) ' below; otherwise, enter any other integer. '
C     WRITE (*,*) ' ENTER an appropriate integer. '
C     READ  (*,*)   READY
C     IF (READY .EQ. 0) THEN
C         WRITE (*,*) ' Set up gftabs.dat and irrtab.dat now. '
C         WRITE (*,*) ' Exit GENIN. '
C         STOP
C     ENDIF
C
C     WRITE (*,*)  'If the number of bit per word, excluding sign'
C     WRITE (*,*)  'is not 31, enter the integer 0 below '
C     WRITE (*,*)  'and fix the parameter NBITS everywhere; '
C     WRITE (*,*)  'otherwise, enter a positive integer below.'
C     WRITE (*,*)  'ENTER an appropriate integer. '
C     READ  (*,*)   READY
C     IF (READY .EQ. 0) THEN
C         WRITE(*,*) 'Fix NBITS.'
C         WRITE(*,*) 'Exit GENIN.'
C         STOP
C     ENDIF
C
C     WRITE (*,*) ' Output file name is OUTFIL.DAT'
C     OPEN (unit = OUTUNT, file = 'OUTFIL.DAT', status = 'UNKNOWN')
C
C     *****  OPEN statements are system-dependent.
C     *****  Therefore the statement above may have to be
C     *****  modified for use at your computer center.
C
C
C   5 WRITE (*,*) ' Choose a test integral (1 to 4) or 0 to quit :'
C     READ (*,*) NUM
C     IF (NUM.LE.0) THEN
C       WRITE (*,*) ' End of program GENIN'
C       CLOSE (OUTUNT)
C       STOP
C     ENDIF
C     IF (NUM.GT.4) THEN
C       WRITE (*,*) ' No such test integral'
C       GOTO 5
C     ENDIF
C
C       *****  Each test integral is parameterized by
C       *****  its dimension.
C
C  10 WRITE (*,*) ' Enter dimension :'
      READ (*,*) DIMEN
      IF (DIMEN.GT.MAXDIM) THEN
        WRITE (*,*) ' Dimension may not exceed', MAXDIM
        STOP
      ENDIF
C
C  12 WRITE (*,*) ' Choose a prime or prime-power base .'
C     WRITE (*,*) ' The asymptotically-optimal base '
C     WRITE (*,*) ' for this dimension is OPTBAS(DIMEN) = '
C     WRITE (*,*)   OPTBAS(DIMEN)
C     WRITE (*,*) ' This base may not be empirically optimal. '
C     WRITE (*,*) ' Enter BASE: '
C
      READ (*,*) BASE
      IF (CHARAC(BASE) .EQ. 0) THEN
        WRITE (*,*) ' Out of range or bad value :  try again'
        STOP
      ENDIF
C
C
C        *****  The sequence length is the number of
C        *****  quasi-random points used to estimate the
C        *****  integral, excluding warm-up.
C        *****  The number of initial quasi-random points
C        *****  deleted during warm-up is given by SKIP,
C        *****  chosen below.
C
C        *****  Some users may wish to rewrite
C        *****  the driver to test a [heuristic] "convergence"
C        *****  criterion, stopping the generation of points
C        *****  when it is passed or when a specified number of
C        *****  points have been generated  -- whichever occurs
C        *****  first.
C
C  15 WRITE (*,*) ' Choose sequence length :'
C     WRITE (*,*) ' A power of the base is recommended; e.g.,  '
C     WRITE (*,*) BASE ** POWER(BASE)
C     WRITE (*,*) BASE ** ((POWER(BASE) + 1))
C     WRITE (*,*) BASE ** ((POWER(BASE) + 2)) 
C     WRITE (*,*) BASE ** ((POWER(BASE) + 3))
C     WRITE (*,*) ' Enter SEQLEN '
      READ (*,*) SEQLEN
      IF (SEQLEN.LE.0) THEN
        WRITE (*,*) ' Length must be strictly positive'
        STOP
      ENDIF
C
C  20 WRITE (*,*) ' Choose the number of values to skip.'
C     WRITE (*,*) ' One possibility is given by the heuristic '
C     WRITE (*,*) ' formula SKIP = BASE ** POWER(BASE) '
C     WRITE (*,*) ' when BASE <= MAXBAS, = 10000 otherwise '
C     IF (BASE .LE. MAXBAS) THEN
C         SKIP = BASE ** POWER(BASE)
C       ELSE
C         SKIP = 10000
C     ENDIF
C     WRITE (*,*) ' Numerically, this equals ', SKIP
C     WRITE (8,*) ' Enter SKIP (not necessarily the value above) : '
C
      READ (*,*) SKIP
      IF (SKIP.LT.0) THEN
        WRITE (*,*) ' Number must be nonnegative'
        STOP
      ENDIF
C
C
C
      CALL INLO (DIMEN, BASE, SKIP)
C     WRITE (*,*) ' GENIN :  Initialization complete'
C
C Write title and the exact value of the integral
C
C     WRITE (OUTUNT,27) NUM
C  27 FORMAT (/,' Test integral ',I2)
C     WRITE (OUTUNT,28) DIMEN, BASE, SEQLEN, SKIP
C  28 FORMAT (/,' Dimension ',I6,',    Base ', I9,
C    1 /,' Sequence ',I8,',    Skipped ',I7)
C     WRITE (OUTUNT,29) EXACTF(NUM, DIMEN)
C  29 FORMAT (/,' Correct value is ',G16.7)
C     WRITE (OUTUNT,30)
C  30 FORMAT(/,'      Iteration     Estimate of integral',/)
C
C Now estimate the integral
C
C     WRITE (*,*)  ' The odd-looking iteration numbers '
C     WRITE (*,*)  ' in the output are powers of the base.  '
C
C     PBASE = BASE
C     SUM = 0
C     STEP = 500
      DO 100 I = 1, SEQLEN
        CALL GOLO (QUASI)
        WRITE(*,*) (QUASI(J),J=1,DIMEN)
C       SUM = SUM + TESTF(NUM, DIMEN, QUASI)
C       IF (MOD(I,STEP).EQ.0) THEN
C         WRITE (OUTUNT,99) I, SUM/I
C       ENDIF
C       IF (MOD(I,PBASE) .EQ. 0)  THEN
C             WRITE (OUTUNT,99) I, SUM/I
C             PBASE = PBASE * BASE
C
C         This finds the next power of the base.
C         There is reason to believe that convergenence
C         properties of the sequence of estimates is
C         better along the subsequence corrsponding to
C         powers of the base.
C
C       ENDIF
C  99     FORMAT (I12,G24.7)
C         IF (I .EQ. 5000) STEP = 1000
C         IF (I .EQ. 10000) STEP = 5000
  100 CONTINUE
C
C     WRITE (*,*) ' GENIN :  iteration ends'
C     GOTO 5
C
      END
C
C     ***** end of PROGRAM GENIN

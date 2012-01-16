      PROGRAM GFPLYS
C
C This version :  12 December 1991
C
C   The program calculates irreducible polynomials for various
C   finite fields, and writes them out to the file 'irrtabs.dat'.
C   Finite field arithmetic is carried out with the help of
C   precalculated addition and multiplication tables found on
C   the file 'gftabs.dat'.  The format of the irreducible polys on
C   the output file is
C       Q
C       d1   a1  a2 ... a(d1)
C       d2   b1  b2 ... b(d2)
C           etc
C   where Q is the order of the field, d1 is the degree of the
C   first irreducible poly, a1, a2, ..., a(d1) are its
C   coefficients, and so on.
C
C   The file 'gftabs.dat' is read on unit 1.  GFPLYS and the
C   associated subroutines assume that GFARIT has been run previously
C   to put the required data in this file.  GFPLYS writes its
C   output on file 'irrtabs.dat' using unit 2.
C
C   GFPLYS and its associated subroutines are listed below.
C   An asterisk indicates a subroutine also used by GFARIT
C   and by GENIN.  An ampersand indicates a subroutine also
C   used by GFARIT but not by GENIN.
C
C      GFPLYS
C         IRRED       (reads unit 1, writes unit 2)
C         CHARAC *
C         SETFLD *
C         ITOP   &
C         PTOI   &
C         PLYMUL *
C         FIND   
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
      INTEGER I
C
C   gftabs.dat should NOT be opened here, because it is opened in setfld.f
C
C     OPEN (UNIT=1, FILE='gftabs.dat', STATUS = 'old')
      OPEN (UNIT=2, FILE='irrtabs.dat', STATUS = 'unknown')
C
C      *****  OPEN statements are system dependent
C
      WRITE (*,*) ' This is GFPLYS'
      DO 10 I = 2, MAXQ
        CALL IRRED(I)
        WRITE (*,*) ' GFPLYS :  Case ', I, ' done'
   10 CONTINUE
C
      END
C
C     *****  end of PROGRAM GFPLYS
      SUBROUTINE IRRED (QIN)
C
C   This version :  12 December 1991
C
C
C   This routine reads the file gftabs.dat from unit 1
C   and writes the file irrtabs.dat onto unit 2.
C   GFPLYS opens units 1 and 2.
C   SETFLD opens unit 1, reads gftabs.dat, and then closes unit 1
C   on each call.
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
      INTEGER NPOLS, MAXS
      PARAMETER (NPOLS=25, MAXS=400)
      LOGICAL SEIVE(MAXS)
      INTEGER MONPOL(MAXS)
C
C   The parameter NPOLS says how many irreducible polys are to
C   be calculated for each field.
C   We find the irreducible polys using a seive.  Parameter
C   MAXS gives the size of this seive.  Array MONPOL holds monic
C   polys, array SEIVE says whether the poly is still OK.
C
      INTEGER QIN, I, J, K, N, PTOI, FIND, CHARAC
      INTEGER PI(-1:MAXDEG), PJ(-1:MAXDEG), PK(-1:MAXDEG)
C
      IF (QIN .LE. 1 .OR. QIN .GT. MAXQ) THEN
        WRITE (*,*) ' IRRED :  Bad value of Q'
        RETURN
      ENDIF
C
      P = CHARAC(QIN)
C
C   If no field of order QIN exists, there is nothing to do.
C
      IF (P .EQ. 0) RETURN
C
C   Otherwise, set up the field arithmetic tables.
C
      CALL SETFLD (QIN)
C
C   Set up seive containing only monic polys
C
      I = 0
      J = 1
      K = Q
      DO 50 N = 1, MAXS
        I = I + 1
        IF (I .EQ. J) THEN
          I = K
          J = 2 * K
          K = Q * K
        ENDIF
        MONPOL(N) = I
        SEIVE(N) = .TRUE.
   50 CONTINUE
C
C   Write out irreducible polys as they are found
C
      N = 0
      WRITE (2, 900) QIN
  900 FORMAT (20I3)
      DO 200 I = 1, MAXS
        IF (SEIVE(I)) THEN
          CALL ITOP (MONPOL(I), PI)
          K = PI(DEG)
          WRITE (2, 900) K, (PI(J), J = 0, K)
          N = N + 1
          IF (N .EQ. NPOLS) RETURN
C
          DO 100 J = I, MAXS
            CALL ITOP (MONPOL(J), PJ)
            CALL PLYMUL (PI, PJ, PK)
            K = FIND (PTOI (PK), MONPOL, J, MAXS)
            IF (K .NE. 0) SEIVE(K) = .FALSE.
  100     CONTINUE
        ENDIF
  200 CONTINUE
C
      WRITE (*,*) ' IRRED :  Seive too small'
      WRITE (*,*) ' Only', N, ' irreducible polys were found'
      RETURN
C
      END
C
C     *****   end of SUBROUTINE IRRED
      INTEGER FUNCTION FIND (N, TAB, I, MAXTAB)
C
C   This version :  12 December 1991
C
      INTEGER N, I, MAXTAB, TAB(MAXTAB), J
C
C   Look up N in ordered TAB(I) to TAB(MAXTAB)
C
      FIND = 0
      IF (N .GT. TAB(MAXTAB)) RETURN
      DO 10 J = I, MAXTAB
        IF (TAB(J) .EQ. N) THEN
          FIND = J
          RETURN
        ENDIF
   10 CONTINUE
      RETURN
      END
C
C     *****  end of INTEGER FUNCTION FIND

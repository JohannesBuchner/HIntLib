      PROGRAM GFARIT
C
C   This version :  12 December 1991
C
C   The program calculates addition and multiplication tables
C   for arithmetic in finite fields, and writes them out to
C   the file 'gftabs.dat'.  Tables are only calculated for fields
C   of prime-power order Q, the other cases being trivial.
C   For each value of Q, the file contains first Q, then the
C   addition table, and lastly the multiplication table.
C
C
C   GFARIT and its associated subroutines below must be
C   run to set up the file gftabs.dat [on unit 1].  
C   After gftabs.dat has been set up, run GFPLYS with its
C   associated subroutine to set up the file irrtabs.dat
C   [on unit 2]; this requires
C   reading gftabs.dat from unit 1.  The files gftabs.dat
C   and irrtabs.dat can be saved for future use.  Thus, each
C   user [or set of users with access to these files] needs to
C   run the respective sets of programs associated with GFARIT
C   and GFPLYS just once.  This must be done before running the
C   set of programs associated with GENIN.  The set of programs
C   tailored for base 2, using the driver GENIN2, requires
C   neither gftabs.dat nor irrtabs.dat, hence neither GFARIT 
C   nor GFPLYS.
C
C   Below we list [the main program] GFARIT and its associated
C   subroutines.  An asterisk indicates a subroutine also used
C   by GFPLYS and by GENIN.  An ampersand denotes a subroutine
C   also used by GFPLYS but not by GENIN.
C
C       GFARIT
C          GFTAB           (writes unit 1)
C          CHARAC *
C          SETFLD *
C          ITOP   &
C          PTOI   &
C          PLYADD 
C          PLYMUL *
C          PLYDIV
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
      INTEGER I
C
      OPEN (UNIT=1, FILE='gftabs.dat', STATUS = 'unknown')
C
C       ******  OPEN statements are system dependent
C
      WRITE (*,*) ' This is GFARIT'
      DO 10 I = 2, MAXQ
        CALL GFTAB(I)
        WRITE (*,*) ' GFARIT :  Case ', I, ' done'
   10 CONTINUE
C
      WRITE (*,*) ' GFARIT ends'
      STOP
      END
C
C     *****  end of PROGRAM GFARIT
      SUBROUTINE GFTAB (QIN)
C
C   This version :  12 December 1991
C
C
C    This routine writes gftabs.dat onto unit 1.
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
C   The common /FIELD/ will be set up to work modulo P, the
C   characteristic of field QIN.  We construct the tables for
C   the field of order QIN in GFADD and GFMUL.
C
      INTEGER QIN, I, J, PTOI, CHARAC
      INTEGER PI(-1:MAXDEG), PJ(-1:MAXDEG), PK(-1:MAXDEG)
      INTEGER GFADD(0:MAXQ-1, 0:MAXQ-1), GFMUL(0:MAXQ-1, 0:MAXQ-1)
C
C   IRRPLY holds irreducible polynomials for constructing
C   prime-power fields.  IRRPLY(-2,I) says which field this
C   row is used for, and then the rest of the row is a
C   polynomial (with the degree in IRRPLY(-1,I) as usual).
C   The chosen irreducible poly is copied into MODPLY for use.
C
      INTEGER IRRPLY(8, -2:MAXDEG), MODPLY(-1:MAXDEG)
      SAVE IRRPLY
C
      DATA (IRRPLY(1,J), J=-2,2) /4, 2, 1, 1, 1/
      DATA (IRRPLY(2,J), J=-2,3) /8, 3, 1, 1, 0, 1/
      DATA (IRRPLY(3,J), J=-2,2) /9, 2, 1, 0, 1/
      DATA (IRRPLY(4,J), J=-2,4) /16, 4, 1, 1, 0, 0, 1/
      DATA (IRRPLY(5,J), J=-2,2) /25, 2, 2, 0, 1/
      DATA (IRRPLY(6,J), J=-2,3) /27, 3, 1, 2, 0, 1/
      DATA (IRRPLY(7,J), J=-2,5) /32, 5, 1, 0, 1, 0, 0, 1/
      DATA (IRRPLY(8,J), J=-2,2) /49, 2, 1, 0, 1/
C
      IF (QIN .LE. 1 .OR. QIN .GT. MAXQ) THEN
        WRITE (*,*) ' ARITH :  Bad value of Q'
        RETURN
      ENDIF
C
      P = CHARAC(QIN)
C
C   If QIN is not a prime-power, we are not interested.
C
      IF (P .EQ. 0 .OR. P .EQ. QIN) RETURN
C
C   Otherwise, we set up the elements of the common /FIELD/
C   ready to do arithmetic mod P, the characteristic of QIN.
C
      CALL SETFLD (P)
C
C   Next find a suitable irreducible polynomial and copy it
C   to array MODPLY.
C
      I = 1
   20 IF (IRRPLY(I,-2) .NE. QIN) THEN
        I = I + 1
        GOTO 20
      ENDIF
      DO 30 J = -1, IRRPLY(I, DEG)
   30   MODPLY(J) = IRRPLY(I, J)
      DO 40 J = IRRPLY(I,DEG)+1, MAXDEG
   40   MODPLY(J) = 0
C
C   Deal with the trivial cases ...
C
      DO 60 I = 0, QIN-1
        GFADD(I,0) = I
        GFADD(0,I) = I
        GFMUL(I,0) = 0
   60   GFMUL(0,I) = 0
      DO 70 I = 1, QIN-1
        GFMUL(I,1) = I
   70   GFMUL(1,I) = I
C
C   ... and now deal with the rest.  Each integer from 1 to QIN-1
C   is treated as a polynomial with coefficients handled mod P.
C   Multiplication of polynomials is mod MODPLY.
C
      DO 80 I = 1, QIN-1
        CALL ITOP (I, PI)
        DO 80 J = 1, I
          CALL ITOP (J, PJ)
          CALL PLYADD (PI, PJ, PK)
          GFADD(I,J) = PTOI (PK)
          GFADD(J,I) = GFADD(I,J)
          IF (I .GT. 1 .AND. J .GT. 1) THEN
            CALL PLYMUL (PI, PJ, PK)
            CALL PLYDIV (PK, MODPLY, PJ, PK)
            GFMUL(I,J) = PTOI (PK)
            GFMUL(J,I) = GFMUL(I,J)
          ENDIF
   80 CONTINUE
C
C   Write the newly-calculated tables out to unit 1
C
      WRITE (1, 900) QIN
C
C   This is the table gftabs.dat.
C   SETFLD opens unit 1, reads gftabs.dat, and then closes unit 1.
C
      DO 100 I = 0, QIN-1
  100   WRITE (1, 900) (GFADD(I,J), J = 0, QIN-1)
      DO 110 I = 0, QIN-1
  110   WRITE (1, 900) (GFMUL(I,J), J = 0, QIN-1)
  900 FORMAT (20I3)
C
      RETURN
      END
C
C     *****  end of SUBROUTINE GFTAB

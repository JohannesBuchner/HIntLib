      SUBROUTINE CALCV (PX, B, V, MAXV)
C
C   This version :  12 February 1991
C
C   See the general comments on implementing Niederreiter's
C   low-discrepancy sequences.
C
C   This program calculates the values of the constants V(J,R) as
C   described in BFN section 3.3.  It is called from either CALCC or
C   CALCC2.  The values transmitted through common /FIELD/ determine
C   which field we are working in.
C
C INPUT :
C   PX is the appropriate irreducible polynomial for the dimension
C     currently being considered.  The degree of PX will be called E.
C   B is the polynomial defined in section 2.3 of BFN.  On entry to
C     the subroutine, the degree of B implicitly defines the parameter
C     J of section 3.3, by degree(B) = E*(J-1).
C   MAXV gives the dimension of the array V.
C   On entry, we suppose that the common block /FIELD/ has been set
C     up correctly (using SETFLD).
C
C OUTPUT :
C   On output B has been multiplied by PX, so its degree is now E*J.
C   V contains the values required.
C
C USES :
C   The subroutine PLYMUL is used to multiply polynomials.
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
      INTEGER MAXV, E, I, J, KJ, M, BIGM, R, TERM
      INTEGER PX(-1:MAXDEG), B(-1:MAXDEG), V(0:MAXV)
      INTEGER H(-1:MAXDEG)
C
      INTEGER ARBIT, NONZER
      ARBIT() = 1
C
C   We use ARBIT() to indicate where the user can place
C   an arbitrary element of the field of order Q, while NONZER
C   shows where he must put an arbitrary non-zero element
C   of the same field.  For the code,
C   this means 0 <= ARBIT < Q and 0 < NONZER < Q.  Within these
C   limits, the user can do what he likes.  ARBIT is declared as
C   a function as a reminder that a different arbitrary value may
C   be returned each time ARBIT is referenced.
C
C    BIGM is the M used in section 3.3.
C    It differs from the [little] m used in section 2.3,
C    denoted here by M.
C
      NONZER = 1
C
      E = PX(DEG)
C
C   The poly H is PX**(J-1), which is the value of B on arrival.
C   In section 3.3, the values of Hi are defined with a minus sign :
C   don't forget this if you use them later !
C
      DO 10 I = -1, B(DEG)
   10   H(I) = B(I)
      BIGM = H(DEG)
C
C   Now multiply B by PX so B becomes PX**J.
C   In section 2.3, the values of Bi are defined with a minus sign :
C   don't forget this if you use them later !
C
      CALL PLYMUL (PX, B, B)
      M = B(DEG)
C
C   We don't use J explicitly anywhere, but here it is just in case.
C
      J = M / E
C
C   Now choose a value of Kj as defined in section 3.3.
C   We must have 0 <= Kj < E*J = M.
C   The limit condition on Kj does not seem very relevant
C   in this program.
C
      KJ = BIGM
C
C   Now choose values of V in accordance with the conditions in
C   section 3.3
C
      DO 20 R = 0, KJ-1
   20   V(R) = 0
      V(KJ) = 1
C
      IF (KJ .LT. BIGM) THEN
C
        TERM = SUB (0, H(KJ))
C
        DO 30 R = KJ+1, BIGM-1
          V(R) = ARBIT()
C
C         Check the condition of section 3.3,
C         remembering that the H's have the opposite sign.
C
          TERM = SUB (TERM, MUL (H(R), V(R)))
   30   CONTINUE
C
C         Now V(BIGM) is anything but TERM
C
          V(BIGM) = ADD (NONZER, TERM)
C
        DO 40 R = BIGM+1, M-1
   40     V(R) = ARBIT()
C
      ELSE
C       This is the case KJ .GE. BIGM
C
        DO 50 R = KJ+1, M-1
   50     V(R) = ARBIT()
C
      ENDIF
C
C   Calculate the remaining V's using the recursion of section 2.3,
C   remembering that the B's have the opposite sign.
C
      DO 70 R = 0, MAXV-M
        TERM = 0
        DO 60 I = 0, M-1
          TERM = SUB (TERM, MUL (B(I), V(R+I)))
   60   CONTINUE
        V(R+M) = TERM
   70 CONTINUE
C
      RETURN
      END
C
C     ***** end of SUBROUTINE CALCV

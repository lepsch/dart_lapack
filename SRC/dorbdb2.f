      SUBROUTINE DORBDB2( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LWORK, M, P, Q, LDX11, LDX21;
*     ..
*     .. Array Arguments ..
      double             PHI(*), THETA(*);
      double             TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), X11(LDX11,*), X21(LDX21,*);
*     ..
*
*  ====================================================================
*
*     .. Parameters ..
      double             NEGONE, ONE;
      PARAMETER          ( NEGONE = -1.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      double             C, S;
      int                CHILDINFO, I, ILARF, IORBDB5, LLARF, LORBDB5, LWORKMIN, LWORKOPT;
      bool               LQUERY;
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFGP, DORBDB5, DROT, DSCAL, XERBLA
*     ..
*     .. External Functions ..
      double             DNRM2;
      EXTERNAL           DNRM2
*     ..
*     .. Intrinsic Function ..
      // INTRINSIC ATAN2, COS, MAX, SIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test input arguments
*
      INFO = 0
      LQUERY = LWORK .EQ. -1
*
      IF( M .LT. 0 ) THEN
         INFO = -1
      ELSE IF( P .LT. 0 .OR. P .GT. M-P ) THEN
         INFO = -2
      ELSE IF( Q .LT. 0 .OR. Q .LT. P .OR. M-Q .LT. P ) THEN
         INFO = -3
      ELSE IF( LDX11 .LT. MAX( 1, P ) ) THEN
         INFO = -5
      ELSE IF( LDX21 .LT. MAX( 1, M-P ) ) THEN
         INFO = -7
      END IF
*
*     Compute workspace
*
      IF( INFO .EQ. 0 ) THEN
         ILARF = 2
         LLARF = MAX( P-1, M-P, Q-1 )
         IORBDB5 = 2
         LORBDB5 = Q-1
         LWORKOPT = MAX( ILARF+LLARF-1, IORBDB5+LORBDB5-1 )
         LWORKMIN = LWORKOPT
         WORK(1) = LWORKOPT
         IF( LWORK .LT. LWORKMIN .AND. .NOT.LQUERY ) THEN
           INFO = -14
         END IF
      END IF
      IF( INFO .NE. 0 ) THEN
         CALL XERBLA( 'DORBDB2', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Reduce rows 1, ..., P of X11 and X21
*
      DO I = 1, P
*
         IF( I .GT. 1 ) THEN
            CALL DROT( Q-I+1, X11(I,I), LDX11, X21(I-1,I), LDX21, C, S )
         END IF
         CALL DLARFGP( Q-I+1, X11(I,I), X11(I,I+1), LDX11, TAUQ1(I) )
         C = X11(I,I)
         X11(I,I) = ONE
         CALL DLARF( 'R', P-I, Q-I+1, X11(I,I), LDX11, TAUQ1(I), X11(I+1,I), LDX11, WORK(ILARF) )          CALL DLARF( 'R', M-P-I+1, Q-I+1, X11(I,I), LDX11, TAUQ1(I), X21(I,I), LDX21, WORK(ILARF) )          S = SQRT( DNRM2( P-I, X11(I+1,I), 1 )**2 + DNRM2( M-P-I+1, X21(I,I), 1 )**2 )
         THETA(I) = ATAN2( S, C )
*
         CALL DORBDB5( P-I, M-P-I+1, Q-I, X11(I+1,I), 1, X21(I,I), 1, X11(I+1,I+1), LDX11, X21(I,I+1), LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO )
         CALL DSCAL( P-I, NEGONE, X11(I+1,I), 1 )
         CALL DLARFGP( M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) )
         IF( I .LT. P ) THEN
            CALL DLARFGP( P-I, X11(I+1,I), X11(I+2,I), 1, TAUP1(I) )
            PHI(I) = ATAN2( X11(I+1,I), X21(I,I) )
            C = COS( PHI(I) )
            S = SIN( PHI(I) )
            X11(I+1,I) = ONE
            CALL DLARF( 'L', P-I, Q-I, X11(I+1,I), 1, TAUP1(I), X11(I+1,I+1), LDX11, WORK(ILARF) )
         END IF
         X21(I,I) = ONE
         CALL DLARF( 'L', M-P-I+1, Q-I, X21(I,I), 1, TAUP2(I), X21(I,I+1), LDX21, WORK(ILARF) )
*
      END DO
*
*     Reduce the bottom-right portion of X21 to the identity matrix
*
      DO I = P + 1, Q
         CALL DLARFGP( M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) )
         X21(I,I) = ONE
         CALL DLARF( 'L', M-P-I+1, Q-I, X21(I,I), 1, TAUP2(I), X21(I,I+1), LDX21, WORK(ILARF) )
      END DO
*
      RETURN
*
*     End of DORBDB2
*
      END

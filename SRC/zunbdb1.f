      SUBROUTINE ZUNBDB1( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LWORK, M, P, Q, LDX11, LDX21;
      // ..
      // .. Array Arguments ..
      double             PHI(*), THETA(*);
      COMPLEX*16         TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), X11(LDX11,*), X21(LDX21,*)
      // ..

*  ====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = (1.0D0,0.0D0) )
      // ..
      // .. Local Scalars ..
      double             C, S;
      int                CHILDINFO, I, ILARF, IORBDB5, LLARF, LORBDB5, LWORKMIN, LWORKOPT;
      bool               LQUERY;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLARF, ZLARFGP, ZUNBDB5, ZDROT, XERBLA
      // EXTERNAL ZLACGV
      // ..
      // .. External Functions ..
      double             DZNRM2;
      // EXTERNAL DZNRM2
      // ..
      // .. Intrinsic Function ..
      // INTRINSIC ATAN2, COS, MAX, SIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test input arguments

      INFO = 0
      LQUERY = LWORK .EQ. -1

      IF( M .LT. 0 ) THEN
         INFO = -1
      ELSE IF( P .LT. Q .OR. M-P .LT. Q ) THEN
         INFO = -2
      ELSE IF( Q .LT. 0 .OR. M-Q .LT. Q ) THEN
         INFO = -3
      ELSE IF( LDX11 .LT. MAX( 1, P ) ) THEN
         INFO = -5
      ELSE IF( LDX21 .LT. MAX( 1, M-P ) ) THEN
         INFO = -7
      END IF

      // Compute workspace

      IF( INFO .EQ. 0 ) THEN
         ILARF = 2
         LLARF = MAX( P-1, M-P-1, Q-1 )
         IORBDB5 = 2
         LORBDB5 = Q-2
         LWORKOPT = MAX( ILARF+LLARF-1, IORBDB5+LORBDB5-1 )
         LWORKMIN = LWORKOPT
         WORK(1) = LWORKOPT
         IF( LWORK .LT. LWORKMIN .AND. .NOT.LQUERY ) THEN
           INFO = -14
         END IF
      END IF
      IF( INFO .NE. 0 ) THEN
         CALL XERBLA( 'ZUNBDB1', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Reduce columns 1, ..., Q of X11 and X21

      DO I = 1, Q

         CALL ZLARFGP( P-I+1, X11(I,I), X11(I+1,I), 1, TAUP1(I) )
         CALL ZLARFGP( M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) )
         THETA(I) = ATAN2( DBLE( X21(I,I) ), DBLE( X11(I,I) ) )
         C = COS( THETA(I) )
         S = SIN( THETA(I) )
         X11(I,I) = ONE
         X21(I,I) = ONE
         CALL ZLARF( 'L', P-I+1, Q-I, X11(I,I), 1, DCONJG(TAUP1(I)), X11(I,I+1), LDX11, WORK(ILARF) )          CALL ZLARF( 'L', M-P-I+1, Q-I, X21(I,I), 1, DCONJG(TAUP2(I)), X21(I,I+1), LDX21, WORK(ILARF) )

         IF( I .LT. Q ) THEN
            CALL ZDROT( Q-I, X11(I,I+1), LDX11, X21(I,I+1), LDX21, C, S )
            CALL ZLACGV( Q-I, X21(I,I+1), LDX21 )
            CALL ZLARFGP( Q-I, X21(I,I+1), X21(I,I+2), LDX21, TAUQ1(I) )
            S = DBLE( X21(I,I+1) )
            X21(I,I+1) = ONE
            CALL ZLARF( 'R', P-I, Q-I, X21(I,I+1), LDX21, TAUQ1(I), X11(I+1,I+1), LDX11, WORK(ILARF) )             CALL ZLARF( 'R', M-P-I, Q-I, X21(I,I+1), LDX21, TAUQ1(I), X21(I+1,I+1), LDX21, WORK(ILARF) )
            CALL ZLACGV( Q-I, X21(I,I+1), LDX21 )
            C = SQRT( DZNRM2( P-I, X11(I+1,I+1), 1 )**2 + DZNRM2( M-P-I, X21(I+1,I+1), 1 )**2 )
            PHI(I) = ATAN2( S, C )
            CALL ZUNBDB5( P-I, M-P-I, Q-I-1, X11(I+1,I+1), 1, X21(I+1,I+1), 1, X11(I+1,I+2), LDX11, X21(I+1,I+2), LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO )
         END IF

      END DO

      RETURN

      // End of ZUNBDB1

      END

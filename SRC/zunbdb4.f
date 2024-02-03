      SUBROUTINE ZUNBDB4( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, TAUP1, TAUP2, TAUQ1, PHANTOM, WORK, LWORK, INFO )
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
      COMPLEX*16         PHANTOM(*), TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), X11(LDX11,*), X21(LDX21,*)
*     ..
*
*  ====================================================================
*
*     .. Parameters ..
      COMPLEX*16         NEGONE, ONE, ZERO
      PARAMETER          ( NEGONE = (-1.0D0,0.0D0), ONE = (1.0D0,0.0D0), ZERO = (0.0D0,0.0D0) )
*     ..
*     .. Local Scalars ..
      double             C, S;
      int                CHILDINFO, I, ILARF, IORBDB5, J, LLARF, LORBDB5, LWORKMIN, LWORKOPT;
      bool               LQUERY;
*     ..
*     .. External Subroutines ..
      // EXTERNAL ZLARF, ZLARFGP, ZUNBDB5, ZDROT, ZSCAL, ZLACGV, XERBLA
*     ..
*     .. External Functions ..
      double             DZNRM2;
      // EXTERNAL DZNRM2
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
      ELSE IF( P .LT. M-Q .OR. M-P .LT. M-Q ) THEN
         INFO = -2
      ELSE IF( Q .LT. M-Q .OR. Q .GT. M ) THEN
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
         LLARF = MAX( Q-1, P-1, M-P-1 )
         IORBDB5 = 2
         LORBDB5 = Q
         LWORKOPT = ILARF + LLARF - 1
         LWORKOPT = MAX( LWORKOPT, IORBDB5 + LORBDB5 - 1 )
         LWORKMIN = LWORKOPT
         WORK(1) = LWORKOPT
         IF( LWORK .LT. LWORKMIN .AND. .NOT.LQUERY ) THEN
           INFO = -14
         END IF
      END IF
      IF( INFO .NE. 0 ) THEN
         CALL XERBLA( 'ZUNBDB4', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Reduce columns 1, ..., M-Q of X11 and X21
*
      DO I = 1, M-Q
*
         IF( I .EQ. 1 ) THEN
            DO J = 1, M
               PHANTOM(J) = ZERO
            END DO
            CALL ZUNBDB5( P, M-P, Q, PHANTOM(1), 1, PHANTOM(P+1), 1, X11, LDX11, X21, LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO )
            CALL ZSCAL( P, NEGONE, PHANTOM(1), 1 )
            CALL ZLARFGP( P, PHANTOM(1), PHANTOM(2), 1, TAUP1(1) )
            CALL ZLARFGP( M-P, PHANTOM(P+1), PHANTOM(P+2), 1, TAUP2(1) )
            THETA(I) = ATAN2( DBLE( PHANTOM(1) ), DBLE( PHANTOM(P+1) ) )
            C = COS( THETA(I) )
            S = SIN( THETA(I) )
            PHANTOM(1) = ONE
            PHANTOM(P+1) = ONE
            CALL ZLARF( 'L', P, Q, PHANTOM(1), 1, DCONJG(TAUP1(1)), X11, LDX11, WORK(ILARF) )             CALL ZLARF( 'L', M-P, Q, PHANTOM(P+1), 1, DCONJG(TAUP2(1)), X21, LDX21, WORK(ILARF) )
         ELSE
            CALL ZUNBDB5( P-I+1, M-P-I+1, Q-I+1, X11(I,I-1), 1, X21(I,I-1), 1, X11(I,I), LDX11, X21(I,I), LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO )
            CALL ZSCAL( P-I+1, NEGONE, X11(I,I-1), 1 )
            CALL ZLARFGP( P-I+1, X11(I,I-1), X11(I+1,I-1), 1, TAUP1(I) )
            CALL ZLARFGP( M-P-I+1, X21(I,I-1), X21(I+1,I-1), 1, TAUP2(I) )
            THETA(I) = ATAN2( DBLE( X11(I,I-1) ), DBLE( X21(I,I-1) ) )
            C = COS( THETA(I) )
            S = SIN( THETA(I) )
            X11(I,I-1) = ONE
            X21(I,I-1) = ONE
            CALL ZLARF( 'L', P-I+1, Q-I+1, X11(I,I-1), 1, DCONJG(TAUP1(I)), X11(I,I), LDX11, WORK(ILARF) )             CALL ZLARF( 'L', M-P-I+1, Q-I+1, X21(I,I-1), 1, DCONJG(TAUP2(I)), X21(I,I), LDX21, WORK(ILARF) )
         END IF
*
         CALL ZDROT( Q-I+1, X11(I,I), LDX11, X21(I,I), LDX21, S, -C )
         CALL ZLACGV( Q-I+1, X21(I,I), LDX21 )
         CALL ZLARFGP( Q-I+1, X21(I,I), X21(I,I+1), LDX21, TAUQ1(I) )
         C = DBLE( X21(I,I) )
         X21(I,I) = ONE
         CALL ZLARF( 'R', P-I, Q-I+1, X21(I,I), LDX21, TAUQ1(I), X11(I+1,I), LDX11, WORK(ILARF) )          CALL ZLARF( 'R', M-P-I, Q-I+1, X21(I,I), LDX21, TAUQ1(I), X21(I+1,I), LDX21, WORK(ILARF) )
         CALL ZLACGV( Q-I+1, X21(I,I), LDX21 )
         IF( I .LT. M-Q ) THEN
            S = SQRT( DZNRM2( P-I, X11(I+1,I), 1 )**2 + DZNRM2( M-P-I, X21(I+1,I), 1 )**2 )
            PHI(I) = ATAN2( S, C )
         END IF
*
      END DO
*
*     Reduce the bottom-right portion of X11 to [ I 0 ]
*
      DO I = M - Q + 1, P
         CALL ZLACGV( Q-I+1, X11(I,I), LDX11 )
         CALL ZLARFGP( Q-I+1, X11(I,I), X11(I,I+1), LDX11, TAUQ1(I) )
         X11(I,I) = ONE
         CALL ZLARF( 'R', P-I, Q-I+1, X11(I,I), LDX11, TAUQ1(I), X11(I+1,I), LDX11, WORK(ILARF) )          CALL ZLARF( 'R', Q-P, Q-I+1, X11(I,I), LDX11, TAUQ1(I), X21(M-Q+1,I), LDX21, WORK(ILARF) )
         CALL ZLACGV( Q-I+1, X11(I,I), LDX11 )
      END DO
*
*     Reduce the bottom-right portion of X21 to [ 0 I ]
*
      DO I = P + 1, Q
         CALL ZLACGV( Q-I+1, X21(M-Q+I-P,I), LDX21 )
         CALL ZLARFGP( Q-I+1, X21(M-Q+I-P,I), X21(M-Q+I-P,I+1), LDX21, TAUQ1(I) )
         X21(M-Q+I-P,I) = ONE
         CALL ZLARF( 'R', Q-I, Q-I+1, X21(M-Q+I-P,I), LDX21, TAUQ1(I), X21(M-Q+I-P+1,I), LDX21, WORK(ILARF) )
         CALL ZLACGV( Q-I+1, X21(M-Q+I-P,I), LDX21 )
      END DO
*
      RETURN
*
*     End of ZUNBDB4
*
      END

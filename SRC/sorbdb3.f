      SUBROUTINE SORBDB3( M, P, Q, X11, LDX11, X21, LDX21, THETA, PHI, TAUP1, TAUP2, TAUQ1, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LWORK, M, P, Q, LDX11, LDX21;
      // ..
      // .. Array Arguments ..
      REAL               PHI(*), THETA(*)
      REAL               TAUP1(*), TAUP2(*), TAUQ1(*), WORK(*), X11(LDX11,*), X21(LDX21,*)
      // ..

*  ====================================================================

      // .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
      // ..
      // .. Local Scalars ..
      REAL               C, S
      int                CHILDINFO, I, ILARF, IORBDB5, LLARF, LORBDB5, LWORKMIN, LWORKOPT;
      bool               LQUERY;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARF, SLARFGP, SORBDB5, SROT, XERBLA
      // ..
      // .. External Functions ..
      REAL               SNRM2
      // EXTERNAL SNRM2
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
      ELSE IF( 2*P .LT. M .OR. P .GT. M ) THEN
         INFO = -2
      ELSE IF( Q .LT. M-P .OR. M-Q .LT. M-P ) THEN
         INFO = -3
      ELSE IF( LDX11 .LT. MAX( 1, P ) ) THEN
         INFO = -5
      ELSE IF( LDX21 .LT. MAX( 1, M-P ) ) THEN
         INFO = -7
      END IF

      // Compute workspace

      IF( INFO .EQ. 0 ) THEN
         ILARF = 2
         LLARF = MAX( P, M-P-1, Q-1 )
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
         CALL XERBLA( 'SORBDB3', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Reduce rows 1, ..., M-P of X11 and X21

      DO I = 1, M-P

         IF( I .GT. 1 ) THEN
            CALL SROT( Q-I+1, X11(I-1,I), LDX11, X21(I,I), LDX11, C, S )
         END IF

         CALL SLARFGP( Q-I+1, X21(I,I), X21(I,I+1), LDX21, TAUQ1(I) )
         S = X21(I,I)
         X21(I,I) = ONE
         CALL SLARF( 'R', P-I+1, Q-I+1, X21(I,I), LDX21, TAUQ1(I), X11(I,I), LDX11, WORK(ILARF) )          CALL SLARF( 'R', M-P-I, Q-I+1, X21(I,I), LDX21, TAUQ1(I), X21(I+1,I), LDX21, WORK(ILARF) )          C = SQRT( SNRM2( P-I+1, X11(I,I), 1 )**2 + SNRM2( M-P-I, X21(I+1,I), 1 )**2 )
         THETA(I) = ATAN2( S, C )

         CALL SORBDB5( P-I+1, M-P-I, Q-I, X11(I,I), 1, X21(I+1,I), 1, X11(I,I+1), LDX11, X21(I+1,I+1), LDX21, WORK(IORBDB5), LORBDB5, CHILDINFO )
         CALL SLARFGP( P-I+1, X11(I,I), X11(I+1,I), 1, TAUP1(I) )
         IF( I .LT. M-P ) THEN
            CALL SLARFGP( M-P-I, X21(I+1,I), X21(I+2,I), 1, TAUP2(I) )
            PHI(I) = ATAN2( X21(I+1,I), X11(I,I) )
            C = COS( PHI(I) )
            S = SIN( PHI(I) )
            X21(I+1,I) = ONE
            CALL SLARF( 'L', M-P-I, Q-I, X21(I+1,I), 1, TAUP2(I), X21(I+1,I+1), LDX21, WORK(ILARF) )
         END IF
         X11(I,I) = ONE
         CALL SLARF( 'L', P-I+1, Q-I, X11(I,I), 1, TAUP1(I), X11(I,I+1), LDX11, WORK(ILARF) )

      END DO

      // Reduce the bottom-right portion of X11 to the identity matrix

      DO I = M-P + 1, Q
         CALL SLARFGP( P-I+1, X11(I,I), X11(I+1,I), 1, TAUP1(I) )
         X11(I,I) = ONE
         CALL SLARF( 'L', P-I+1, Q-I, X11(I,I), 1, TAUP1(I), X11(I,I+1), LDX11, WORK(ILARF) )
      END DO

      RETURN

      // End of SORBDB3

      END

      SUBROUTINE SORBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, THETA, PHI, TAUP1, TAUP2, TAUQ1, TAUQ2, WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      String             SIGNS, TRANS;
      int                INFO, LDX11, LDX12, LDX21, LDX22, LWORK, M, P, Q;
      // ..
      // .. Array Arguments ..
      REAL               PHI( * ), THETA( * )
      REAL               TAUP1( * ), TAUP2( * ), TAUQ1( * ), TAUQ2( * ), WORK( * ), X11( LDX11, * ), X12( LDX12, * ), X21( LDX21, * ), X22( LDX22, * )
      // ..
*
*  ====================================================================
*
      // .. Parameters ..
      REAL               REALONE
      PARAMETER          ( REALONE = 1.0E0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
      // ..
      // .. Local Scalars ..
      bool               COLMAJOR, LQUERY;
      int                I, LWORKMIN, LWORKOPT;
      REAL               Z1, Z2, Z3, Z4
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SLARF, SLARFGP, SSCAL, XERBLA
      // ..
      // .. External Functions ..
      REAL               SNRM2
      bool               LSAME;
      // EXTERNAL SNRM2, LSAME
      // ..
      // .. Intrinsic Functions
      // INTRINSIC ATAN2, COS, MAX, SIN
      // ..
      // .. Executable Statements ..
*
      // Test input arguments
*
      INFO = 0
      COLMAJOR = .NOT. LSAME( TRANS, 'T' )
      IF( .NOT. LSAME( SIGNS, 'O' ) ) THEN
         Z1 = REALONE
         Z2 = REALONE
         Z3 = REALONE
         Z4 = REALONE
      ELSE
         Z1 = REALONE
         Z2 = -REALONE
         Z3 = REALONE
         Z4 = -REALONE
      END IF
      LQUERY = LWORK .EQ. -1
*
      IF( M .LT. 0 ) THEN
         INFO = -3
      ELSE IF( P .LT. 0 .OR. P .GT. M ) THEN
         INFO = -4
      ELSE IF( Q .LT. 0 .OR. Q .GT. P .OR. Q .GT. M-P .OR. Q .GT. M-Q ) THEN
         INFO = -5
      ELSE IF( COLMAJOR .AND. LDX11 .LT. MAX( 1, P ) ) THEN
         INFO = -7
      ELSE IF( .NOT.COLMAJOR .AND. LDX11 .LT. MAX( 1, Q ) ) THEN
         INFO = -7
      ELSE IF( COLMAJOR .AND. LDX12 .LT. MAX( 1, P ) ) THEN
         INFO = -9
      ELSE IF( .NOT.COLMAJOR .AND. LDX12 .LT. MAX( 1, M-Q ) ) THEN
         INFO = -9
      ELSE IF( COLMAJOR .AND. LDX21 .LT. MAX( 1, M-P ) ) THEN
         INFO = -11
      ELSE IF( .NOT.COLMAJOR .AND. LDX21 .LT. MAX( 1, Q ) ) THEN
         INFO = -11
      ELSE IF( COLMAJOR .AND. LDX22 .LT. MAX( 1, M-P ) ) THEN
         INFO = -13
      ELSE IF( .NOT.COLMAJOR .AND. LDX22 .LT. MAX( 1, M-Q ) ) THEN
         INFO = -13
      END IF
*
      // Compute workspace
*
      IF( INFO .EQ. 0 ) THEN
         LWORKOPT = M - Q
         LWORKMIN = M - Q
         WORK(1) = LWORKOPT
         IF( LWORK .LT. LWORKMIN .AND. .NOT. LQUERY ) THEN
            INFO = -21
         END IF
      END IF
      IF( INFO .NE. 0 ) THEN
         CALL XERBLA( 'xORBDB', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
      // Handle column-major and row-major separately
*
      IF( COLMAJOR ) THEN
*
         // Reduce columns 1, ..., Q of X11, X12, X21, and X22
*
         DO I = 1, Q
*
            IF( I .EQ. 1 ) THEN
               CALL SSCAL( P-I+1, Z1, X11(I,I), 1 )
            ELSE
               CALL SSCAL( P-I+1, Z1*COS(PHI(I-1)), X11(I,I), 1 )
               CALL SAXPY( P-I+1, -Z1*Z3*Z4*SIN(PHI(I-1)), X12(I,I-1), 1, X11(I,I), 1 )
            END IF
            IF( I .EQ. 1 ) THEN
               CALL SSCAL( M-P-I+1, Z2, X21(I,I), 1 )
            ELSE
               CALL SSCAL( M-P-I+1, Z2*COS(PHI(I-1)), X21(I,I), 1 )
               CALL SAXPY( M-P-I+1, -Z2*Z3*Z4*SIN(PHI(I-1)), X22(I,I-1), 1, X21(I,I), 1 )
            END IF
*
            THETA(I) = ATAN2( SNRM2( M-P-I+1, X21(I,I), 1 ), SNRM2( P-I+1, X11(I,I), 1 ) )
*
            IF( P .GT. I ) THEN
               CALL SLARFGP( P-I+1, X11(I,I), X11(I+1,I), 1, TAUP1(I) )
            ELSE IF( P .EQ. I ) THEN
               CALL SLARFGP( P-I+1, X11(I,I), X11(I,I), 1, TAUP1(I) )
            END IF
            X11(I,I) = ONE
            IF ( M-P .GT. I ) THEN
               CALL SLARFGP( M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) )
            ELSE IF ( M-P .EQ. I ) THEN
               CALL SLARFGP( M-P-I+1, X21(I,I), X21(I,I), 1, TAUP2(I) )
            END IF
            X21(I,I) = ONE
*
            IF ( Q .GT. I ) THEN
               CALL SLARF( 'L', P-I+1, Q-I, X11(I,I), 1, TAUP1(I), X11(I,I+1), LDX11, WORK )
            END IF
            IF ( M-Q+1 .GT. I ) THEN
               CALL SLARF( 'L', P-I+1, M-Q-I+1, X11(I,I), 1, TAUP1(I), X12(I,I), LDX12, WORK )
            END IF
            IF ( Q .GT. I ) THEN
               CALL SLARF( 'L', M-P-I+1, Q-I, X21(I,I), 1, TAUP2(I), X21(I,I+1), LDX21, WORK )
            END IF
            IF ( M-Q+1 .GT. I ) THEN
               CALL SLARF( 'L', M-P-I+1, M-Q-I+1, X21(I,I), 1, TAUP2(I), X22(I,I), LDX22, WORK )
            END IF
*
            IF( I .LT. Q ) THEN
               CALL SSCAL( Q-I, -Z1*Z3*SIN(THETA(I)), X11(I,I+1), LDX11 )                CALL SAXPY( Q-I, Z2*Z3*COS(THETA(I)), X21(I,I+1), LDX21, X11(I,I+1), LDX11 )
            END IF
            CALL SSCAL( M-Q-I+1, -Z1*Z4*SIN(THETA(I)), X12(I,I), LDX12 )
            CALL SAXPY( M-Q-I+1, Z2*Z4*COS(THETA(I)), X22(I,I), LDX22, X12(I,I), LDX12 )
*
            IF( I .LT. Q ) PHI(I) = ATAN2( SNRM2( Q-I, X11(I,I+1), LDX11 ), SNRM2( M-Q-I+1, X12(I,I), LDX12 ) )
*
            IF( I .LT. Q ) THEN
               IF ( Q-I .EQ. 1 ) THEN
                  CALL SLARFGP( Q-I, X11(I,I+1), X11(I,I+1), LDX11, TAUQ1(I) )
               ELSE
                  CALL SLARFGP( Q-I, X11(I,I+1), X11(I,I+2), LDX11, TAUQ1(I) )
               END IF
               X11(I,I+1) = ONE
            END IF
            IF ( Q+I-1 .LT. M ) THEN
               IF ( M-Q .EQ. I ) THEN
                  CALL SLARFGP( M-Q-I+1, X12(I,I), X12(I,I), LDX12, TAUQ2(I) )
               ELSE
                  CALL SLARFGP( M-Q-I+1, X12(I,I), X12(I,I+1), LDX12, TAUQ2(I) )
               END IF
            END IF
            X12(I,I) = ONE
*
            IF( I .LT. Q ) THEN
               CALL SLARF( 'R', P-I, Q-I, X11(I,I+1), LDX11, TAUQ1(I), X11(I+1,I+1), LDX11, WORK )                CALL SLARF( 'R', M-P-I, Q-I, X11(I,I+1), LDX11, TAUQ1(I), X21(I+1,I+1), LDX21, WORK )
            END IF
            IF ( P .GT. I ) THEN
               CALL SLARF( 'R', P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X12(I+1,I), LDX12, WORK )
            END IF
            IF ( M-P .GT. I ) THEN
               CALL SLARF( 'R', M-P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X22(I+1,I), LDX22, WORK )
            END IF
*
         END DO
*
         // Reduce columns Q + 1, ..., P of X12, X22
*
         DO I = Q + 1, P
*
            CALL SSCAL( M-Q-I+1, -Z1*Z4, X12(I,I), LDX12 )
            IF ( I .GE. M-Q ) THEN
               CALL SLARFGP( M-Q-I+1, X12(I,I), X12(I,I), LDX12, TAUQ2(I) )
            ELSE
               CALL SLARFGP( M-Q-I+1, X12(I,I), X12(I,I+1), LDX12, TAUQ2(I) )
            END IF
            X12(I,I) = ONE
*
            IF ( P .GT. I ) THEN
               CALL SLARF( 'R', P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X12(I+1,I), LDX12, WORK )
            END IF
            IF( M-P-Q .GE. 1 ) CALL SLARF( 'R', M-P-Q, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X22(Q+1,I), LDX22, WORK )
*
         END DO
*
         // Reduce columns P + 1, ..., M - Q of X12, X22
*
         DO I = 1, M - P - Q
*
            CALL SSCAL( M-P-Q-I+1, Z2*Z4, X22(Q+I,P+I), LDX22 )
            IF ( I .EQ. M-P-Q ) THEN
               CALL SLARFGP( M-P-Q-I+1, X22(Q+I,P+I), X22(Q+I,P+I), LDX22, TAUQ2(P+I) )
            ELSE
               CALL SLARFGP( M-P-Q-I+1, X22(Q+I,P+I), X22(Q+I,P+I+1), LDX22, TAUQ2(P+I) )
            END IF
            X22(Q+I,P+I) = ONE
            IF ( I .LT. M-P-Q ) THEN
               CALL SLARF( 'R', M-P-Q-I, M-P-Q-I+1, X22(Q+I,P+I), LDX22, TAUQ2(P+I), X22(Q+I+1,P+I), LDX22, WORK )
            END IF
*
         END DO
*
      ELSE
*
         // Reduce columns 1, ..., Q of X11, X12, X21, X22
*
         DO I = 1, Q
*
            IF( I .EQ. 1 ) THEN
               CALL SSCAL( P-I+1, Z1, X11(I,I), LDX11 )
            ELSE
               CALL SSCAL( P-I+1, Z1*COS(PHI(I-1)), X11(I,I), LDX11 )
               CALL SAXPY( P-I+1, -Z1*Z3*Z4*SIN(PHI(I-1)), X12(I-1,I), LDX12, X11(I,I), LDX11 )
            END IF
            IF( I .EQ. 1 ) THEN
               CALL SSCAL( M-P-I+1, Z2, X21(I,I), LDX21 )
            ELSE
               CALL SSCAL( M-P-I+1, Z2*COS(PHI(I-1)), X21(I,I), LDX21 )
               CALL SAXPY( M-P-I+1, -Z2*Z3*Z4*SIN(PHI(I-1)), X22(I-1,I), LDX22, X21(I,I), LDX21 )
            END IF
*
            THETA(I) = ATAN2( SNRM2( M-P-I+1, X21(I,I), LDX21 ), SNRM2( P-I+1, X11(I,I), LDX11 ) )
*
            CALL SLARFGP( P-I+1, X11(I,I), X11(I,I+1), LDX11, TAUP1(I) )
            X11(I,I) = ONE
            IF ( I .EQ. M-P ) THEN
               CALL SLARFGP( M-P-I+1, X21(I,I), X21(I,I), LDX21, TAUP2(I) )
            ELSE
               CALL SLARFGP( M-P-I+1, X21(I,I), X21(I,I+1), LDX21, TAUP2(I) )
            END IF
            X21(I,I) = ONE
*
            IF ( Q .GT. I ) THEN
               CALL SLARF( 'R', Q-I, P-I+1, X11(I,I), LDX11, TAUP1(I), X11(I+1,I), LDX11, WORK )
            END IF
            IF ( M-Q+1 .GT. I ) THEN
               CALL SLARF( 'R', M-Q-I+1, P-I+1, X11(I,I), LDX11, TAUP1(I), X12(I,I), LDX12, WORK )
            END IF
            IF ( Q .GT. I ) THEN
               CALL SLARF( 'R', Q-I, M-P-I+1, X21(I,I), LDX21, TAUP2(I), X21(I+1,I), LDX21, WORK )
            END IF
            IF ( M-Q+1 .GT. I ) THEN
               CALL SLARF( 'R', M-Q-I+1, M-P-I+1, X21(I,I), LDX21, TAUP2(I), X22(I,I), LDX22, WORK )
            END IF
*
            IF( I .LT. Q ) THEN
               CALL SSCAL( Q-I, -Z1*Z3*SIN(THETA(I)), X11(I+1,I), 1 )
               CALL SAXPY( Q-I, Z2*Z3*COS(THETA(I)), X21(I+1,I), 1, X11(I+1,I), 1 )
            END IF
            CALL SSCAL( M-Q-I+1, -Z1*Z4*SIN(THETA(I)), X12(I,I), 1 )
            CALL SAXPY( M-Q-I+1, Z2*Z4*COS(THETA(I)), X22(I,I), 1, X12(I,I), 1 )
*
            IF( I .LT. Q ) PHI(I) = ATAN2( SNRM2( Q-I, X11(I+1,I), 1 ), SNRM2( M-Q-I+1, X12(I,I), 1 ) )
*
            IF( I .LT. Q ) THEN
               IF ( Q-I .EQ. 1) THEN
                  CALL SLARFGP( Q-I, X11(I+1,I), X11(I+1,I), 1, TAUQ1(I) )
               ELSE
                  CALL SLARFGP( Q-I, X11(I+1,I), X11(I+2,I), 1, TAUQ1(I) )
               END IF
               X11(I+1,I) = ONE
            END IF
            IF ( M-Q .GT. I ) THEN
               CALL SLARFGP( M-Q-I+1, X12(I,I), X12(I+1,I), 1, TAUQ2(I) )
            ELSE
               CALL SLARFGP( M-Q-I+1, X12(I,I), X12(I,I), 1, TAUQ2(I) )
            END IF
            X12(I,I) = ONE
*
            IF( I .LT. Q ) THEN
               CALL SLARF( 'L', Q-I, P-I, X11(I+1,I), 1, TAUQ1(I), X11(I+1,I+1), LDX11, WORK )                CALL SLARF( 'L', Q-I, M-P-I, X11(I+1,I), 1, TAUQ1(I), X21(I+1,I+1), LDX21, WORK )
            END IF
            CALL SLARF( 'L', M-Q-I+1, P-I, X12(I,I), 1, TAUQ2(I), X12(I,I+1), LDX12, WORK )
            IF ( M-P-I .GT. 0 ) THEN
               CALL SLARF( 'L', M-Q-I+1, M-P-I, X12(I,I), 1, TAUQ2(I), X22(I,I+1), LDX22, WORK )
            END IF
*
         END DO
*
         // Reduce columns Q + 1, ..., P of X12, X22
*
         DO I = Q + 1, P
*
            CALL SSCAL( M-Q-I+1, -Z1*Z4, X12(I,I), 1 )
            CALL SLARFGP( M-Q-I+1, X12(I,I), X12(I+1,I), 1, TAUQ2(I) )
            X12(I,I) = ONE
*
            IF ( P .GT. I ) THEN
               CALL SLARF( 'L', M-Q-I+1, P-I, X12(I,I), 1, TAUQ2(I), X12(I,I+1), LDX12, WORK )
            END IF
            IF( M-P-Q .GE. 1 ) CALL SLARF( 'L', M-Q-I+1, M-P-Q, X12(I,I), 1, TAUQ2(I), X22(I,Q+1), LDX22, WORK )
*
         END DO
*
         // Reduce columns P + 1, ..., M - Q of X12, X22
*
         DO I = 1, M - P - Q
*
            CALL SSCAL( M-P-Q-I+1, Z2*Z4, X22(P+I,Q+I), 1 )
            IF ( M-P-Q .EQ. I ) THEN
               CALL SLARFGP( M-P-Q-I+1, X22(P+I,Q+I), X22(P+I,Q+I), 1, TAUQ2(P+I) )
               X22(P+I,Q+I) = ONE
            ELSE
               CALL SLARFGP( M-P-Q-I+1, X22(P+I,Q+I), X22(P+I+1,Q+I), 1, TAUQ2(P+I) )
               X22(P+I,Q+I) = ONE
               CALL SLARF( 'L', M-P-Q-I+1, M-P-Q-I, X22(P+I,Q+I), 1, TAUQ2(P+I), X22(P+I,Q+I+1), LDX22, WORK )
            END IF
*
*
         END DO
*
      END IF
*
      RETURN
*
      // End of SORBDB
*
      END

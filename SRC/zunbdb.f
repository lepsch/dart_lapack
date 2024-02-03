      SUBROUTINE ZUNBDB( TRANS, SIGNS, M, P, Q, X11, LDX11, X12, LDX12, X21, LDX21, X22, LDX22, THETA, PHI, TAUP1, TAUP2, TAUQ1, TAUQ2, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIGNS, TRANS;
      int                INFO, LDX11, LDX12, LDX21, LDX22, LWORK, M, P, Q;
      // ..
      // .. Array Arguments ..
      double             PHI( * ), THETA( * );
      COMPLEX*16         TAUP1( * ), TAUP2( * ), TAUQ1( * ), TAUQ2( * ), WORK( * ), X11( LDX11, * ), X12( LDX12, * ), X21( LDX21, * ), X22( LDX22, * )
      // ..

*  ====================================================================

      // .. Parameters ..
      double             REALONE;
      const              REALONE = 1.0D0 ;
      COMPLEX*16         ONE
      const              ONE = (1.0D0,0.0D0) ;
      // ..
      // .. Local Scalars ..
      bool               COLMAJOR, LQUERY;
      int                I, LWORKMIN, LWORKOPT;
      double             Z1, Z2, Z3, Z4;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZLARF, ZLARFGP, ZSCAL, XERBLA
      // EXTERNAL ZLACGV

      // ..
      // .. External Functions ..
      double             DZNRM2;
      bool               LSAME;
      // EXTERNAL DZNRM2, LSAME
      // ..
      // .. Intrinsic Functions
      // INTRINSIC ATAN2, COS, MAX, MIN, SIN
      // INTRINSIC DCMPLX, DCONJG
      // ..
      // .. Executable Statements ..

      // Test input arguments

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

      // Compute workspace

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

      // Handle column-major and row-major separately

      IF( COLMAJOR ) THEN

         // Reduce columns 1, ..., Q of X11, X12, X21, and X22

         DO I = 1, Q

            IF( I .EQ. 1 ) THEN
               CALL ZSCAL( P-I+1, DCMPLX( Z1, 0.0D0 ), X11(I,I), 1 )
            ELSE
               CALL ZSCAL( P-I+1, DCMPLX( Z1*COS(PHI(I-1)), 0.0D0 ), X11(I,I), 1 )                CALL ZAXPY( P-I+1, DCMPLX( -Z1*Z3*Z4*SIN(PHI(I-1)), 0.0D0 ), X12(I,I-1), 1, X11(I,I), 1 )
            END IF
            IF( I .EQ. 1 ) THEN
               CALL ZSCAL( M-P-I+1, DCMPLX( Z2, 0.0D0 ), X21(I,I), 1 )
            ELSE
               CALL ZSCAL( M-P-I+1, DCMPLX( Z2*COS(PHI(I-1)), 0.0D0 ), X21(I,I), 1 )                CALL ZAXPY( M-P-I+1, DCMPLX( -Z2*Z3*Z4*SIN(PHI(I-1)), 0.0D0 ), X22(I,I-1), 1, X21(I,I), 1 )
            END IF

            THETA(I) = ATAN2( DZNRM2( M-P-I+1, X21(I,I), 1 ), DZNRM2( P-I+1, X11(I,I), 1 ) )

            IF( P .GT. I ) THEN
               CALL ZLARFGP( P-I+1, X11(I,I), X11(I+1,I), 1, TAUP1(I) )
            ELSE IF ( P .EQ. I ) THEN
               CALL ZLARFGP( P-I+1, X11(I,I), X11(I,I), 1, TAUP1(I) )
            END IF
            X11(I,I) = ONE
            IF ( M-P .GT. I ) THEN
               CALL ZLARFGP( M-P-I+1, X21(I,I), X21(I+1,I), 1, TAUP2(I) )
            ELSE IF ( M-P .EQ. I ) THEN
               CALL ZLARFGP( M-P-I+1, X21(I,I), X21(I,I), 1, TAUP2(I) )
            END IF
            X21(I,I) = ONE

            IF ( Q .GT. I ) THEN
               CALL ZLARF( 'L', P-I+1, Q-I, X11(I,I), 1, DCONJG(TAUP1(I)), X11(I,I+1), LDX11, WORK )                CALL ZLARF( 'L', M-P-I+1, Q-I, X21(I,I), 1, DCONJG(TAUP2(I)), X21(I,I+1), LDX21, WORK )
            END IF
            IF ( M-Q+1 .GT. I ) THEN
               CALL ZLARF( 'L', P-I+1, M-Q-I+1, X11(I,I), 1, DCONJG(TAUP1(I)), X12(I,I), LDX12, WORK )                CALL ZLARF( 'L', M-P-I+1, M-Q-I+1, X21(I,I), 1, DCONJG(TAUP2(I)), X22(I,I), LDX22, WORK )
            END IF

            IF( I .LT. Q ) THEN
               CALL ZSCAL( Q-I, DCMPLX( -Z1*Z3*SIN(THETA(I)), 0.0D0 ), X11(I,I+1), LDX11 )                CALL ZAXPY( Q-I, DCMPLX( Z2*Z3*COS(THETA(I)), 0.0D0 ), X21(I,I+1), LDX21, X11(I,I+1), LDX11 )
            END IF
            CALL ZSCAL( M-Q-I+1, DCMPLX( -Z1*Z4*SIN(THETA(I)), 0.0D0 ), X12(I,I), LDX12 )             CALL ZAXPY( M-Q-I+1, DCMPLX( Z2*Z4*COS(THETA(I)), 0.0D0 ), X22(I,I), LDX22, X12(I,I), LDX12 )

            IF( I .LT. Q ) PHI(I) = ATAN2( DZNRM2( Q-I, X11(I,I+1), LDX11 ), DZNRM2( M-Q-I+1, X12(I,I), LDX12 ) )

            IF( I .LT. Q ) THEN
               CALL ZLACGV( Q-I, X11(I,I+1), LDX11 )
               IF ( I .EQ. Q-1 ) THEN
                  CALL ZLARFGP( Q-I, X11(I,I+1), X11(I,I+1), LDX11, TAUQ1(I) )
               ELSE
                  CALL ZLARFGP( Q-I, X11(I,I+1), X11(I,I+2), LDX11, TAUQ1(I) )
               END IF
               X11(I,I+1) = ONE
            END IF
            IF ( M-Q+1 .GT. I ) THEN
               CALL ZLACGV( M-Q-I+1, X12(I,I), LDX12 )
               IF ( M-Q .EQ. I ) THEN
                  CALL ZLARFGP( M-Q-I+1, X12(I,I), X12(I,I), LDX12, TAUQ2(I) )
               ELSE
                  CALL ZLARFGP( M-Q-I+1, X12(I,I), X12(I,I+1), LDX12, TAUQ2(I) )
               END IF
            END IF
            X12(I,I) = ONE

            IF( I .LT. Q ) THEN
               CALL ZLARF( 'R', P-I, Q-I, X11(I,I+1), LDX11, TAUQ1(I), X11(I+1,I+1), LDX11, WORK )                CALL ZLARF( 'R', M-P-I, Q-I, X11(I,I+1), LDX11, TAUQ1(I), X21(I+1,I+1), LDX21, WORK )
            END IF
            IF ( P .GT. I ) THEN
               CALL ZLARF( 'R', P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X12(I+1,I), LDX12, WORK )
            END IF
            IF ( M-P .GT. I ) THEN
               CALL ZLARF( 'R', M-P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X22(I+1,I), LDX22, WORK )
            END IF

            IF( I .LT. Q ) CALL ZLACGV( Q-I, X11(I,I+1), LDX11 )
            CALL ZLACGV( M-Q-I+1, X12(I,I), LDX12 )

         END DO

         // Reduce columns Q + 1, ..., P of X12, X22

         DO I = Q + 1, P

            CALL ZSCAL( M-Q-I+1, DCMPLX( -Z1*Z4, 0.0D0 ), X12(I,I), LDX12 )
            CALL ZLACGV( M-Q-I+1, X12(I,I), LDX12 )
            IF ( I .GE. M-Q ) THEN
               CALL ZLARFGP( M-Q-I+1, X12(I,I), X12(I,I), LDX12, TAUQ2(I) )
            ELSE
               CALL ZLARFGP( M-Q-I+1, X12(I,I), X12(I,I+1), LDX12, TAUQ2(I) )
            END IF
            X12(I,I) = ONE

            IF ( P .GT. I ) THEN
               CALL ZLARF( 'R', P-I, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X12(I+1,I), LDX12, WORK )
            END IF
            IF( M-P-Q .GE. 1 ) CALL ZLARF( 'R', M-P-Q, M-Q-I+1, X12(I,I), LDX12, TAUQ2(I), X22(Q+1,I), LDX22, WORK )

            CALL ZLACGV( M-Q-I+1, X12(I,I), LDX12 )

         END DO

         // Reduce columns P + 1, ..., M - Q of X12, X22

         DO I = 1, M - P - Q

            CALL ZSCAL( M-P-Q-I+1, DCMPLX( Z2*Z4, 0.0D0 ), X22(Q+I,P+I), LDX22 )
            CALL ZLACGV( M-P-Q-I+1, X22(Q+I,P+I), LDX22 )
            CALL ZLARFGP( M-P-Q-I+1, X22(Q+I,P+I), X22(Q+I,P+I+1), LDX22, TAUQ2(P+I) )
            X22(Q+I,P+I) = ONE
            CALL ZLARF( 'R', M-P-Q-I, M-P-Q-I+1, X22(Q+I,P+I), LDX22, TAUQ2(P+I), X22(Q+I+1,P+I), LDX22, WORK )

            CALL ZLACGV( M-P-Q-I+1, X22(Q+I,P+I), LDX22 )

         END DO

      ELSE

         // Reduce columns 1, ..., Q of X11, X12, X21, X22

         DO I = 1, Q

            IF( I .EQ. 1 ) THEN
               CALL ZSCAL( P-I+1, DCMPLX( Z1, 0.0D0 ), X11(I,I), LDX11 )
            ELSE
               CALL ZSCAL( P-I+1, DCMPLX( Z1*COS(PHI(I-1)), 0.0D0 ), X11(I,I), LDX11 )                CALL ZAXPY( P-I+1, DCMPLX( -Z1*Z3*Z4*SIN(PHI(I-1)), 0.0D0 ), X12(I-1,I), LDX12, X11(I,I), LDX11 )
            END IF
            IF( I .EQ. 1 ) THEN
               CALL ZSCAL( M-P-I+1, DCMPLX( Z2, 0.0D0 ), X21(I,I), LDX21 )
            ELSE
               CALL ZSCAL( M-P-I+1, DCMPLX( Z2*COS(PHI(I-1)), 0.0D0 ), X21(I,I), LDX21 )                CALL ZAXPY( M-P-I+1, DCMPLX( -Z2*Z3*Z4*SIN(PHI(I-1)), 0.0D0 ), X22(I-1,I), LDX22, X21(I,I), LDX21 )
            END IF

            THETA(I) = ATAN2( DZNRM2( M-P-I+1, X21(I,I), LDX21 ), DZNRM2( P-I+1, X11(I,I), LDX11 ) )

            CALL ZLACGV( P-I+1, X11(I,I), LDX11 )
            CALL ZLACGV( M-P-I+1, X21(I,I), LDX21 )

            CALL ZLARFGP( P-I+1, X11(I,I), X11(I,I+1), LDX11, TAUP1(I) )
            X11(I,I) = ONE
            IF ( I .EQ. M-P ) THEN
               CALL ZLARFGP( M-P-I+1, X21(I,I), X21(I,I), LDX21, TAUP2(I) )
            ELSE
               CALL ZLARFGP( M-P-I+1, X21(I,I), X21(I,I+1), LDX21, TAUP2(I) )
            END IF
            X21(I,I) = ONE

            CALL ZLARF( 'R', Q-I, P-I+1, X11(I,I), LDX11, TAUP1(I), X11(I+1,I), LDX11, WORK )             CALL ZLARF( 'R', M-Q-I+1, P-I+1, X11(I,I), LDX11, TAUP1(I), X12(I,I), LDX12, WORK )             CALL ZLARF( 'R', Q-I, M-P-I+1, X21(I,I), LDX21, TAUP2(I), X21(I+1,I), LDX21, WORK )             CALL ZLARF( 'R', M-Q-I+1, M-P-I+1, X21(I,I), LDX21, TAUP2(I), X22(I,I), LDX22, WORK )

            CALL ZLACGV( P-I+1, X11(I,I), LDX11 )
            CALL ZLACGV( M-P-I+1, X21(I,I), LDX21 )

            IF( I .LT. Q ) THEN
               CALL ZSCAL( Q-I, DCMPLX( -Z1*Z3*SIN(THETA(I)), 0.0D0 ), X11(I+1,I), 1 )                CALL ZAXPY( Q-I, DCMPLX( Z2*Z3*COS(THETA(I)), 0.0D0 ), X21(I+1,I), 1, X11(I+1,I), 1 )
            END IF
            CALL ZSCAL( M-Q-I+1, DCMPLX( -Z1*Z4*SIN(THETA(I)), 0.0D0 ), X12(I,I), 1 )             CALL ZAXPY( M-Q-I+1, DCMPLX( Z2*Z4*COS(THETA(I)), 0.0D0 ), X22(I,I), 1, X12(I,I), 1 )

            IF( I .LT. Q ) PHI(I) = ATAN2( DZNRM2( Q-I, X11(I+1,I), 1 ), DZNRM2( M-Q-I+1, X12(I,I), 1 ) )

            IF( I .LT. Q ) THEN
               CALL ZLARFGP( Q-I, X11(I+1,I), X11(I+2,I), 1, TAUQ1(I) )
               X11(I+1,I) = ONE
            END IF
            CALL ZLARFGP( M-Q-I+1, X12(I,I), X12(I+1,I), 1, TAUQ2(I) )
            X12(I,I) = ONE

            IF( I .LT. Q ) THEN
               CALL ZLARF( 'L', Q-I, P-I, X11(I+1,I), 1, DCONJG(TAUQ1(I)), X11(I+1,I+1), LDX11, WORK )                CALL ZLARF( 'L', Q-I, M-P-I, X11(I+1,I), 1, DCONJG(TAUQ1(I)), X21(I+1,I+1), LDX21, WORK )
            END IF
            CALL ZLARF( 'L', M-Q-I+1, P-I, X12(I,I), 1, DCONJG(TAUQ2(I)), X12(I,I+1), LDX12, WORK )
            IF ( M-P .GT. I ) THEN
               CALL ZLARF( 'L', M-Q-I+1, M-P-I, X12(I,I), 1, DCONJG(TAUQ2(I)), X22(I,I+1), LDX22, WORK )
            END IF

         END DO

         // Reduce columns Q + 1, ..., P of X12, X22

         DO I = Q + 1, P

            CALL ZSCAL( M-Q-I+1, DCMPLX( -Z1*Z4, 0.0D0 ), X12(I,I), 1 )
            CALL ZLARFGP( M-Q-I+1, X12(I,I), X12(I+1,I), 1, TAUQ2(I) )
            X12(I,I) = ONE

            IF ( P .GT. I ) THEN
               CALL ZLARF( 'L', M-Q-I+1, P-I, X12(I,I), 1, DCONJG(TAUQ2(I)), X12(I,I+1), LDX12, WORK )
            END IF
            IF( M-P-Q .GE. 1 ) CALL ZLARF( 'L', M-Q-I+1, M-P-Q, X12(I,I), 1, DCONJG(TAUQ2(I)), X22(I,Q+1), LDX22, WORK )

         END DO

         // Reduce columns P + 1, ..., M - Q of X12, X22

         DO I = 1, M - P - Q

            CALL ZSCAL( M-P-Q-I+1, DCMPLX( Z2*Z4, 0.0D0 ), X22(P+I,Q+I), 1 )             CALL ZLARFGP( M-P-Q-I+1, X22(P+I,Q+I), X22(P+I+1,Q+I), 1, TAUQ2(P+I) )
            X22(P+I,Q+I) = ONE

            IF ( M-P-Q .NE. I ) THEN
               CALL ZLARF( 'L', M-P-Q-I+1, M-P-Q-I, X22(P+I,Q+I), 1, DCONJG(TAUQ2(P+I)), X22(P+I,Q+I+1), LDX22, WORK )
            END IF

         END DO

      END IF

      RETURN

      // End of ZUNBDB

      }

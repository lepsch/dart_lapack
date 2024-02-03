      SUBROUTINE DCSDTS( M, P, Q, X, XF, LDX, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, THETA, IWORK, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDX, LDU1, LDU2, LDV1T, LDV2T, LWORK, M, P, Q;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             RESULT( 15 ), RWORK( * ), THETA( * );
      double             U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ), V2T( LDV2T, * ), WORK( LWORK ), X( LDX, * ), XF( LDX, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             REALONE, REALZERO;
      const              REALONE = 1.0D0, REALZERO = 0.0D0 ;
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      double             PIOVER2;
      const     PIOVER2 = 1.57079632679489661923132169163975144210D0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, R;
      double             EPS2, RESID, ULP, ULPINV;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE, DLANSY;
      // EXTERNAL DLAMCH, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DLASET, DORCSD, DORCSD2BY1, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC COS, DBLE, MAX, MIN, SIN
      // ..
      // .. Executable Statements ..

      ULP = DLAMCH( 'Precision' )
      ULPINV = REALONE / ULP

      // The first half of the routine checks the 2-by-2 CSD

      CALL DLASET( 'Full', M, M, ZERO, ONE, WORK, LDX )
      CALL DSYRK( 'Upper', 'Conjugate transpose', M, M, -ONE, X, LDX, ONE, WORK, LDX )
      IF (M.GT.0) THEN
         EPS2 = MAX( ULP, DLANGE( '1', M, M, WORK, LDX, RWORK ) / DBLE( M ) )
      } else {
         EPS2 = ULP
      END IF
      R = MIN( P, M-P, Q, M-Q )

      // Copy the matrix X to the array XF.

      CALL DLACPY( 'Full', M, M, X, LDX, XF, LDX )

      // Compute the CSD

      CALL DORCSD( 'Y', 'Y', 'Y', 'Y', 'N', 'D', M, P, Q, XF(1,1), LDX, XF(1,Q+1), LDX, XF(P+1,1), LDX, XF(P+1,Q+1), LDX, THETA, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, WORK, LWORK, IWORK, INFO )

      // Compute XF := diag(U1,U2)'*X*diag(V1,V2) - [D11 D12; D21 D22]

      CALL DLACPY( 'Full', M, M, X, LDX, XF, LDX )

      CALL DGEMM( 'No transpose', 'Conjugate transpose', P, Q, Q, ONE, XF, LDX, V1T, LDV1T, ZERO, WORK, LDX )

      CALL DGEMM( 'Conjugate transpose', 'No transpose', P, Q, P, ONE, U1, LDU1, WORK, LDX, ZERO, XF, LDX )

      DO I = 1, MIN(P,Q)-R
         XF(I,I) = XF(I,I) - ONE
      END DO
      DO I = 1, R
         XF(MIN(P,Q)-R+I,MIN(P,Q)-R+I) = XF(MIN(P,Q)-R+I,MIN(P,Q)-R+I) - COS(THETA(I))
      END DO

      CALL DGEMM( 'No transpose', 'Conjugate transpose', P, M-Q, M-Q, ONE, XF(1,Q+1), LDX, V2T, LDV2T, ZERO, WORK, LDX )

      CALL DGEMM( 'Conjugate transpose', 'No transpose', P, M-Q, P, ONE, U1, LDU1, WORK, LDX, ZERO, XF(1,Q+1), LDX )

      DO I = 1, MIN(P,M-Q)-R
         XF(P-I+1,M-I+1) = XF(P-I+1,M-I+1) + ONE
      END DO
      DO I = 1, R
         XF(P-(MIN(P,M-Q)-R)+1-I,M-(MIN(P,M-Q)-R)+1-I) = XF(P-(MIN(P,M-Q)-R)+1-I,M-(MIN(P,M-Q)-R)+1-I) + SIN(THETA(R-I+1))
      END DO

      CALL DGEMM( 'No transpose', 'Conjugate transpose', M-P, Q, Q, ONE, XF(P+1,1), LDX, V1T, LDV1T, ZERO, WORK, LDX )

      CALL DGEMM( 'Conjugate transpose', 'No transpose', M-P, Q, M-P, ONE, U2, LDU2, WORK, LDX, ZERO, XF(P+1,1), LDX )

      DO I = 1, MIN(M-P,Q)-R
         XF(M-I+1,Q-I+1) = XF(M-I+1,Q-I+1) - ONE
      END DO
      DO I = 1, R
         XF(M-(MIN(M-P,Q)-R)+1-I,Q-(MIN(M-P,Q)-R)+1-I) = XF(M-(MIN(M-P,Q)-R)+1-I,Q-(MIN(M-P,Q)-R)+1-I) - SIN(THETA(R-I+1))
      END DO

      CALL DGEMM( 'No transpose', 'Conjugate transpose', M-P, M-Q, M-Q, ONE, XF(P+1,Q+1), LDX, V2T, LDV2T, ZERO, WORK, LDX )

      CALL DGEMM( 'Conjugate transpose', 'No transpose', M-P, M-Q, M-P, ONE, U2, LDU2, WORK, LDX, ZERO, XF(P+1,Q+1), LDX )

      DO I = 1, MIN(M-P,M-Q)-R
         XF(P+I,Q+I) = XF(P+I,Q+I) - ONE
      END DO
      DO I = 1, R
         XF(P+(MIN(M-P,M-Q)-R)+I,Q+(MIN(M-P,M-Q)-R)+I) = XF(P+(MIN(M-P,M-Q)-R)+I,Q+(MIN(M-P,M-Q)-R)+I) - COS(THETA(I))
      END DO

      // Compute norm( U1'*X11*V1 - D11 ) / ( MAX(1,P,Q)*EPS2 ) .

      RESID = DLANGE( '1', P, Q, XF, LDX, RWORK )
      RESULT( 1 ) = ( RESID / DBLE(MAX(1,P,Q)) ) / EPS2

      // Compute norm( U1'*X12*V2 - D12 ) / ( MAX(1,P,M-Q)*EPS2 ) .

      RESID = DLANGE( '1', P, M-Q, XF(1,Q+1), LDX, RWORK )
      RESULT( 2 ) = ( RESID / DBLE(MAX(1,P,M-Q)) ) / EPS2

      // Compute norm( U2'*X21*V1 - D21 ) / ( MAX(1,M-P,Q)*EPS2 ) .

      RESID = DLANGE( '1', M-P, Q, XF(P+1,1), LDX, RWORK )
      RESULT( 3 ) = ( RESID / DBLE(MAX(1,M-P,Q)) ) / EPS2

      // Compute norm( U2'*X22*V2 - D22 ) / ( MAX(1,M-P,M-Q)*EPS2 ) .

      RESID = DLANGE( '1', M-P, M-Q, XF(P+1,Q+1), LDX, RWORK )
      RESULT( 4 ) = ( RESID / DBLE(MAX(1,M-P,M-Q)) ) / EPS2

      // Compute I - U1'*U1

      CALL DLASET( 'Full', P, P, ZERO, ONE, WORK, LDU1 )
      CALL DSYRK( 'Upper', 'Conjugate transpose', P, P, -ONE, U1, LDU1, ONE, WORK, LDU1 )

      // Compute norm( I - U'*U ) / ( MAX(1,P) * ULP ) .

      RESID = DLANSY( '1', 'Upper', P, WORK, LDU1, RWORK )
      RESULT( 5 ) = ( RESID / DBLE(MAX(1,P)) ) / ULP

      // Compute I - U2'*U2

      CALL DLASET( 'Full', M-P, M-P, ZERO, ONE, WORK, LDU2 )
      CALL DSYRK( 'Upper', 'Conjugate transpose', M-P, M-P, -ONE, U2, LDU2, ONE, WORK, LDU2 )

      // Compute norm( I - U2'*U2 ) / ( MAX(1,M-P) * ULP ) .

      RESID = DLANSY( '1', 'Upper', M-P, WORK, LDU2, RWORK )
      RESULT( 6 ) = ( RESID / DBLE(MAX(1,M-P)) ) / ULP

      // Compute I - V1T*V1T'

      CALL DLASET( 'Full', Q, Q, ZERO, ONE, WORK, LDV1T )
      CALL DSYRK( 'Upper', 'No transpose', Q, Q, -ONE, V1T, LDV1T, ONE, WORK, LDV1T )

      // Compute norm( I - V1T*V1T' ) / ( MAX(1,Q) * ULP ) .

      RESID = DLANSY( '1', 'Upper', Q, WORK, LDV1T, RWORK )
      RESULT( 7 ) = ( RESID / DBLE(MAX(1,Q)) ) / ULP

      // Compute I - V2T*V2T'

      CALL DLASET( 'Full', M-Q, M-Q, ZERO, ONE, WORK, LDV2T )
      CALL DSYRK( 'Upper', 'No transpose', M-Q, M-Q, -ONE, V2T, LDV2T, ONE, WORK, LDV2T )

      // Compute norm( I - V2T*V2T' ) / ( MAX(1,M-Q) * ULP ) .

      RESID = DLANSY( '1', 'Upper', M-Q, WORK, LDV2T, RWORK )
      RESULT( 8 ) = ( RESID / DBLE(MAX(1,M-Q)) ) / ULP

      // Check sorting

      RESULT( 9 ) = REALZERO
      DO I = 1, R
         IF( THETA(I).LT.REALZERO .OR. THETA(I).GT.PIOVER2 ) THEN
            RESULT( 9 ) = ULPINV
         END IF
         IF( I.GT.1 ) THEN
            IF ( THETA(I).LT.THETA(I-1) ) THEN
               RESULT( 9 ) = ULPINV
            END IF
         END IF
      END DO

      // The second half of the routine checks the 2-by-1 CSD

      CALL DLASET( 'Full', Q, Q, ZERO, ONE, WORK, LDX )
      CALL DSYRK( 'Upper', 'Conjugate transpose', Q, M, -ONE, X, LDX, ONE, WORK, LDX )
      IF( M.GT.0 ) THEN
         EPS2 = MAX( ULP, DLANGE( '1', Q, Q, WORK, LDX, RWORK ) / DBLE( M ) )
      } else {
         EPS2 = ULP
      END IF
      R = MIN( P, M-P, Q, M-Q )

      // Copy the matrix [ X11; X21 ] to the array XF.

      CALL DLACPY( 'Full', M, Q, X, LDX, XF, LDX )

      // Compute the CSD

      CALL DORCSD2BY1( 'Y', 'Y', 'Y', M, P, Q, XF(1,1), LDX, XF(P+1,1), LDX, THETA, U1, LDU1, U2, LDU2, V1T, LDV1T, WORK, LWORK, IWORK, INFO )

      // Compute [X11;X21] := diag(U1,U2)'*[X11;X21]*V1 - [D11;D21]

      CALL DGEMM( 'No transpose', 'Conjugate transpose', P, Q, Q, ONE, X, LDX, V1T, LDV1T, ZERO, WORK, LDX )

      CALL DGEMM( 'Conjugate transpose', 'No transpose', P, Q, P, ONE, U1, LDU1, WORK, LDX, ZERO, X, LDX )

      DO I = 1, MIN(P,Q)-R
         X(I,I) = X(I,I) - ONE
      END DO
      DO I = 1, R
         X(MIN(P,Q)-R+I,MIN(P,Q)-R+I) = X(MIN(P,Q)-R+I,MIN(P,Q)-R+I) - COS(THETA(I))
      END DO

      CALL DGEMM( 'No transpose', 'Conjugate transpose', M-P, Q, Q, ONE, X(P+1,1), LDX, V1T, LDV1T, ZERO, WORK, LDX )

      CALL DGEMM( 'Conjugate transpose', 'No transpose', M-P, Q, M-P, ONE, U2, LDU2, WORK, LDX, ZERO, X(P+1,1), LDX )

      DO I = 1, MIN(M-P,Q)-R
         X(M-I+1,Q-I+1) = X(M-I+1,Q-I+1) - ONE
      END DO
      DO I = 1, R
         X(M-(MIN(M-P,Q)-R)+1-I,Q-(MIN(M-P,Q)-R)+1-I) = X(M-(MIN(M-P,Q)-R)+1-I,Q-(MIN(M-P,Q)-R)+1-I) - SIN(THETA(R-I+1))
      END DO

      // Compute norm( U1'*X11*V1 - D11 ) / ( MAX(1,P,Q)*EPS2 ) .

      RESID = DLANGE( '1', P, Q, X, LDX, RWORK )
      RESULT( 10 ) = ( RESID / DBLE(MAX(1,P,Q)) ) / EPS2

      // Compute norm( U2'*X21*V1 - D21 ) / ( MAX(1,M-P,Q)*EPS2 ) .

      RESID = DLANGE( '1', M-P, Q, X(P+1,1), LDX, RWORK )
      RESULT( 11 ) = ( RESID / DBLE(MAX(1,M-P,Q)) ) / EPS2

      // Compute I - U1'*U1

      CALL DLASET( 'Full', P, P, ZERO, ONE, WORK, LDU1 )
      CALL DSYRK( 'Upper', 'Conjugate transpose', P, P, -ONE, U1, LDU1, ONE, WORK, LDU1 )

      // Compute norm( I - U1'*U1 ) / ( MAX(1,P) * ULP ) .

      RESID = DLANSY( '1', 'Upper', P, WORK, LDU1, RWORK )
      RESULT( 12 ) = ( RESID / DBLE(MAX(1,P)) ) / ULP

      // Compute I - U2'*U2

      CALL DLASET( 'Full', M-P, M-P, ZERO, ONE, WORK, LDU2 )
      CALL DSYRK( 'Upper', 'Conjugate transpose', M-P, M-P, -ONE, U2, LDU2, ONE, WORK, LDU2 )

      // Compute norm( I - U2'*U2 ) / ( MAX(1,M-P) * ULP ) .

      RESID = DLANSY( '1', 'Upper', M-P, WORK, LDU2, RWORK )
      RESULT( 13 ) = ( RESID / DBLE(MAX(1,M-P)) ) / ULP

      // Compute I - V1T*V1T'

      CALL DLASET( 'Full', Q, Q, ZERO, ONE, WORK, LDV1T )
      CALL DSYRK( 'Upper', 'No transpose', Q, Q, -ONE, V1T, LDV1T, ONE, WORK, LDV1T )

      // Compute norm( I - V1T*V1T' ) / ( MAX(1,Q) * ULP ) .

      RESID = DLANSY( '1', 'Upper', Q, WORK, LDV1T, RWORK )
      RESULT( 14 ) = ( RESID / DBLE(MAX(1,Q)) ) / ULP

      // Check sorting

      RESULT( 15 ) = REALZERO
      DO I = 1, R
         IF( THETA(I).LT.REALZERO .OR. THETA(I).GT.PIOVER2 ) THEN
            RESULT( 15 ) = ULPINV
         END IF
         IF( I.GT.1 ) THEN
            IF ( THETA(I).LT.THETA(I-1) ) THEN
               RESULT( 15 ) = ULPINV
            END IF
         END IF
      END DO

      RETURN

      // End of DCSDTS

      }

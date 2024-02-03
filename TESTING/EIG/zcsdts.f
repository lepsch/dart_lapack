      void zcsdts(M, P, Q, X, XF, LDX, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, THETA, IWORK, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDX, LDU1, LDU2, LDV1T, LDV2T, LWORK, M, P, Q;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             RESULT( 15 ), RWORK( * ), THETA( * );
      Complex         U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ), V2T( LDV2T, * ), WORK( LWORK ), X( LDX, * ), XF( LDX, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             REALONE, REALZERO;
      const              REALONE = 1.0, REALZERO = 0.0 ;
      Complex         ZERO, ONE;
      const              ZERO = (0.0,0.0), ONE = (1.0,0.0) ;
      double             PIOVER2;
      const     PIOVER2 = 1.57079632679489661923132169163975144210 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, R;
      double             EPS2, RESID, ULP, ULPINV;
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE, ZLANHE;
      // EXTERNAL DLAMCH, ZLANGE, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZHERK, ZLACPY, ZLASET, ZUNCSD, ZUNCSD2BY1
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC COS, DBLE, DCMPLX, MAX, MIN, REAL, SIN
      // ..
      // .. Executable Statements ..

      ULP = DLAMCH( 'Precision' );
      ULPINV = REALONE / ULP;

      // The first half of the routine checks the 2-by-2 CSD

      zlaset('Full', M, M, ZERO, ONE, WORK, LDX );
      zherk('Upper', 'Conjugate transpose', M, M, -REALONE, X, LDX, REALONE, WORK, LDX );
      if (M > 0) {
         EPS2 = max( ULP, ZLANGE( '1', M, M, WORK, LDX, RWORK ) / DBLE( M ) );
      } else {
         EPS2 = ULP;
      }
      R = min( P, M-P, Q, M-Q );

      // Copy the matrix X to the array XF.

      zlacpy('Full', M, M, X, LDX, XF, LDX );

      // Compute the CSD

      zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'D', M, P, Q, XF(1,1), LDX, XF(1,Q+1), LDX, XF(P+1,1), LDX, XF(P+1,Q+1), LDX, THETA, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, WORK, LWORK, RWORK, 17*(R+2), IWORK, INFO );

      // Compute XF := diag(U1,U2)'*X*diag(V1,V2) - [D11 D12; D21 D22]

      zlacpy('Full', M, M, X, LDX, XF, LDX );

      zgemm('No transpose', 'Conjugate transpose', P, Q, Q, ONE, XF, LDX, V1T, LDV1T, ZERO, WORK, LDX );

      zgemm('Conjugate transpose', 'No transpose', P, Q, P, ONE, U1, LDU1, WORK, LDX, ZERO, XF, LDX );

      DO I = 1, min(P,Q)-R;
         XF(I,I) = XF(I,I) - ONE;
      }
      for (I = 1; I <= R; I++) {
         XF(min(P,Q)-R+I,min(P,Q)-R+I) = XF(min(P,Q)-R+I,min(P,Q)-R+I) - DCMPLX( COS(THETA(I)), 0.0 );
      }

      zgemm('No transpose', 'Conjugate transpose', P, M-Q, M-Q, ONE, XF(1,Q+1), LDX, V2T, LDV2T, ZERO, WORK, LDX );

      zgemm('Conjugate transpose', 'No transpose', P, M-Q, P, ONE, U1, LDU1, WORK, LDX, ZERO, XF(1,Q+1), LDX );

      DO I = 1, min(P,M-Q)-R;
         XF(P-I+1,M-I+1) = XF(P-I+1,M-I+1) + ONE;
      }
      for (I = 1; I <= R; I++) {
         XF(P-(min(P,M-Q)-R)+1-I,M-(min(P,M-Q)-R)+1-I) = XF(P-(min(P,M-Q)-R)+1-I,M-(min(P,M-Q)-R)+1-I) + DCMPLX( SIN(THETA(R-I+1)), 0.0 );
      }

      zgemm('No transpose', 'Conjugate transpose', M-P, Q, Q, ONE, XF(P+1,1), LDX, V1T, LDV1T, ZERO, WORK, LDX );

      zgemm('Conjugate transpose', 'No transpose', M-P, Q, M-P, ONE, U2, LDU2, WORK, LDX, ZERO, XF(P+1,1), LDX );

      DO I = 1, min(M-P,Q)-R;
         XF(M-I+1,Q-I+1) = XF(M-I+1,Q-I+1) - ONE;
      }
      for (I = 1; I <= R; I++) {
         XF(M-(min(M-P,Q)-R)+1-I,Q-(min(M-P,Q)-R)+1-I) = XF(M-(min(M-P,Q)-R)+1-I,Q-(min(M-P,Q)-R)+1-I) - DCMPLX( SIN(THETA(R-I+1)), 0.0 );
      }

      zgemm('No transpose', 'Conjugate transpose', M-P, M-Q, M-Q, ONE, XF(P+1,Q+1), LDX, V2T, LDV2T, ZERO, WORK, LDX );

      zgemm('Conjugate transpose', 'No transpose', M-P, M-Q, M-P, ONE, U2, LDU2, WORK, LDX, ZERO, XF(P+1,Q+1), LDX );

      DO I = 1, min(M-P,M-Q)-R;
         XF(P+I,Q+I) = XF(P+I,Q+I) - ONE;
      }
      for (I = 1; I <= R; I++) {
         XF(P+(min(M-P,M-Q)-R)+I,Q+(min(M-P,M-Q)-R)+I) = XF(P+(min(M-P,M-Q)-R)+I,Q+(min(M-P,M-Q)-R)+I) - DCMPLX( COS(THETA(I)), 0.0 );
      }

      // Compute norm( U1'*X11*V1 - D11 ) / ( max(1,P,Q)*EPS2 ) .

      RESID = ZLANGE( '1', P, Q, XF, LDX, RWORK );
      RESULT( 1 ) = ( RESID / REAL(max(1,P,Q)) ) / EPS2;

      // Compute norm( U1'*X12*V2 - D12 ) / ( max(1,P,M-Q)*EPS2 ) .

      RESID = ZLANGE( '1', P, M-Q, XF(1,Q+1), LDX, RWORK );
      RESULT( 2 ) = ( RESID / REAL(max(1,P,M-Q)) ) / EPS2;

      // Compute norm( U2'*X21*V1 - D21 ) / ( max(1,M-P,Q)*EPS2 ) .

      RESID = ZLANGE( '1', M-P, Q, XF(P+1,1), LDX, RWORK );
      RESULT( 3 ) = ( RESID / REAL(max(1,M-P,Q)) ) / EPS2;

      // Compute norm( U2'*X22*V2 - D22 ) / ( max(1,M-P,M-Q)*EPS2 ) .

      RESID = ZLANGE( '1', M-P, M-Q, XF(P+1,Q+1), LDX, RWORK );
      RESULT( 4 ) = ( RESID / REAL(max(1,M-P,M-Q)) ) / EPS2;

      // Compute I - U1'*U1

      zlaset('Full', P, P, ZERO, ONE, WORK, LDU1 );
      zherk('Upper', 'Conjugate transpose', P, P, -REALONE, U1, LDU1, REALONE, WORK, LDU1 );

      // Compute norm( I - U'*U ) / ( max(1,P) * ULP ) .

      RESID = ZLANHE( '1', 'Upper', P, WORK, LDU1, RWORK );
      RESULT( 5 ) = ( RESID / REAL(max(1,P)) ) / ULP;

      // Compute I - U2'*U2

      zlaset('Full', M-P, M-P, ZERO, ONE, WORK, LDU2 );
      zherk('Upper', 'Conjugate transpose', M-P, M-P, -REALONE, U2, LDU2, REALONE, WORK, LDU2 );

      // Compute norm( I - U2'*U2 ) / ( max(1,M-P) * ULP ) .

      RESID = ZLANHE( '1', 'Upper', M-P, WORK, LDU2, RWORK );
      RESULT( 6 ) = ( RESID / REAL(max(1,M-P)) ) / ULP;

      // Compute I - V1T*V1T'

      zlaset('Full', Q, Q, ZERO, ONE, WORK, LDV1T );
      zherk('Upper', 'No transpose', Q, Q, -REALONE, V1T, LDV1T, REALONE, WORK, LDV1T );

      // Compute norm( I - V1T*V1T' ) / ( max(1,Q) * ULP ) .

      RESID = ZLANHE( '1', 'Upper', Q, WORK, LDV1T, RWORK );
      RESULT( 7 ) = ( RESID / REAL(max(1,Q)) ) / ULP;

      // Compute I - V2T*V2T'

      zlaset('Full', M-Q, M-Q, ZERO, ONE, WORK, LDV2T );
      zherk('Upper', 'No transpose', M-Q, M-Q, -REALONE, V2T, LDV2T, REALONE, WORK, LDV2T );

      // Compute norm( I - V2T*V2T' ) / ( max(1,M-Q) * ULP ) .

      RESID = ZLANHE( '1', 'Upper', M-Q, WORK, LDV2T, RWORK );
      RESULT( 8 ) = ( RESID / REAL(max(1,M-Q)) ) / ULP;

      // Check sorting

      RESULT( 9 ) = REALZERO;
      for (I = 1; I <= R; I++) {
         if ( THETA(I) < REALZERO || THETA(I) > PIOVER2 ) {
            RESULT( 9 ) = ULPINV;
         }
         if ( I > 1) {
            if ( THETA(I) < THETA(I-1) ) {
               RESULT( 9 ) = ULPINV;
            }
         }
      }

      // The second half of the routine checks the 2-by-1 CSD

      zlaset('Full', Q, Q, ZERO, ONE, WORK, LDX );
      zherk('Upper', 'Conjugate transpose', Q, M, -REALONE, X, LDX, REALONE, WORK, LDX );
      if (M > 0) {
         EPS2 = max( ULP, ZLANGE( '1', Q, Q, WORK, LDX, RWORK ) / DBLE( M ) );
      } else {
         EPS2 = ULP;
      }
      R = min( P, M-P, Q, M-Q );

      // Copy the matrix X to the array XF.

      zlacpy('Full', M, M, X, LDX, XF, LDX );

      // Compute the CSD

      zuncsd2by1('Y', 'Y', 'Y', M, P, Q, XF(1,1), LDX, XF(P+1,1), LDX, THETA, U1, LDU1, U2, LDU2, V1T, LDV1T, WORK, LWORK, RWORK, 17*(R+2), IWORK, INFO );

      // Compute [X11;X21] := diag(U1,U2)'*[X11;X21]*V1 - [D11;D21]

      zgemm('No transpose', 'Conjugate transpose', P, Q, Q, ONE, X, LDX, V1T, LDV1T, ZERO, WORK, LDX );

      zgemm('Conjugate transpose', 'No transpose', P, Q, P, ONE, U1, LDU1, WORK, LDX, ZERO, X, LDX );

      DO I = 1, min(P,Q)-R;
         X(I,I) = X(I,I) - ONE;
      }
      for (I = 1; I <= R; I++) {
         X(min(P,Q)-R+I,min(P,Q)-R+I) = X(min(P,Q)-R+I,min(P,Q)-R+I) - DCMPLX( COS(THETA(I)), 0.0 );
      }

      zgemm('No transpose', 'Conjugate transpose', M-P, Q, Q, ONE, X(P+1,1), LDX, V1T, LDV1T, ZERO, WORK, LDX );

      zgemm('Conjugate transpose', 'No transpose', M-P, Q, M-P, ONE, U2, LDU2, WORK, LDX, ZERO, X(P+1,1), LDX );

      DO I = 1, min(M-P,Q)-R;
         X(M-I+1,Q-I+1) = X(M-I+1,Q-I+1) - ONE;
      }
      for (I = 1; I <= R; I++) {
         X(M-(min(M-P,Q)-R)+1-I,Q-(min(M-P,Q)-R)+1-I) = X(M-(min(M-P,Q)-R)+1-I,Q-(min(M-P,Q)-R)+1-I) - DCMPLX( SIN(THETA(R-I+1)), 0.0 );
      }

      // Compute norm( U1'*X11*V1 - D11 ) / ( max(1,P,Q)*EPS2 ) .

      RESID = ZLANGE( '1', P, Q, X, LDX, RWORK );
      RESULT( 10 ) = ( RESID / REAL(max(1,P,Q)) ) / EPS2;

      // Compute norm( U2'*X21*V1 - D21 ) / ( max(1,M-P,Q)*EPS2 ) .

      RESID = ZLANGE( '1', M-P, Q, X(P+1,1), LDX, RWORK );
      RESULT( 11 ) = ( RESID / REAL(max(1,M-P,Q)) ) / EPS2;

      // Compute I - U1'*U1

      zlaset('Full', P, P, ZERO, ONE, WORK, LDU1 );
      zherk('Upper', 'Conjugate transpose', P, P, -REALONE, U1, LDU1, REALONE, WORK, LDU1 );

      // Compute norm( I - U'*U ) / ( max(1,P) * ULP ) .

      RESID = ZLANHE( '1', 'Upper', P, WORK, LDU1, RWORK );
      RESULT( 12 ) = ( RESID / REAL(max(1,P)) ) / ULP;

      // Compute I - U2'*U2

      zlaset('Full', M-P, M-P, ZERO, ONE, WORK, LDU2 );
      zherk('Upper', 'Conjugate transpose', M-P, M-P, -REALONE, U2, LDU2, REALONE, WORK, LDU2 );

      // Compute norm( I - U2'*U2 ) / ( max(1,M-P) * ULP ) .

      RESID = ZLANHE( '1', 'Upper', M-P, WORK, LDU2, RWORK );
      RESULT( 13 ) = ( RESID / REAL(max(1,M-P)) ) / ULP;

      // Compute I - V1T*V1T'

      zlaset('Full', Q, Q, ZERO, ONE, WORK, LDV1T );
      zherk('Upper', 'No transpose', Q, Q, -REALONE, V1T, LDV1T, REALONE, WORK, LDV1T );

      // Compute norm( I - V1T*V1T' ) / ( max(1,Q) * ULP ) .

      RESID = ZLANHE( '1', 'Upper', Q, WORK, LDV1T, RWORK );
      RESULT( 14 ) = ( RESID / REAL(max(1,Q)) ) / ULP;

      // Check sorting

      RESULT( 15 ) = REALZERO;
      for (I = 1; I <= R; I++) {
         if ( THETA(I) < REALZERO || THETA(I) > PIOVER2 ) {
            RESULT( 15 ) = ULPINV;
         }
         if ( I > 1) {
            if ( THETA(I) < THETA(I-1) ) {
               RESULT( 15 ) = ULPINV;
            }
         }
      }

      return;
      }

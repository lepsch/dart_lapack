      void ccsdts(final int M, final int P, final int Q, final int X, final int XF, final int LDX, final Matrix<double> U1_, final int LDU1, final Matrix<double> U2_, final int LDU2, final Matrix<double> V1T_, final int LDV1T, final Matrix<double> V2T_, final int LDV2T, final int THETA, final Array<int> IWORK_, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final int RESULT,) {
  final U1 = U1_.dim();
  final U2 = U2_.dim();
  final V1T = V1T_.dim();
  final V2T = V2T_.dim();
  final IWORK = IWORK_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDX, LDU1, LDU2, LDV1T, LDV2T, LWORK, M, P, Q;
      int                IWORK( * );
      double               RESULT( 15 ), RWORK( * ), THETA( * );
      Complex            U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ), V2T( LDV2T, * ), WORK( LWORK ), X( LDX, * ), XF( LDX, * );
      // ..

      double               REALONE, REALZERO;
      const              REALONE = 1.0, REALZERO = 0.0 ;
      Complex            ZERO, ONE;
      const              ZERO = (0.0,0.0), ONE = (1.0,0.0) ;
      double               PIOVER2;
      const     PIOVER2 = 1.57079632679489661923132169163975144210 ;
      int                I, INFO, R;
      double               EPS2, RESID, ULP, ULPINV;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, CLANGE, CLANHE;
      // EXTERNAL SLAMCH, CLANGE, CLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHERK, CLACPY, CLASET, CUNCSD, CUNCSD2BY1
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, COS, MAX, MIN, REAL, SIN

      ULP = SLAMCH( 'Precision' );
      ULPINV = REALONE / ULP;

      // The first half of the routine checks the 2-by-2 CSD

      claset('Full', M, M, ZERO, ONE, WORK, LDX );
      cherk('Upper', 'Conjugate transpose', M, M, -REALONE, X, LDX, REALONE, WORK, LDX );
      if (M > 0) {
         EPS2 = max( ULP, CLANGE( '1', M, M, WORK, LDX, RWORK ) / REAL( M ) );
      } else {
         EPS2 = ULP;
      }
      R = min( P, M-P, Q, M-Q );

      // Copy the matrix X to the array XF.

      clacpy('Full', M, M, X, LDX, XF, LDX );

      // Compute the CSD

      cuncsd('Y', 'Y', 'Y', 'Y', 'N', 'D', M, P, Q, XF(1,1), LDX, XF(1,Q+1), LDX, XF(P+1,1), LDX, XF(P+1,Q+1), LDX, THETA, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, WORK, LWORK, RWORK, 17*(R+2), IWORK, INFO );

      // Compute XF := diag(U1,U2)'*X*diag(V1,V2) - [D11 D12; D21 D22]

      clacpy('Full', M, M, X, LDX, XF, LDX );

      cgemm('No transpose', 'Conjugate transpose', P, Q, Q, ONE, XF, LDX, V1T, LDV1T, ZERO, WORK, LDX );

      cgemm('Conjugate transpose', 'No transpose', P, Q, P, ONE, U1, LDU1, WORK, LDX, ZERO, XF, LDX );

      for (I = 1; I <= min(P,Q)-R; I++) {
         XF[I][I] = XF(I,I) - ONE;
      }
      for (I = 1; I <= R; I++) {
         XF[min(P,Q)-R+I,min(P,Q)-R+I] = XF(min(P,Q)-R+I,min(P,Q)-R+I) - CMPLX( COS(THETA(I)), 0.0 );
      }

      cgemm('No transpose', 'Conjugate transpose', P, M-Q, M-Q, ONE, XF(1,Q+1), LDX, V2T, LDV2T, ZERO, WORK, LDX );

      cgemm('Conjugate transpose', 'No transpose', P, M-Q, P, ONE, U1, LDU1, WORK, LDX, ZERO, XF(1,Q+1), LDX );

      for (I = 1; I <= min(P,M-Q)-R; I++) {
         XF[P-I+1][M-I+1] = XF(P-I+1,M-I+1) + ONE;
      }
      for (I = 1; I <= R; I++) {
         XF[P-(min(P,M-Q)-R)+1-I,M-(min(P,M-Q)-R)+1-I] = XF(P-(min(P,M-Q)-R)+1-I,M-(min(P,M-Q)-R)+1-I) + CMPLX( SIN(THETA(R-I+1)), 0.0 );
      }

      cgemm('No transpose', 'Conjugate transpose', M-P, Q, Q, ONE, XF(P+1,1), LDX, V1T, LDV1T, ZERO, WORK, LDX );

      cgemm('Conjugate transpose', 'No transpose', M-P, Q, M-P, ONE, U2, LDU2, WORK, LDX, ZERO, XF(P+1,1), LDX );

      for (I = 1; I <= min(M-P,Q)-R; I++) {
         XF[M-I+1][Q-I+1] = XF(M-I+1,Q-I+1) - ONE;
      }
      for (I = 1; I <= R; I++) {
         XF[M-(min(M-P,Q)-R)+1-I,Q-(min(M-P,Q)-R)+1-I] = XF(M-(min(M-P,Q)-R)+1-I,Q-(min(M-P,Q)-R)+1-I) - CMPLX( SIN(THETA(R-I+1)), 0.0 );
      }

      cgemm('No transpose', 'Conjugate transpose', M-P, M-Q, M-Q, ONE, XF(P+1,Q+1), LDX, V2T, LDV2T, ZERO, WORK, LDX );

      cgemm('Conjugate transpose', 'No transpose', M-P, M-Q, M-P, ONE, U2, LDU2, WORK, LDX, ZERO, XF(P+1,Q+1), LDX );

      for (I = 1; I <= min(M-P,M-Q)-R; I++) {
         XF[P+I][Q+I] = XF(P+I,Q+I) - ONE;
      }
      for (I = 1; I <= R; I++) {
         XF[P+(min(M-P,M-Q)-R)+I,Q+(min(M-P,M-Q)-R)+I] = XF(P+(min(M-P,M-Q)-R)+I,Q+(min(M-P,M-Q)-R)+I) - CMPLX( COS(THETA(I)), 0.0 );
      }

      // Compute norm( U1'*X11*V1 - D11 ) / ( max(1,P,Q)*EPS2 ) .

      RESID = CLANGE( '1', P, Q, XF, LDX, RWORK );
      RESULT[1] = ( RESID / REAL(max(1,P,Q)) ) / EPS2;

      // Compute norm( U1'*X12*V2 - D12 ) / ( max(1,P,M-Q)*EPS2 ) .

      RESID = CLANGE( '1', P, M-Q, XF(1,Q+1), LDX, RWORK );
      RESULT[2] = ( RESID / REAL(max(1,P,M-Q)) ) / EPS2;

      // Compute norm( U2'*X21*V1 - D21 ) / ( max(1,M-P,Q)*EPS2 ) .

      RESID = CLANGE( '1', M-P, Q, XF(P+1,1), LDX, RWORK );
      RESULT[3] = ( RESID / REAL(max(1,M-P,Q)) ) / EPS2;

      // Compute norm( U2'*X22*V2 - D22 ) / ( max(1,M-P,M-Q)*EPS2 ) .

      RESID = CLANGE( '1', M-P, M-Q, XF(P+1,Q+1), LDX, RWORK );
      RESULT[4] = ( RESID / REAL(max(1,M-P,M-Q)) ) / EPS2;

      // Compute I - U1'*U1

      claset('Full', P, P, ZERO, ONE, WORK, LDU1 );
      cherk('Upper', 'Conjugate transpose', P, P, -REALONE, U1, LDU1, REALONE, WORK, LDU1 );

      // Compute norm( I - U'*U ) / ( max(1,P) * ULP ) .

      RESID = CLANHE( '1', 'Upper', P, WORK, LDU1, RWORK );
      RESULT[5] = ( RESID / REAL(max(1,P)) ) / ULP;

      // Compute I - U2'*U2

      claset('Full', M-P, M-P, ZERO, ONE, WORK, LDU2 );
      cherk('Upper', 'Conjugate transpose', M-P, M-P, -REALONE, U2, LDU2, REALONE, WORK, LDU2 );

      // Compute norm( I - U2'*U2 ) / ( max(1,M-P) * ULP ) .

      RESID = CLANHE( '1', 'Upper', M-P, WORK, LDU2, RWORK );
      RESULT[6] = ( RESID / REAL(max(1,M-P)) ) / ULP;

      // Compute I - V1T*V1T'

      claset('Full', Q, Q, ZERO, ONE, WORK, LDV1T );
      cherk('Upper', 'No transpose', Q, Q, -REALONE, V1T, LDV1T, REALONE, WORK, LDV1T );

      // Compute norm( I - V1T*V1T' ) / ( max(1,Q) * ULP ) .

      RESID = CLANHE( '1', 'Upper', Q, WORK, LDV1T, RWORK );
      RESULT[7] = ( RESID / REAL(max(1,Q)) ) / ULP;

      // Compute I - V2T*V2T'

      claset('Full', M-Q, M-Q, ZERO, ONE, WORK, LDV2T );
      cherk('Upper', 'No transpose', M-Q, M-Q, -REALONE, V2T, LDV2T, REALONE, WORK, LDV2T );

      // Compute norm( I - V2T*V2T' ) / ( max(1,M-Q) * ULP ) .

      RESID = CLANHE( '1', 'Upper', M-Q, WORK, LDV2T, RWORK );
      RESULT[8] = ( RESID / REAL(max(1,M-Q)) ) / ULP;

      // Check sorting

      RESULT[9] = REALZERO;
      for (I = 1; I <= R; I++) {
         if ( THETA(I) < REALZERO || THETA(I) > PIOVER2 ) {
            RESULT[9] = ULPINV;
         }
         if ( I > 1) {
            if ( THETA(I) < THETA(I-1) ) {
               RESULT[9] = ULPINV;
            }
         }
      }

      // The second half of the routine checks the 2-by-1 CSD

      claset('Full', Q, Q, ZERO, ONE, WORK, LDX );
      cherk('Upper', 'Conjugate transpose', Q, M, -REALONE, X, LDX, REALONE, WORK, LDX );
      if (M > 0) {
         EPS2 = max( ULP, CLANGE( '1', Q, Q, WORK, LDX, RWORK ) / REAL( M ) );
      } else {
         EPS2 = ULP;
      }
      R = min( P, M-P, Q, M-Q );

      // Copy the matrix X to the array XF.

      clacpy('Full', M, Q, X, LDX, XF, LDX );

      // Compute the CSD

      cuncsd2by1('Y', 'Y', 'Y', M, P, Q, XF(1,1), LDX, XF(P+1,1), LDX, THETA, U1, LDU1, U2, LDU2, V1T, LDV1T, WORK, LWORK, RWORK, 17*(R+2), IWORK, INFO );

      // Compute [X11;X21] := diag(U1,U2)'*[X11;X21]*V1 - [D11;D21]

      cgemm('No transpose', 'Conjugate transpose', P, Q, Q, ONE, X, LDX, V1T, LDV1T, ZERO, WORK, LDX );

      cgemm('Conjugate transpose', 'No transpose', P, Q, P, ONE, U1, LDU1, WORK, LDX, ZERO, X, LDX );

      for (I = 1; I <= min(P,Q)-R; I++) {
         X[I][I] = X(I,I) - ONE;
      }
      for (I = 1; I <= R; I++) {
         X[min(P,Q)-R+I,min(P,Q)-R+I] = X(min(P,Q)-R+I,min(P,Q)-R+I) - CMPLX( COS(THETA(I)), 0.0 );
      }

      cgemm('No transpose', 'Conjugate transpose', M-P, Q, Q, ONE, X(P+1,1), LDX, V1T, LDV1T, ZERO, WORK, LDX );

      cgemm('Conjugate transpose', 'No transpose', M-P, Q, M-P, ONE, U2, LDU2, WORK, LDX, ZERO, X(P+1,1), LDX );

      for (I = 1; I <= min(M-P,Q)-R; I++) {
         X[M-I+1][Q-I+1] = X(M-I+1,Q-I+1) - ONE;
      }
      for (I = 1; I <= R; I++) {
         X[M-(min(M-P,Q)-R)+1-I,Q-(min(M-P,Q)-R)+1-I] = X(M-(min(M-P,Q)-R)+1-I,Q-(min(M-P,Q)-R)+1-I) - CMPLX( SIN(THETA(R-I+1)), 0.0 );
      }

      // Compute norm( U1'*X11*V1 - D11 ) / ( max(1,P,Q)*EPS2 ) .

      RESID = CLANGE( '1', P, Q, X, LDX, RWORK );
      RESULT[10] = ( RESID / REAL(max(1,P,Q)) ) / EPS2;

      // Compute norm( U2'*X21*V1 - D21 ) / ( max(1,M-P,Q)*EPS2 ) .

      RESID = CLANGE( '1', M-P, Q, X(P+1,1), LDX, RWORK );
      RESULT[11] = ( RESID / REAL(max(1,M-P,Q)) ) / EPS2;

      // Compute I - U1'*U1

      claset('Full', P, P, ZERO, ONE, WORK, LDU1 );
      cherk('Upper', 'Conjugate transpose', P, P, -REALONE, U1, LDU1, REALONE, WORK, LDU1 );

      // Compute norm( I - U1'*U1 ) / ( max(1,P) * ULP ) .

      RESID = CLANHE( '1', 'Upper', P, WORK, LDU1, RWORK );
      RESULT[12] = ( RESID / REAL(max(1,P)) ) / ULP;

      // Compute I - U2'*U2

      claset('Full', M-P, M-P, ZERO, ONE, WORK, LDU2 );
      cherk('Upper', 'Conjugate transpose', M-P, M-P, -REALONE, U2, LDU2, REALONE, WORK, LDU2 );

      // Compute norm( I - U2'*U2 ) / ( max(1,M-P) * ULP ) .

      RESID = CLANHE( '1', 'Upper', M-P, WORK, LDU2, RWORK );
      RESULT[13] = ( RESID / REAL(max(1,M-P)) ) / ULP;

      // Compute I - V1T*V1T'

      claset('Full', Q, Q, ZERO, ONE, WORK, LDV1T );
      cherk('Upper', 'No transpose', Q, Q, -REALONE, V1T, LDV1T, REALONE, WORK, LDV1T );

      // Compute norm( I - V1T*V1T' ) / ( max(1,Q) * ULP ) .

      RESID = CLANHE( '1', 'Upper', Q, WORK, LDV1T, RWORK );
      RESULT[14] = ( RESID / REAL(max(1,Q)) ) / ULP;

      // Check sorting

      RESULT[15] = REALZERO;
      for (I = 1; I <= R; I++) {
         if ( THETA(I) < REALZERO || THETA(I) > PIOVER2 ) {
            RESULT[15] = ULPINV;
         }
         if ( I > 1) {
            if ( THETA(I) < THETA(I-1) ) {
               RESULT[15] = ULPINV;
            }
         }
      }

      }
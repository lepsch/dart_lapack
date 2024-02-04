      void cunt01(ROWCOL, M, N, U, LDU, WORK, LWORK, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             ROWCOL;
      int                LDU, LWORK, M, N;
      REAL               RESID;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * );
      COMPLEX            U( LDU, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      String             TRANSU;
      int                I, J, K, LDWORK, MNMIN;
      REAL               EPS;
      COMPLEX            TMP, ZDUM;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- REAL               CLANSY, SLAMCH;
      //- COMPLEX            CDOTC;
      // EXTERNAL LSAME, CLANSY, SLAMCH, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHERK, CLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, CMPLX, MAX, MIN, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( REAL( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      RESID = ZERO;

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      EPS = SLAMCH( 'Precision' );
      if ( M < N || ( M == N && LSAME( ROWCOL, 'R' ) ) ) {
         TRANSU = 'N';
         K = N;
      } else {
         TRANSU = 'C';
         K = M;
      }
      MNMIN = min( M, N );

      if ( ( MNMIN+1 )*MNMIN <= LWORK ) {
         LDWORK = MNMIN;
      } else {
         LDWORK = 0;
      }
      if ( LDWORK > 0 ) {

         // Compute I - U*U' or I - U'*U.

         claset('Upper', MNMIN, MNMIN, CMPLX( ZERO ), CMPLX( ONE ), WORK, LDWORK );
         cherk('Upper', TRANSU, MNMIN, K, -ONE, U, LDU, ONE, WORK, LDWORK );

         // Compute norm( I - U*U' ) / ( K * EPS ) .

         RESID = CLANSY( '1', 'Upper', MNMIN, WORK, LDWORK, RWORK );
         RESID = ( RESID / REAL( K ) ) / EPS;
      } else if ( TRANSU == 'C' ) {

         // Find the maximum element in abs( I - U'*U ) / ( m * EPS )

         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               if ( I != J ) {
                  TMP = ZERO;
               } else {
                  TMP = ONE;
               }
               TMP = TMP - CDOTC( M, U( 1, I ), 1, U( 1, J ), 1 );
               RESID = max( RESID, CABS1( TMP ) );
            } // 10
         } // 20
         RESID = ( RESID / REAL( M ) ) / EPS;
      } else {

         // Find the maximum element in abs( I - U*U' ) / ( n * EPS )

         for (J = 1; J <= M; J++) { // 40
            for (I = 1; I <= J; I++) { // 30
               if ( I != J ) {
                  TMP = ZERO;
               } else {
                  TMP = ONE;
               }
               TMP = TMP - CDOTC( N, U( J, 1 ), LDU, U( I, 1 ), LDU );
               RESID = max( RESID, CABS1( TMP ) );
            } // 30
         } // 40
         RESID = ( RESID / REAL( N ) ) / EPS;
      }
      return;
      }

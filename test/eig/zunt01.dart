      void zunt01(ROWCOL, M, N, U, LDU, WORK, LWORK, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             ROWCOL;
      int                LDU, LWORK, M, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      Complex         U( LDU, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      String             TRANSU;
      int                I, J, K, LDWORK, MNMIN;
      double             EPS;
      Complex         TMP, ZDUM;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, ZLANSY;
      //- Complex         ZDOTC;
      // EXTERNAL lsame, DLAMCH, ZLANSY, ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZHERK, ZLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DIMAG, MAX, MIN
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[ZDUM] = ( ZDUM.toDouble() ).abs() + ( DIMAG( ZDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      RESID = ZERO;

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      EPS = DLAMCH( 'Precision' );
      if ( M < N || ( M == N && lsame( ROWCOL, 'R' ) ) ) {
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

         zlaset('Upper', MNMIN, MNMIN, DCMPLX( ZERO ), DCMPLX( ONE ), WORK, LDWORK );
         zherk('Upper', TRANSU, MNMIN, K, -ONE, U, LDU, ONE, WORK, LDWORK );

         // Compute norm( I - U*U' ) / ( K * EPS ) .

         RESID = ZLANSY( '1', 'Upper', MNMIN, WORK, LDWORK, RWORK );
         RESID = ( RESID / K.toDouble() ) / EPS;
      } else if ( TRANSU == 'C' ) {

         // Find the maximum element in abs( I - U'*U ) / ( m * EPS )

         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               if ( I != J ) {
                  TMP = ZERO;
               } else {
                  TMP = ONE;
               }
               TMP = TMP - ZDOTC( M, U( 1, I ), 1, U( 1, J ), 1 );
               RESID = max( RESID, CABS1( TMP ) );
            } // 10
         } // 20
         RESID = ( RESID / M.toDouble() ) / EPS;
      } else {

         // Find the maximum element in abs( I - U*U' ) / ( n * EPS )

         for (J = 1; J <= M; J++) { // 40
            for (I = 1; I <= J; I++) { // 30
               if ( I != J ) {
                  TMP = ZERO;
               } else {
                  TMP = ONE;
               }
               TMP = TMP - ZDOTC( N, U( J, 1 ), LDU, U( I, 1 ), LDU );
               RESID = max( RESID, CABS1( TMP ) );
            } // 30
         } // 40
         RESID = ( RESID / N.toDouble() ) / EPS;
      }
      return;
      }
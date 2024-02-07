      void dort01(ROWCOL, M, N, U, LDU, WORK, LWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             ROWCOL;
      int                LDU, LWORK, M, N;
      double             RESID;
      double             U( LDU, * ), WORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      String             TRANSU;
      int                I, J, K, LDWORK, MNMIN;
      double             EPS, TMP;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DDOT, DLAMCH, DLANSY;
      // EXTERNAL lsame, DDOT, DLAMCH, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASET, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN

      RESID = ZERO;

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      EPS = dlamch( 'Precision' );
      if ( M < N || ( M == N && lsame( ROWCOL, 'R' ) ) ) {
         TRANSU = 'N';
         K = N;
      } else {
         TRANSU = 'T';
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

         dlaset('Upper', MNMIN, MNMIN, ZERO, ONE, WORK, LDWORK );
         dsyrk('Upper', TRANSU, MNMIN, K, -ONE, U, LDU, ONE, WORK, LDWORK );

         // Compute norm( I - U*U' ) / ( K * EPS ) .

         RESID = dlansy( '1', 'Upper', MNMIN, WORK, LDWORK, WORK( LDWORK*MNMIN+1 ) );
         RESID = ( RESID / K.toDouble() ) / EPS;
      } else if ( TRANSU == 'T' ) {

         // Find the maximum element in abs( I - U'*U ) / ( m * EPS )

         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               if ( I != J ) {
                  TMP = ZERO;
               } else {
                  TMP = ONE;
               }
               TMP = TMP - ddot( M, U( 1, I ), 1, U( 1, J ), 1 );
               RESID = max( RESID, ( TMP ).abs() );
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
               TMP = TMP - ddot( N, U( J, 1 ), LDU, U( I, 1 ), LDU );
               RESID = max( RESID, ( TMP ).abs() );
            } // 30
         } // 40
         RESID = ( RESID / N.toDouble() ) / EPS;
      }
      return;
      }

      SUBROUTINE ZUNT01( ROWCOL, M, N, U, LDU, WORK, LWORK, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             ROWCOL;
      int                LDU, LWORK, M, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             RWORK( * );
      COMPLEX*16         U( LDU, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      String             TRANSU;
      int                I, J, K, LDWORK, MNMIN;
      double             EPS;
      COMPLEX*16         TMP, ZDUM
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANSY;
      COMPLEX*16         ZDOTC
      // EXTERNAL LSAME, DLAMCH, ZLANSY, ZDOTC
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
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      RESID = ZERO

      // Quick return if possible

      if (M.LE.0 .OR. N.LE.0) RETURN;

      EPS = DLAMCH( 'Precision' )
      if ( M.LT.N .OR. ( M == N .AND. LSAME( ROWCOL, 'R' ) ) ) {
         TRANSU = 'N'
         K = N
      } else {
         TRANSU = 'C'
         K = M
      }
      MNMIN = MIN( M, N )

      if ( ( MNMIN+1 )*MNMIN.LE.LWORK ) {
         LDWORK = MNMIN
      } else {
         LDWORK = 0
      }
      if ( LDWORK.GT.0 ) {

         // Compute I - U*U' or I - U'*U.

         zlaset('Upper', MNMIN, MNMIN, DCMPLX( ZERO ), DCMPLX( ONE ), WORK, LDWORK );
         zherk('Upper', TRANSU, MNMIN, K, -ONE, U, LDU, ONE, WORK, LDWORK );

         // Compute norm( I - U*U' ) / ( K * EPS ) .

         RESID = ZLANSY( '1', 'Upper', MNMIN, WORK, LDWORK, RWORK )
         RESID = ( RESID / DBLE( K ) ) / EPS
      } else if ( TRANSU == 'C' ) {

         // Find the maximum element in abs( I - U'*U ) / ( m * EPS )

         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               if ( I.NE.J ) {
                  TMP = ZERO
               } else {
                  TMP = ONE
               }
               TMP = TMP - ZDOTC( M, U( 1, I ), 1, U( 1, J ), 1 )
               RESID = MAX( RESID, CABS1( TMP ) )
            } // 10
         } // 20
         RESID = ( RESID / DBLE( M ) ) / EPS
      } else {

         // Find the maximum element in abs( I - U*U' ) / ( n * EPS )

         for (J = 1; J <= M; J++) { // 40
            for (I = 1; I <= J; I++) { // 30
               if ( I.NE.J ) {
                  TMP = ZERO
               } else {
                  TMP = ONE
               }
               TMP = TMP - ZDOTC( N, U( J, 1 ), LDU, U( I, 1 ), LDU )
               RESID = MAX( RESID, CABS1( TMP ) )
            } // 30
         } // 40
         RESID = ( RESID / DBLE( N ) ) / EPS
      }
      RETURN

      // End of ZUNT01

      }

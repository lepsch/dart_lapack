      void zspt01(UPLO, N, A, AFAC, IPIV, C, LDC, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDC, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double             RWORK( * );
      Complex         A( * ), AFAC( * ), C( LDC, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, JC;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, ZLANSP, ZLANSY;
      // EXTERNAL lsame, DLAMCH, ZLANSP, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASET, ZLAVSP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Determine EPS and the norm of A.

      EPS = dlamch( 'Epsilon' );
      ANORM = ZLANSP( '1', UPLO, N, A, RWORK );

      // Initialize C to the identity matrix.

      zlaset('Full', N, N, CZERO, CONE, C, LDC );

      // Call ZLAVSP to form the product D * U' (or D * L' ).

      zlavsp(UPLO, 'Transpose', 'Non-unit', N, N, AFAC, IPIV, C, LDC, INFO );

      // Call ZLAVSP again to multiply by U ( or L ).

      zlavsp(UPLO, 'No transpose', 'Unit', N, N, AFAC, IPIV, C, LDC, INFO );

      // Compute the difference  C - A .

      if ( lsame( UPLO, 'U' ) ) {
         JC = 0;
         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               C[I, J] = C( I, J ) - A( JC+I );
            } // 10
            JC = JC + J;
         } // 20
      } else {
         JC = 1;
         for (J = 1; J <= N; J++) { // 40
            for (I = J; I <= N; I++) { // 30
               C[I, J] = C( I, J ) - A( JC+I-J );
            } // 30
            JC = JC + N - J + 1;
         } // 40
      }

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = ZLANSY( '1', UPLO, N, C, LDC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / N.toDouble() ) / ANORM ) / EPS;
      }

      return;
      }

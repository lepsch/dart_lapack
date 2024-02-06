      void sspt01(UPLO, N, A, AFAC, IPIV, C, LDC, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDC, N;
      double               RESID;
      int                IPIV( * );
      double               A( * ), AFAC( * ), C( LDC, * ), RWORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, INFO, J, JC;
      double               ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH, SLANSP, SLANSY;
      // EXTERNAL lsame, SLAMCH, SLANSP, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLAVSP, SLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Determine EPS and the norm of A.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = SLANSP( '1', UPLO, N, A, RWORK );

      // Initialize C to the identity matrix.

      slaset('Full', N, N, ZERO, ONE, C, LDC );

      // Call SLAVSP to form the product D * U' (or D * L' ).

      slavsp(UPLO, 'Transpose', 'Non-unit', N, N, AFAC, IPIV, C, LDC, INFO );

      // Call SLAVSP again to multiply by U ( or L ).

      slavsp(UPLO, 'No transpose', 'Unit', N, N, AFAC, IPIV, C, LDC, INFO );

      // Compute the difference  C - A .

      if ( lsame( UPLO, 'U' ) ) {
         JC = 0;
         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               C[I][J] = C( I, J ) - A( JC+I );
            } // 10
            JC = JC + J;
         } // 20
      } else {
         JC = 1;
         for (J = 1; J <= N; J++) { // 40
            for (I = J; I <= N; I++) { // 30
               C[I][J] = C( I, J ) - A( JC+I-J );
            } // 30
            JC = JC + N - J + 1;
         } // 40
      }

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = SLANSY( '1', UPLO, N, C, LDC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS;
      }

      return;
      }

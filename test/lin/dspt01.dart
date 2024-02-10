      void dspt01(final int UPLO, final int N, final int A, final int AFAC, final Array<int> IPIV, final Matrix<double> C, final int LDC, final Array<double> RWORK, final int RESID) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDC, N;
      double             RESID;
      int                IPIV( * );
      double             A( * ), AFAC( * ), C( LDC, * ), RWORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, INFO, J, JC;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DLANSP, DLANSY;
      // EXTERNAL lsame, DLAMCH, DLANSP, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASET, DLAVSP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Determine EPS and the norm of A.

      EPS = dlamch( 'Epsilon' );
      ANORM = dlansp( '1', UPLO, N, A, RWORK );

      // Initialize C to the identity matrix.

      dlaset('Full', N, N, ZERO, ONE, C, LDC );

      // Call DLAVSP to form the product D * U' (or D * L' ).

      dlavsp(UPLO, 'Transpose', 'Non-unit', N, N, AFAC, IPIV, C, LDC, INFO );

      // Call DLAVSP again to multiply by U ( or L ).

      dlavsp(UPLO, 'No transpose', 'Unit', N, N, AFAC, IPIV, C, LDC, INFO );

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

      RESID = dlansy( '1', UPLO, N, C, LDC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / N.toDouble() ) / ANORM ) / EPS;
      }

      }

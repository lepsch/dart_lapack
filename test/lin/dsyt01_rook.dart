      void dsyt01_rook(final int UPLO, final int N, final Matrix<double> A, final int LDA, final Matrix<double> AFAC, final int LDAFAC, final Array<int> IPIV, final Matrix<double> C, final int LDC, final Array<double> RWORK, final int RESID) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      double             RESID;
      int                IPIV( * );
      double             A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ), RWORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, INFO, J;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DLANSY;
      // EXTERNAL lsame, DLAMCH, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASET, DLAVSY_ROOK
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
      ANORM = dlansy( '1', UPLO, N, A, LDA, RWORK );

      // Initialize C to the identity matrix.

      dlaset('Full', N, N, ZERO, ONE, C, LDC );

      // Call DLAVSY_ROOK to form the product D * U' (or D * L' ).

      dlavsy_rook(UPLO, 'Transpose', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // Call DLAVSY_ROOK again to multiply by U (or L ).

      dlavsy_rook(UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // Compute the difference  C - A .

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 20
            for (I = 1; I <= J; I++) { // 10
               C[I][J] = C( I, J ) - A( I, J );
            } // 10
         } // 20
      } else {
         for (J = 1; J <= N; J++) { // 40
            for (I = J; I <= N; I++) { // 30
               C[I][J] = C( I, J ) - A( I, J );
            } // 30
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

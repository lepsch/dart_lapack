      void zhet01_rook(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> AFAC_, final int LDAFAC, final Array<int> IPIV_, final Matrix<double> C_, final int LDC, final Array<double> RWORK_, final int RESID,) {
  final A = A_.having();
  final AFAC = AFAC_.having();
  final IPIV = IPIV_.having();
  final C = C_.having();
  final RWORK = RWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      double             RESID;
      int                IPIV( * );
      double             RWORK( * );
      Complex         A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                I, INFO, J;
      double             ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             ZLANHE, DLAMCH;
      // EXTERNAL lsame, ZLANHE, DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASET, ZLAVHE_ROOK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DIMAG, DBLE

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Determine EPS and the norm of A.

      EPS = dlamch( 'Epsilon' );
      ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK );

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      for (J = 1; J <= N; J++) { // 10
         if ( DIMAG( AFAC( J, J ) ) != ZERO ) {
            RESID = ONE / EPS;
            return;
         }
      } // 10

      // Initialize C to the identity matrix.

      zlaset('Full', N, N, CZERO, CONE, C, LDC );

      // Call ZLAVHE_ROOK to form the product D * U' (or D * L' ).

      zlavhe_rook(UPLO, 'Conjugate', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // Call ZLAVHE_ROOK again to multiply by U (or L ).

      zlavhe_rook(UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // Compute the difference  C - A .

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 30
            for (I = 1; I <= J - 1; I++) { // 20
               C[I][J] = C( I, J ) - A( I, J );
            } // 20
            C[J][J] = C( J, J ) - (A( J, J )).toDouble();
         } // 30
      } else {
         for (J = 1; J <= N; J++) { // 50
            C[J][J] = C( J, J ) - (A( J, J )).toDouble();
            for (I = J + 1; I <= N; I++) { // 40
               C[I][J] = C( I, J ) - A( I, J );
            } // 40
         } // 50
      }

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = ZLANHE( '1', UPLO, N, C, LDC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID/N.toDouble() )/ANORM ) / EPS;
      }

      }

      void chet01_3(UPLO, N, final Matrix<double> A, final int LDA, final Matrix<double> AFAC, final int LDAFAC, E, IPIV, final Matrix<double> C, final int LDC, final Array<double> RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      double               RESID;
      int                IPIV( * );
      double               RWORK( * );
      Complex            A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ), E( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                I, INFO, J;
      double               ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANHE, SLAMCH;
      // EXTERNAL lsame, CLANHE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASET, CLAVHE_ROOK, CSYCONVF_ROOK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, REAL

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // a) Revert to multipliers of L

      csyconvf_rook(UPLO, 'R', N, AFAC, LDAFAC, E, IPIV, INFO );

      // 1) Determine EPS and the norm of A.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK );

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      for (J = 1; J <= N; J++) {
         if ( AIMAG( AFAC( J, J ) ) != ZERO ) {
            RESID = ONE / EPS;
            return;
         }
      }

      // 2) Initialize C to the identity matrix.

      claset('Full', N, N, CZERO, CONE, C, LDC );

      // 3) Call CLAVHE_ROOK to form the product D * U' (or D * L' ).

      clavhe_rook(UPLO, 'Conjugate', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // 4) Call ZLAVHE_RK again to multiply by U (or L ).

      clavhe_rook(UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // 5) Compute the difference  C - A .

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= J - 1; I++) {
               C[I][J] = C( I, J ) - A( I, J );
            }
            C[J][J] = C( J, J ) - double( A( J, J ) );
         }
      } else {
         for (J = 1; J <= N; J++) {
            C[J][J] = C( J, J ) - double( A( J, J ) );
            for (I = J + 1; I <= N; I++) {
               C[I][J] = C( I, J ) - A( I, J );
            }
         }
      }

      // 6) Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = CLANHE( '1', UPLO, N, C, LDC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID/REAL( N ) )/ANORM ) / EPS;
      }

      // b) Convert to factor of L (or U)

      csyconvf_rook(UPLO, 'C', N, AFAC, LDAFAC, E, IPIV, INFO );

      }

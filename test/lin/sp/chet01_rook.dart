      void chet01_rook(UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      double               RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      double               RWORK( * );
      Complex            A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J;
      double               ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANHE, SLAMCH;
      // EXTERNAL lsame, CLANHE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASET, CLAVHE_ROOK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Determine EPS and the norm of A.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK );

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      for (J = 1; J <= N; J++) { // 10
         if ( AIMAG( AFAC( J, J ) ) != ZERO ) {
            RESID = ONE / EPS;
            return;
         }
      } // 10

      // Initialize C to the identity matrix.

      claset('Full', N, N, CZERO, CONE, C, LDC );

      // Call CLAVHE_ROOK to form the product D * U' (or D * L' ).

      clavhe_rook(UPLO, 'Conjugate', 'Non-unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // Call CLAVHE_ROOK again to multiply by U (or L ).

      clavhe_rook(UPLO, 'No transpose', 'Unit', N, N, AFAC, LDAFAC, IPIV, C, LDC, INFO );

      // Compute the difference  C - A .

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 30
            for (I = 1; I <= J - 1; I++) { // 20
               C[I, J] = C( I, J ) - A( I, J );
            } // 20
            C[J, J] = C( J, J ) - double( A( J, J ) );
         } // 30
      } else {
         for (J = 1; J <= N; J++) { // 50
            C[J, J] = C( J, J ) - double( A( J, J ) );
            for (I = J + 1; I <= N; I++) { // 40
               C[I, J] = C( I, J ) - A( I, J );
            } // 40
         } // 50
      }

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = CLANHE( '1', UPLO, N, C, LDC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID/REAL( N ) )/ANORM ) / EPS;
      }

      return;
      }

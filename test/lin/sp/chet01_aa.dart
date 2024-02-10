      void chet01_aa(final int UPLO, final int N, final Matrix<double> A, final int LDA, final Matrix<double> AFAC, final int LDAFAC, final Array<int> IPIV, final Matrix<double> C, final int LDC, final Array<double> RWORK, final int RESID) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      double               RESID;
      int                IPIV( * );
      double               RWORK( * );
      Complex            A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * );
      // ..

      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE  = ( 1.0, 0.0 ) ;
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, J;
      double               ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH, CLANHE;
      // EXTERNAL lsame, SLAMCH, CLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASET, CLAVHE
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
      ANORM = CLANHE( '1', UPLO, N, A, LDA, RWORK );

      // Initialize C to the tridiagonal matrix T.

      claset('Full', N, N, CZERO, CZERO, C, LDC );
      clacpy('F', 1, N, AFAC( 1, 1 ), LDAFAC+1, C( 1, 1 ), LDC+1 );
      if ( N > 1 ) {
         if ( lsame( UPLO, 'U' ) ) {
            clacpy('F', 1, N-1, AFAC( 1, 2 ), LDAFAC+1, C( 1, 2 ), LDC+1 );
            clacpy('F', 1, N-1, AFAC( 1, 2 ), LDAFAC+1, C( 2, 1 ), LDC+1 );
            clacgv(N-1, C( 2, 1 ), LDC+1 );
         } else {
            clacpy('F', 1, N-1, AFAC( 2, 1 ), LDAFAC+1, C( 1, 2 ), LDC+1 );
            clacpy('F', 1, N-1, AFAC( 2, 1 ), LDAFAC+1, C( 2, 1 ), LDC+1 );
            clacgv(N-1, C( 1, 2 ), LDC+1 );
         }

         // Call CTRMM to form the product U' * D (or L * D ).

         if ( lsame( UPLO, 'U' ) ) {
            ctrmm('Left', UPLO, 'Conjugate transpose', 'Unit', N-1, N, CONE, AFAC( 1, 2 ), LDAFAC, C( 2, 1 ), LDC );
         } else {
            ctrmm('Left', UPLO, 'No transpose', 'Unit', N-1, N, CONE, AFAC( 2, 1 ), LDAFAC, C( 2, 1 ), LDC );
         }

         // Call CTRMM again to multiply by U (or L ).

         if ( lsame( UPLO, 'U' ) ) {
            ctrmm('Right', UPLO, 'No transpose', 'Unit', N, N-1, CONE, AFAC( 1, 2 ), LDAFAC, C( 1, 2 ), LDC );
         } else {
            ctrmm('Right', UPLO, 'Conjugate transpose', 'Unit', N, N-1, CONE, AFAC( 2, 1 ), LDAFAC, C( 1, 2 ), LDC );
         }
      }

      // Apply hermitian pivots

      for (J = N; J >= 1; J--) {
         I = IPIV( J );
         if (I != J) cswap( N, C( J, 1 ), LDC, C( I, 1 ), LDC );
      }
      for (J = N; J >= 1; J--) {
         I = IPIV( J );
         if (I != J) cswap( N, C( 1, J ), 1, C( 1, I ), 1 );
      }


      // Compute the difference  C - A .

      if ( lsame( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= J; I++) {
               C[I][J] = C( I, J ) - A( I, J );
            }
         }
      } else {
         for (J = 1; J <= N; J++) {
            for (I = J; I <= N; I++) {
               C[I][J] = C( I, J ) - A( I, J );
            }
         }
      }

      // Compute norm( C - A ) / ( N * norm(A) * EPS )

      RESID = CLANHE( '1', UPLO, N, C, LDC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS;
      }

      }

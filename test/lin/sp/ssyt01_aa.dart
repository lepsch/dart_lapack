      void ssyt01_aa(UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                LDA, LDAFAC, LDC, N;
      double               RESID;
      int                IPIV( * );
      double               A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ), RWORK( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, J;
      double               ANORM, EPS;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH, SLANSY;
      // EXTERNAL lsame, SLAMCH, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLASET, SLAVSY, SSWAP, STRMM, SLACPY
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
      ANORM = SLANSY( '1', UPLO, N, A, LDA, RWORK );

      // Initialize C to the tridiagonal matrix T.

      slaset('Full', N, N, ZERO, ZERO, C, LDC );
      slacpy('F', 1, N, AFAC( 1, 1 ), LDAFAC+1, C( 1, 1 ), LDC+1 );
      if ( N > 1 ) {
         if ( lsame( UPLO, 'U' ) ) {
            slacpy('F', 1, N-1, AFAC( 1, 2 ), LDAFAC+1, C( 1, 2 ), LDC+1 );
            slacpy('F', 1, N-1, AFAC( 1, 2 ), LDAFAC+1, C( 2, 1 ), LDC+1 );
         } else {
            slacpy('F', 1, N-1, AFAC( 2, 1 ), LDAFAC+1, C( 1, 2 ), LDC+1 );
            slacpy('F', 1, N-1, AFAC( 2, 1 ), LDAFAC+1, C( 2, 1 ), LDC+1 );
         }

         // Call STRMM to form the product U' * D (or L * D ).

         if ( lsame( UPLO, 'U' ) ) {
            strmm('Left', UPLO, 'Transpose', 'Unit', N-1, N, ONE, AFAC( 1, 2 ), LDAFAC, C( 2, 1 ), LDC );
         } else {
            strmm('Left', UPLO, 'No transpose', 'Unit', N-1, N, ONE, AFAC( 2, 1 ), LDAFAC, C( 2, 1 ), LDC );
         }

         // Call STRMM again to multiply by U (or L ).

         if ( lsame( UPLO, 'U' ) ) {
            strmm('Right', UPLO, 'No transpose', 'Unit', N, N-1, ONE, AFAC( 1, 2 ), LDAFAC, C( 1, 2 ), LDC );
         } else {
            strmm('Right', UPLO, 'Transpose', 'Unit', N, N-1, ONE, AFAC( 2, 1 ), LDAFAC, C( 1, 2 ), LDC );
         }
      }

      // Apply symmetric pivots

      for (J = N; J >= 1; J--) {
         I = IPIV( J );
         if (I != J) sswap( N, C( J, 1 ), LDC, C( I, 1 ), LDC );
      }
      for (J = N; J >= 1; J--) {
         I = IPIV( J );
         if (I != J) sswap( N, C( 1, J ), 1, C( 1, I ), 1 );
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

      RESID = SLANSY( '1', UPLO, N, C, LDC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS;
      }

      return;
      }

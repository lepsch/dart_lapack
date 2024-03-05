      void dppt01(final int UPLO, final int N, final int A, final int AFAC, final Array<double> RWORK_, final int RESID,) {
  final RWORK = RWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                N;
      double             RESID;
      double             A( * ), AFAC( * ), RWORK( * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, K, KC, NPP;
      double             ANORM, EPS, T;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DDOT, DLAMCH, DLANSP;
      // EXTERNAL lsame, DDOT, DLAMCH, DLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSPR, DTPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE

      // Quick exit if N = 0

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = dlamch( 'Epsilon' );
      ANORM = dlansp( '1', UPLO, N, A, RWORK );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute the product U'*U, overwriting U.

      if ( lsame( UPLO, 'U' ) ) {
         KC = ( N*( N-1 ) ) / 2 + 1;
         for (K = N; K >= 1; K--) { // 10

            // Compute the (K,K) element of the result.

            T = ddot( K, AFAC( KC ), 1, AFAC( KC ), 1 );
            AFAC[KC+K-1] = T;

            // Compute the rest of column K.

            if ( K > 1 ) {
               dtpmv('Upper', 'Transpose', 'Non-unit', K-1, AFAC, AFAC( KC ), 1 );
               KC = KC - ( K-1 );
            }
         } // 10

      // Compute the product L*L', overwriting L.

      } else {
         KC = ( N*( N+1 ) ) / 2;
         for (K = N; K >= 1; K--) { // 20

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            if (K < N) dspr( 'Lower', N-K, ONE, AFAC( KC+1 ), 1, AFAC( KC+N-K+1 ) );

            // Scale column K by the diagonal element.

            T = AFAC( KC );
            dscal(N-K+1, T, AFAC( KC ), 1 );

            KC = KC - ( N-K+2 );
         } // 20
      }

      // Compute the difference  L*L' - A (or U'*U - A).

      NPP = N*( N+1 ) / 2;
      for (I = 1; I <= NPP; I++) { // 30
         AFAC[I] = AFAC( I ) - A( I );
      } // 30

      // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

      RESID = dlansp( '1', UPLO, N, AFAC, RWORK );

      RESID = ( ( RESID / N.toDouble() ) / ANORM ) / EPS;

      }

      void cppt01(final int UPLO, final int N, final int A, final int AFAC, final Array<double> RWORK, final int RESID,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                N;
      double               RESID;
      double               RWORK( * );
      Complex            A( * ), AFAC( * );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, K, KC;
      double               ANORM, EPS, TR;
      Complex            TC;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANHP, SLAMCH;
      //- COMPLEX            CDOTC;
      // EXTERNAL lsame, CLANHP, SLAMCH, CDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHPR, CSCAL, CTPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, REAL

      // Quick exit if N = 0

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANHP( '1', UPLO, N, A, RWORK );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Check the imaginary parts of the diagonal elements and return with
      // an error code if any are nonzero.

      KC = 1;
      if ( lsame( UPLO, 'U' ) ) {
         for (K = 1; K <= N; K++) { // 10
            if ( AIMAG( AFAC( KC ) ) != ZERO ) {
               RESID = ONE / EPS;
               return;
            }
            KC = KC + K + 1;
         } // 10
      } else {
         for (K = 1; K <= N; K++) { // 20
            if ( AIMAG( AFAC( KC ) ) != ZERO ) {
               RESID = ONE / EPS;
               return;
            }
            KC = KC + N - K + 1;
         } // 20
      }

      // Compute the product U'*U, overwriting U.

      if ( lsame( UPLO, 'U' ) ) {
         KC = ( N*( N-1 ) ) / 2 + 1;
         for (K = N; K >= 1; K--) { // 30

            // Compute the (K,K) element of the result.

            TR = double( CDOTC( K, AFAC( KC ), 1, AFAC( KC ), 1 ) );
            AFAC[KC+K-1] = TR;

            // Compute the rest of column K.

            if ( K > 1 ) {
               ctpmv('Upper', 'Conjugate', 'Non-unit', K-1, AFAC, AFAC( KC ), 1 );
               KC = KC - ( K-1 );
            }
         } // 30

         // Compute the difference  L*L' - A

         KC = 1;
         for (K = 1; K <= N; K++) { // 50
            for (I = 1; I <= K - 1; I++) { // 40
               AFAC[KC+I-1] = AFAC( KC+I-1 ) - A( KC+I-1 );
            } // 40
            AFAC[KC+K-1] = AFAC( KC+K-1 ) - double( A( KC+K-1 ) );
            KC = KC + K;
         } // 50

      // Compute the product L*L', overwriting L.

      } else {
         KC = ( N*( N+1 ) ) / 2;
         for (K = N; K >= 1; K--) { // 60

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            if (K < N) chpr( 'Lower', N-K, ONE, AFAC( KC+1 ), 1, AFAC( KC+N-K+1 ) );

            // Scale column K by the diagonal element.

            TC = AFAC( KC );
            cscal(N-K+1, TC, AFAC( KC ), 1 );

            KC = KC - ( N-K+2 );
         } // 60

         // Compute the difference  U'*U - A

         KC = 1;
         for (K = 1; K <= N; K++) { // 80
            AFAC[KC] = AFAC( KC ) - double( A( KC ) );
            for (I = K + 1; I <= N; I++) { // 70
               AFAC[KC+I-K] = AFAC( KC+I-K ) - A( KC+I-K );
            } // 70
            KC = KC + N - K + 1;
         } // 80
      }

      // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

      RESID = CLANHP( '1', UPLO, N, AFAC, RWORK );

      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS;

      }

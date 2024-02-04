      void dpot01(UPLO, N, A, LDA, AFAC, LDAFAC, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, N;
      double             RESID;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      double             ANORM, EPS, T;
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- double             DDOT, DLAMCH, DLANSY;
      // EXTERNAL LSAME, DDOT, DLAMCH, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSCAL, DSYR, DTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = DLAMCH( 'Epsilon' );
      ANORM = DLANSY( '1', UPLO, N, A, LDA, RWORK );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute the product U**T * U, overwriting U.

      if ( LSAME( UPLO, 'U' ) ) {
         for (K = N; K >= 1; K--) { // 10

            // Compute the (K,K) element of the result.

            T = DDOT( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 );
            AFAC[K, K] = T;

            // Compute the rest of column K.

            dtrmv('Upper', 'Transpose', 'Non-unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );

         } // 10

      // Compute the product L * L**T, overwriting L.

      } else {
         for (K = N; K >= 1; K--) { // 20

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            if (K+1 <= N) dsyr( 'Lower', N-K, ONE, AFAC( K+1, K ), 1, AFAC( K+1, K+1 ), LDAFAC );

            // Scale column K by the diagonal element.

            T = AFAC( K, K );
            dscal(N-K+1, T, AFAC( K, K ), 1 );

         } // 20
      }

      // Compute the difference L * L**T - A (or U**T * U - A).

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 40
            for (I = 1; I <= J; I++) { // 30
               AFAC[I, J] = AFAC( I, J ) - A( I, J );
            } // 30
         } // 40
      } else {
         for (J = 1; J <= N; J++) { // 60
            for (I = J; I <= N; I++) { // 50
               AFAC[I, J] = AFAC( I, J ) - A( I, J );
            } // 50
         } // 60
      }

      // Compute norm(L*U - A) / ( N * norm(A) * EPS )

      RESID = DLANSY( '1', UPLO, N, AFAC, LDAFAC, RWORK );

      RESID = ( ( RESID / N.toDouble() ) / ANORM ) / EPS;

      return;
      }

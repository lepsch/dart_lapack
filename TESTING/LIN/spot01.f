      SUBROUTINE SPOT01( UPLO, N, A, LDA, AFAC, LDAFAC, RWORK, RESID );

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                LDA, LDAFAC, N;
      REAL               RESID;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      REAL               ANORM, EPS, T;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SDOT, SLAMCH, SLANSY;
      // EXTERNAL LSAME, SDOT, SLAMCH, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSYR, STRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0.

      if ( N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = SLANSY( '1', UPLO, N, A, LDA, RWORK );
      if ( ANORM <= ZERO ) {
         RESID = ONE / EPS;
         return;
      }

      // Compute the product U**T * U, overwriting U.

      if ( LSAME( UPLO, 'U' ) ) {
         DO 10 K = N, 1, -1;

            // Compute the (K,K) element of the result.

            T = SDOT( K, AFAC( 1, K ), 1, AFAC( 1, K ), 1 );
            AFAC( K, K ) = T;

            // Compute the rest of column K.

            strmv('Upper', 'Transpose', 'Non-unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );

         } // 10

      // Compute the product L * L**T, overwriting L.

      } else {
         DO 20 K = N, 1, -1;

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            if (K+1 <= N) CALL SSYR( 'Lower', N-K, ONE, AFAC( K+1, K ), 1, AFAC( K+1, K+1 ), LDAFAC );

            // Scale column K by the diagonal element.

            T = AFAC( K, K );
            sscal(N-K+1, T, AFAC( K, K ), 1 );

         } // 20
      }

      // Compute the difference L * L**T - A (or U**T * U - A).

      if ( LSAME( UPLO, 'U' ) ) {
         for (J = 1; J <= N; J++) { // 40
            for (I = 1; I <= J; I++) { // 30
               AFAC( I, J ) = AFAC( I, J ) - A( I, J );
            } // 30
         } // 40
      } else {
         for (J = 1; J <= N; J++) { // 60
            for (I = J; I <= N; I++) { // 50
               AFAC( I, J ) = AFAC( I, J ) - A( I, J );
            } // 50
         } // 60
      }

      // Compute norm(L*U - A) / ( N * norm(A) * EPS )

      RESID = SLANSY( '1', UPLO, N, AFAC, LDAFAC, RWORK );

      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS;

      return;

      // End of SPOT01

      }

      SUBROUTINE SPPT01( UPLO, N, A, AFAC, RWORK, RESID )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                N;
      REAL               RESID
      // ..
      // .. Array Arguments ..
      REAL               A( * ), AFAC( * ), RWORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, K, KC, NPP;
      REAL               ANORM, EPS, T
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SDOT, SLAMCH, SLANSP
      // EXTERNAL LSAME, SDOT, SLAMCH, SLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL SSCAL, SSPR, STPMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if N = 0

      if ( N.LE.0 ) {
         RESID = ZERO
         RETURN
      }

      // Exit with RESID = 1/EPS if ANORM = 0.

      EPS = SLAMCH( 'Epsilon' )
      ANORM = SLANSP( '1', UPLO, N, A, RWORK )
      if ( ANORM.LE.ZERO ) {
         RESID = ONE / EPS
         RETURN
      }

      // Compute the product U'*U, overwriting U.

      if ( LSAME( UPLO, 'U' ) ) {
         KC = ( N*( N-1 ) ) / 2 + 1
         DO 10 K = N, 1, -1

            // Compute the (K,K) element of the result.

            T = SDOT( K, AFAC( KC ), 1, AFAC( KC ), 1 )
            AFAC( KC+K-1 ) = T

            // Compute the rest of column K.

            if ( K.GT.1 ) {
               stpmv('Upper', 'Transpose', 'Non-unit', K-1, AFAC, AFAC( KC ), 1 );
               KC = KC - ( K-1 )
            }
         } // 10

      // Compute the product L*L', overwriting L.

      } else {
         KC = ( N*( N+1 ) ) / 2
         DO 20 K = N, 1, -1

            // Add a multiple of column K of the factor L to each of
            // columns K+1 through N.

            if (K < N) CALL SSPR( 'Lower', N-K, ONE, AFAC( KC+1 ), 1, AFAC( KC+N-K+1 ) );

            // Scale column K by the diagonal element.

            T = AFAC( KC )
            sscal(N-K+1, T, AFAC( KC ), 1 );

            KC = KC - ( N-K+2 )
         } // 20
      }

      // Compute the difference  L*L' - A (or U'*U - A).

      NPP = N*( N+1 ) / 2
      for (I = 1; I <= NPP; I++) { // 30
         AFAC( I ) = AFAC( I ) - A( I )
      } // 30

      // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

      RESID = SLANSP( '1', UPLO, N, AFAC, RWORK )

      RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS

      RETURN

      // End of SPPT01

      }

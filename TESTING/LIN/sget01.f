      void sget01(M, N, A, LDA, AFAC, LDAFAC, IPIV, RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDAFAC, M, N;
      REAL               RESID;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
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
      //- REAL               SDOT, SLAMCH, SLANGE;
      // EXTERNAL SDOT, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, SLASWP, SSCAL, STRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, REAL
      // ..
      // .. Executable Statements ..

      // Quick exit if M = 0 or N = 0.

      if ( M <= 0 || N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Determine EPS and the norm of A.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = SLANGE( '1', M, N, A, LDA, RWORK );

      // Compute the product L*U and overwrite AFAC with the result.
      // A column at a time of the product is obtained, starting with
      // column N.

      DO 10 K = N, 1, -1;
         if ( K > M ) {
            strmv('Lower', 'No transpose', 'Unit', M, AFAC, LDAFAC, AFAC( 1, K ), 1 );
         } else {

            // Compute elements (K+1:M,K)

            T = AFAC( K, K );
            if ( K+1 <= M ) {
               sscal(M-K, T, AFAC( K+1, K ), 1 );
               sgemv('No transpose', M-K, K-1, ONE, AFAC( K+1, 1 ), LDAFAC, AFAC( 1, K ), 1, ONE, AFAC( K+1, K ), 1 );
            }

            // Compute the (K,K) element

            AFAC( K, K ) = T + SDOT( K-1, AFAC( K, 1 ), LDAFAC, AFAC( 1, K ), 1 );

            // Compute elements (1:K-1,K)

            strmv('Lower', 'No transpose', 'Unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );
         }
      } // 10
      slaswp(N, AFAC, LDAFAC, 1, min( M, N ), IPIV, -1 );

      // Compute the difference  L*U - A  and store in AFAC.

      for (J = 1; J <= N; J++) { // 30
         for (I = 1; I <= M; I++) { // 20
            AFAC( I, J ) = AFAC( I, J ) - A( I, J );
         } // 20
      } // 30

      // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

      RESID = SLANGE( '1', M, N, AFAC, LDAFAC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / REAL( N ) ) / ANORM ) / EPS;
      }

      return;
      }

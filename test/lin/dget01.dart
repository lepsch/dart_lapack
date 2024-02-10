      void dget01(M, N, final Matrix<double> A, final int LDA, final Matrix<double> AFAC, final int LDAFAC, final Array<int> IPIV, final Array<double> RWORK, final int RESID) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDAFAC, M, N;
      double             RESID;
      int                IPIV( * );
      double             A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * );
      // ..

// =====================================================================


      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, J, K;
      double             ANORM, EPS, T;
      // ..
      // .. External Functions ..
      //- double             DDOT, DLAMCH, DLANGE;
      // EXTERNAL DDOT, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMV, DLASWP, DSCAL, DTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MIN

      // Quick exit if M = 0 or N = 0.

      if ( M <= 0 || N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Determine EPS and the norm of A.

      EPS = dlamch( 'Epsilon' );
      ANORM = dlange( '1', M, N, A, LDA, RWORK );

      // Compute the product L*U and overwrite AFAC with the result.
      // A column at a time of the product is obtained, starting with
      // column N.

      for (K = N; K >= 1; K--) { // 10
         if ( K > M ) {
            dtrmv('Lower', 'No transpose', 'Unit', M, AFAC, LDAFAC, AFAC( 1, K ), 1 );
         } else {

            // Compute elements (K+1:M,K)

            T = AFAC( K, K );
            if ( K+1 <= M ) {
               dscal(M-K, T, AFAC( K+1, K ), 1 );
               dgemv('No transpose', M-K, K-1, ONE, AFAC( K+1, 1 ), LDAFAC, AFAC( 1, K ), 1, ONE, AFAC( K+1, K ), 1 );
            }

            // Compute the (K,K) element

            AFAC[K][K] = T + ddot( K-1, AFAC( K, 1 ), LDAFAC, AFAC( 1, K ), 1 );

            // Compute elements (1:K-1,K)

            dtrmv('Lower', 'No transpose', 'Unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );
         }
      } // 10
      dlaswp(N, AFAC, LDAFAC, 1, min( M, N ), IPIV, -1 );

      // Compute the difference  L*U - A  and store in AFAC.

      for (J = 1; J <= N; J++) { // 30
         for (I = 1; I <= M; I++) { // 20
            AFAC[I][J] = AFAC( I, J ) - A( I, J );
         } // 20
      } // 30

      // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

      RESID = dlange( '1', M, N, AFAC, LDAFAC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / N.toDouble() ) / ANORM ) / EPS;
      }

      }

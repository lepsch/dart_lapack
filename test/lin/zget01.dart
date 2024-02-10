      void zget01(M, N, final Matrix<double> A, final int LDA, final Matrix<double> AFAC, final int LDAFAC, final Array<int> IPIV, final Array<double> RWORK, RESID ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDAFAC, M, N;
      double             RESID;
      int                IPIV( * );
      double             RWORK( * );
      Complex         A( LDA, * ), AFAC( LDAFAC, * );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      int                I, J, K;
      double             ANORM, EPS;
      Complex         T;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, ZLANGE;
      //- Complex         ZDOTU;
      // EXTERNAL DLAMCH, ZLANGE, ZDOTU
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMV, ZLASWP, ZSCAL, ZTRMV
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
      ANORM = ZLANGE( '1', M, N, A, LDA, RWORK );

      // Compute the product L*U and overwrite AFAC with the result.
      // A column at a time of the product is obtained, starting with
      // column N.

      for (K = N; K >= 1; K--) { // 10
         if ( K > M ) {
            ztrmv('Lower', 'No transpose', 'Unit', M, AFAC, LDAFAC, AFAC( 1, K ), 1 );
         } else {

            // Compute elements (K+1:M,K)

            T = AFAC( K, K );
            if ( K+1 <= M ) {
               zscal(M-K, T, AFAC( K+1, K ), 1 );
               zgemv('No transpose', M-K, K-1, CONE, AFAC( K+1, 1 ), LDAFAC, AFAC( 1, K ), 1, CONE, AFAC( K+1, K ), 1 );
            }

            // Compute the (K,K) element

            AFAC[K][K] = T + ZDOTU( K-1, AFAC( K, 1 ), LDAFAC, AFAC( 1, K ), 1 );

            // Compute elements (1:K-1,K)

            ztrmv('Lower', 'No transpose', 'Unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );
         }
      } // 10
      zlaswp(N, AFAC, LDAFAC, 1, min( M, N ), IPIV, -1 );

      // Compute the difference  L*U - A  and store in AFAC.

      for (J = 1; J <= N; J++) { // 30
         for (I = 1; I <= M; I++) { // 20
            AFAC[I][J] = AFAC( I, J ) - A( I, J );
         } // 20
      } // 30

      // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

      RESID = ZLANGE( '1', M, N, AFAC, LDAFAC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID / N.toDouble() ) / ANORM ) / EPS;
      }

      }

      void cget01(final int M, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> AFAC_, final int LDAFAC, final Array<int> IPIV_, final Array<double> RWORK_, final int RESID,) {
  final A = A_.dim();
  final AFAC = AFAC_.dim();
  final IPIV = IPIV_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDAFAC, M, N;
      double               RESID;
      int                IPIV( * );
      double               RWORK( * );
      Complex            A( LDA, * ), AFAC( LDAFAC, * );
      // ..

      double               ONE, ZERO;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      int                I, J, K;
      double               ANORM, EPS;
      Complex            T;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, SLAMCH;
      //- COMPLEX            CDOTU;
      // EXTERNAL CLANGE, SLAMCH, CDOTU
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CLASWP, CSCAL, CTRMV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, REAL

      // Quick exit if M = 0 or N = 0.

      if ( M <= 0 || N <= 0 ) {
         RESID = ZERO;
         return;
      }

      // Determine EPS and the norm of A.

      EPS = SLAMCH( 'Epsilon' );
      ANORM = CLANGE( '1', M, N, A, LDA, RWORK );

      // Compute the product L*U and overwrite AFAC with the result.
      // A column at a time of the product is obtained, starting with
      // column N.

      for (K = N; K >= 1; K--) { // 10
         if ( K > M ) {
            ctrmv('Lower', 'No transpose', 'Unit', M, AFAC, LDAFAC, AFAC( 1, K ), 1 );
         } else {

            // Compute elements (K+1:M,K)

            T = AFAC( K, K );
            if ( K+1 <= M ) {
               cscal(M-K, T, AFAC( K+1, K ), 1 );
               cgemv('No transpose', M-K, K-1, CONE, AFAC( K+1, 1 ), LDAFAC, AFAC( 1, K ), 1, CONE, AFAC( K+1, K ), 1 );
            }

            // Compute the (K,K) element

            AFAC[K][K] = T + CDOTU( K-1, AFAC( K, 1 ), LDAFAC, AFAC( 1, K ), 1 );

            // Compute elements (1:K-1,K)

            ctrmv('Lower', 'No transpose', 'Unit', K-1, AFAC, LDAFAC, AFAC( 1, K ), 1 );
         }
      } // 10
      claswp(N, AFAC, LDAFAC, 1, min( M, N ), IPIV, -1 );

      // Compute the difference  L*U - A  and store in AFAC.

      for (J = 1; J <= N; J++) { // 30
         for (I = 1; I <= M; I++) { // 20
            AFAC[I][J] = AFAC( I, J ) - A( I, J );
         } // 20
      } // 30

      // Compute norm( L*U - A ) / ( N * norm(A) * EPS )

      RESID = CLANGE( '1', M, N, AFAC, LDAFAC, RWORK );

      if ( ANORM <= ZERO ) {
         if (RESID != ZERO) RESID = ONE / EPS;
      } else {
         RESID = ( ( RESID/REAL( N ) )/ANORM ) / EPS;
      }

      }

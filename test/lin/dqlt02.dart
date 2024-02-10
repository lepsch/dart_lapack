      void dqlt02(final int M, final int N, final int K, final int A, final int AF, final int Q, final int L, final int LDA, final int TAU, final Array<double> WORK, final int LWORK, final Array<double> RWORK, final int RESULT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M, N;
      double             A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             ROGUE;
      const              ROGUE = -1.0e+10 ;
      int                INFO;
      double             ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLANGE, DLANSY;
      // EXTERNAL DLAMCH, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DLASET, DORGQL, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT

      // Quick return if possible

      if ( M == 0 || N == 0 || K == 0 ) {
         RESULT[1] = ZERO;
         RESULT[2] = ZERO;
         return;
      }

      EPS = dlamch( 'Epsilon' );

      // Copy the last k columns of the factorization to the array Q

      dlaset('Full', M, N, ROGUE, ROGUE, Q, LDA );
      if (K < M) dlacpy( 'Full', M-K, K, AF( 1, N-K+1 ), LDA, Q( 1, N-K+1 ), LDA );
      IF( K > 1 ) dlacpy( 'Upper', K-1, K-1, AF( M-K+1, N-K+2 ), LDA, Q( M-K+1, N-K+2 ), LDA );

      // Generate the last n columns of the matrix Q

     srnamc.SRNAMT = 'DORGQL';
      dorgql(M, N, K, Q, LDA, TAU( N-K+1 ), WORK, LWORK, INFO );

      // Copy L(m-n+1:m,n-k+1:n)

      dlaset('Full', N, K, ZERO, ZERO, L( M-N+1, N-K+1 ), LDA );
      dlacpy('Lower', K, K, AF( M-K+1, N-K+1 ), LDA, L( M-K+1, N-K+1 ), LDA );

      // Compute L(m-n+1:m,n-k+1:n) - Q(1:m,m-n+1:m)' * A(1:m,n-k+1:n)

      dgemm('Transpose', 'No transpose', N, K, M, -ONE, Q, LDA, A( 1, N-K+1 ), LDA, ONE, L( M-N+1, N-K+1 ), LDA );

      // Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = dlange( '1', M, K, A( 1, N-K+1 ), LDA, RWORK );
      RESID = dlange( '1', N, K, L( M-N+1, N-K+1 ), LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / (max( 1, M )).toDouble() ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q'*Q

      dlaset('Full', N, N, ZERO, ONE, L, LDA );
      dsyrk('Upper', 'Transpose', N, M, -ONE, Q, LDA, ONE, L, LDA );

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = dlansy( '1', 'Upper', N, L, LDA, RWORK );

      RESULT[2] = ( RESID / (max( 1, M )).toDouble() ) / EPS;

      }

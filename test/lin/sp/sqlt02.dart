      void sqlt02(M, N, K, A, AF, Q, L, LDA, TAU, final Array<double> WORK, final int LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M, N;
      double               A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               ROGUE;
      const              ROGUE = -1.0e+10 ;
      int                INFO;
      double               ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGE, SLANSY;
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY, SLASET, SORGQL, SSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL
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

      EPS = SLAMCH( 'Epsilon' );

      // Copy the last k columns of the factorization to the array Q

      slaset('Full', M, N, ROGUE, ROGUE, Q, LDA );
      if (K < M) slacpy( 'Full', M-K, K, AF( 1, N-K+1 ), LDA, Q( 1, N-K+1 ), LDA );
      IF( K > 1 ) slacpy( 'Upper', K-1, K-1, AF( M-K+1, N-K+2 ), LDA, Q( M-K+1, N-K+2 ), LDA );

      // Generate the last n columns of the matrix Q

     srnamc.SRNAMT = 'SORGQL';
      sorgql(M, N, K, Q, LDA, TAU( N-K+1 ), WORK, LWORK, INFO );

      // Copy L(m-n+1:m,n-k+1:n)

      slaset('Full', N, K, ZERO, ZERO, L( M-N+1, N-K+1 ), LDA );
      slacpy('Lower', K, K, AF( M-K+1, N-K+1 ), LDA, L( M-K+1, N-K+1 ), LDA );

      // Compute L(m-n+1:m,n-k+1:n) - Q(1:m,m-n+1:m)' * A(1:m,n-k+1:n)

      sgemm('Transpose', 'No transpose', N, K, M, -ONE, Q, LDA, A( 1, N-K+1 ), LDA, ONE, L( M-N+1, N-K+1 ), LDA );

      // Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = SLANGE( '1', M, K, A( 1, N-K+1 ), LDA, RWORK );
      RESID = SLANGE( '1', N, K, L( M-N+1, N-K+1 ), LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / REAL( max( 1, M ) ) ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q'*Q

      slaset('Full', N, N, ZERO, ONE, L, LDA );
      ssyrk('Upper', 'Transpose', N, M, -ONE, Q, LDA, ONE, L, LDA );

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = SLANSY( '1', 'Upper', N, L, LDA, RWORK );

      RESULT[2] = ( RESID / REAL( max( 1, M ) ) ) / EPS;

      }

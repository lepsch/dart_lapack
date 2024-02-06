      void zqlt02(M, N, K, A, AF, Q, L, LDA, TAU, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M, N;
      double             RESULT( * ), RWORK( * );
      Complex         A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         ROGUE;
      const              ROGUE = ( -1.0e+10, -1.0e+10 ) ;
      int                INFO;
      double             ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, ZLANGE, ZLANSY;
      // EXTERNAL DLAMCH, ZLANGE, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZHERK, ZLACPY, ZLASET, ZUNGQL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX
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

      zlaset('Full', M, N, ROGUE, ROGUE, Q, LDA );
      if (K < M) zlacpy( 'Full', M-K, K, AF( 1, N-K+1 ), LDA, Q( 1, N-K+1 ), LDA );
      IF( K > 1 ) zlacpy( 'Upper', K-1, K-1, AF( M-K+1, N-K+2 ), LDA, Q( M-K+1, N-K+2 ), LDA );

      // Generate the last n columns of the matrix Q

     srnamc.SRNAMT = 'ZUNGQL';
      zungql(M, N, K, Q, LDA, TAU( N-K+1 ), WORK, LWORK, INFO );

      // Copy L(m-n+1:m,n-k+1:n)

      zlaset('Full', N, K, DCMPLX( ZERO ), DCMPLX( ZERO ), L( M-N+1, N-K+1 ), LDA )       CALL ZLACPY( 'Lower', K, K, AF( M-K+1, N-K+1 ), LDA, L( M-K+1, N-K+1 ), LDA );

      // Compute L(m-n+1:m,n-k+1:n) - Q(1:m,m-n+1:m)' * A(1:m,n-k+1:n)

      zgemm('Conjugate transpose', 'No transpose', N, K, M, DCMPLX( -ONE ), Q, LDA, A( 1, N-K+1 ), LDA, DCMPLX( ONE ), L( M-N+1, N-K+1 ), LDA );

      // Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = ZLANGE( '1', M, K, A( 1, N-K+1 ), LDA, RWORK );
      RESID = ZLANGE( '1', N, K, L( M-N+1, N-K+1 ), LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / (max( 1, M )).toDouble() ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q'*Q

      zlaset('Full', N, N, DCMPLX( ZERO ), DCMPLX( ONE ), L, LDA );
      zherk('Upper', 'Conjugate transpose', N, M, -ONE, Q, LDA, ONE, L, LDA );

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = ZLANSY( '1', 'Upper', N, L, LDA, RWORK );

      RESULT[2] = ( RESID / (max( 1, M )).toDouble() ) / EPS;

      return;
      }

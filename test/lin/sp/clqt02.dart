      void clqt02(M, N, K, A, AF, Q, L, LDA, TAU, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M, N;
      double               RESULT( * ), RWORK( * );
      Complex            A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            ROGUE;
      const              ROGUE = ( -1.0e+10, -1.0e+10 ) ;
      int                INFO;
      double               ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, CLANSY, SLAMCH;
      // EXTERNAL CLANGE, CLANSY, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHERK, CLACPY, CLASET, CUNGLQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, REAL
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT

      EPS = SLAMCH( 'Epsilon' );

      // Copy the first k rows of the factorization to the array Q

      claset('Full', M, N, ROGUE, ROGUE, Q, LDA );
      clacpy('Upper', K, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA );

      // Generate the first n columns of the matrix Q

     srnamc.SRNAMT = 'CUNGLQ';
      cunglq(M, N, K, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy L(1:k,1:m)

      claset('Full', K, M, CMPLX( ZERO ), CMPLX( ZERO ), L, LDA );
      clacpy('Lower', K, M, AF, LDA, L, LDA );

      // Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)'

      cgemm('No transpose', 'Conjugate transpose', K, M, N, CMPLX( -ONE ), A, LDA, Q, LDA, CMPLX( ONE ), L, LDA );

      // Compute norm( L - A*Q' ) / ( N * norm(A) * EPS ) .

      ANORM = CLANGE( '1', K, N, A, LDA, RWORK );
      RESID = CLANGE( '1', K, M, L, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / REAL( max( 1, N ) ) ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q*Q'

      claset('Full', M, M, CMPLX( ZERO ), CMPLX( ONE ), L, LDA );
      cherk('Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, L, LDA );

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = CLANSY( '1', 'Upper', M, L, LDA, RWORK );

      RESULT[2] = ( RESID / REAL( max( 1, N ) ) ) / EPS;

      return;
      }

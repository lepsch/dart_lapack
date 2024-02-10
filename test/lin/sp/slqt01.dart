      void slqt01(M, N, A, AF, Q, L, LDA, TAU, final Array<double> WORK, final int LWORK, final Array<double> RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LWORK, M, N;
      double               A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               ROGUE;
      const              ROGUE = -1.0e+10 ;
      int                INFO, MINMN;
      double               ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGE, SLANSY;
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGELQF, SGEMM, SLACPY, SLASET, SORGLQ, SSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT

      MINMN = min( M, N );
      EPS = SLAMCH( 'Epsilon' );

      // Copy the matrix A to the array AF.

      slacpy('Full', M, N, A, LDA, AF, LDA );

      // Factorize the matrix A in the array AF.

     srnamc.SRNAMT = 'SGELQF';
      sgelqf(M, N, AF, LDA, TAU, WORK, LWORK, INFO );

      // Copy details of Q

      slaset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      if (N > 1) slacpy( 'Upper', M, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA );

      // Generate the n-by-n matrix Q

     srnamc.SRNAMT = 'SORGLQ';
      sorglq(N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy L

      slaset('Full', M, N, ZERO, ZERO, L, LDA );
      slacpy('Lower', M, N, AF, LDA, L, LDA );

      // Compute L - A*Q'

      sgemm('No transpose', 'Transpose', M, N, N, -ONE, A, LDA, Q, LDA, ONE, L, LDA );

      // Compute norm( L - Q'*A ) / ( N * norm(A) * EPS ) .

      ANORM = SLANGE( '1', M, N, A, LDA, RWORK );
      RESID = SLANGE( '1', M, N, L, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / REAL( max( 1, N ) ) ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q*Q'

      slaset('Full', N, N, ZERO, ONE, L, LDA );
      ssyrk('Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, L, LDA );

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = SLANSY( '1', 'Upper', N, L, LDA, RWORK );

      RESULT[2] = ( RESID / REAL( max( 1, N ) ) ) / EPS;

      }

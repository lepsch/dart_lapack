      void dlqt01(final int M, final int N, final int A, final int AF, final int Q, final int L, final int LDA, final int TAU, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final int RESULT,) {
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LWORK, M, N;
      double             A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             ROGUE;
      const              ROGUE = -1.0e+10 ;
      int                INFO, MINMN;
      double             ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLANGE, DLANSY;
      // EXTERNAL DLAMCH, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGELQF, DGEMM, DLACPY, DLASET, DORGLQ, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT

      MINMN = min( M, N );
      EPS = dlamch( 'Epsilon' );

      // Copy the matrix A to the array AF.

      dlacpy('Full', M, N, A, LDA, AF, LDA );

      // Factorize the matrix A in the array AF.

     srnamc.SRNAMT = 'DGELQF';
      dgelqf(M, N, AF, LDA, TAU, WORK, LWORK, INFO );

      // Copy details of Q

      dlaset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      if (N > 1) dlacpy( 'Upper', M, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA );

      // Generate the n-by-n matrix Q

     srnamc.SRNAMT = 'DORGLQ';
      dorglq(N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy L

      dlaset('Full', M, N, ZERO, ZERO, L, LDA );
      dlacpy('Lower', M, N, AF, LDA, L, LDA );

      // Compute L - A*Q'

      dgemm('No transpose', 'Transpose', M, N, N, -ONE, A, LDA, Q, LDA, ONE, L, LDA );

      // Compute norm( L - Q'*A ) / ( N * norm(A) * EPS ) .

      ANORM = dlange( '1', M, N, A, LDA, RWORK );
      RESID = dlange( '1', M, N, L, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / (max( 1, N )).toDouble() ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q*Q'

      dlaset('Full', N, N, ZERO, ONE, L, LDA );
      dsyrk('Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, L, LDA );

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = dlansy( '1', 'Upper', N, L, LDA, RWORK );

      RESULT[2] = ( RESID / (max( 1, N )).toDouble() ) / EPS;

      }

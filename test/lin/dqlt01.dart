      void dqlt01(final int M, final int N, final int A, final int AF, final int Q, final int L, final int LDA, final int TAU, final Array<double> WORK, final int LWORK, final Array<double> RWORK, final int RESULT) {

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
      // EXTERNAL DGEMM, DGEQLF, DLACPY, DLASET, DORGQL, DSYRK
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

     srnamc.SRNAMT = 'DGEQLF';
      dgeqlf(M, N, AF, LDA, TAU, WORK, LWORK, INFO );

      // Copy details of Q

      dlaset('Full', M, M, ROGUE, ROGUE, Q, LDA );
      if ( M >= N ) {
         if (N < M && N > 0) dlacpy( 'Full', M-N, N, AF, LDA, Q( 1, M-N+1 ), LDA );
         IF( N > 1 ) dlacpy( 'Upper', N-1, N-1, AF( M-N+1, 2 ), LDA, Q( M-N+1, M-N+2 ), LDA );
      } else {
         if (M > 1) dlacpy( 'Upper', M-1, M-1, AF( 1, N-M+2 ), LDA, Q( 1, 2 ), LDA );
      }

      // Generate the m-by-m matrix Q

     srnamc.SRNAMT = 'DORGQL';
      dorgql(M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy L

      dlaset('Full', M, N, ZERO, ZERO, L, LDA );
      if ( M >= N ) {
         if (N > 0) dlacpy( 'Lower', N, N, AF( M-N+1, 1 ), LDA, L( M-N+1, 1 ), LDA );
      } else {
         if (N > M && M > 0) dlacpy( 'Full', M, N-M, AF, LDA, L, LDA );
         IF( M > 0 ) dlacpy( 'Lower', M, M, AF( 1, N-M+1 ), LDA, L( 1, N-M+1 ), LDA );
      }

      // Compute L - Q'*A

      dgemm('Transpose', 'No transpose', M, N, M, -ONE, Q, LDA, A, LDA, ONE, L, LDA );

      // Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = dlange( '1', M, N, A, LDA, RWORK );
      RESID = dlange( '1', M, N, L, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / (max( 1, M )).toDouble() ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q'*Q

      dlaset('Full', M, M, ZERO, ONE, L, LDA );
      dsyrk('Upper', 'Transpose', M, M, -ONE, Q, LDA, ONE, L, LDA );

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = dlansy( '1', 'Upper', M, L, LDA, RWORK );

      RESULT[2] = ( RESID / (max( 1, M )).toDouble() ) / EPS;

      }

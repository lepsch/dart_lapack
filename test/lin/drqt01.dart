      void drqt01(M, N, A, AF, Q, R, LDA, TAU, final Array<double> WORK, final int LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LWORK, M, N;
      double             A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
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
      // EXTERNAL DGEMM, DGERQF, DLACPY, DLASET, DORGRQ, DSYRK
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

     srnamc.SRNAMT = 'DGERQF';
      dgerqf(M, N, AF, LDA, TAU, WORK, LWORK, INFO );

      // Copy details of Q

      dlaset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      if ( M <= N ) {
         if (M > 0 && M < N) dlacpy( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA );
         IF( M > 1 ) dlacpy( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA, Q( N-M+2, N-M+1 ), LDA );
      } else {
         if (N > 1) dlacpy( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA, Q( 2, 1 ), LDA );
      }

      // Generate the n-by-n matrix Q

     srnamc.SRNAMT = 'DORGRQ';
      dorgrq(N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy R

      dlaset('Full', M, N, ZERO, ZERO, R, LDA );
      if ( M <= N ) {
         if (M > 0) dlacpy( 'Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ), LDA );
      } else {
         if (M > N && N > 0) dlacpy( 'Full', M-N, N, AF, LDA, R, LDA );
         IF( N > 0 ) dlacpy( 'Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ), LDA );
      }

      // Compute R - A*Q'

      dgemm('No transpose', 'Transpose', M, N, N, -ONE, A, LDA, Q, LDA, ONE, R, LDA );

      // Compute norm( R - Q'*A ) / ( N * norm(A) * EPS ) .

      ANORM = dlange( '1', M, N, A, LDA, RWORK );
      RESID = dlange( '1', M, N, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / (max( 1, N )).toDouble() ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q*Q'

      dlaset('Full', N, N, ZERO, ONE, R, LDA );
      dsyrk('Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = dlansy( '1', 'Upper', N, R, LDA, RWORK );

      RESULT[2] = ( RESID / (max( 1, N )).toDouble() ) / EPS;

      }

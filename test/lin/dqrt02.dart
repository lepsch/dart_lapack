      void dqrt02(final int M, final int N, final int K, final int A, final int AF, final int Q, final int R, final int LDA, final int TAU, final Array<double> WORK, final int LWORK, final Array<double> RWORK, final int RESULT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M, N;
      double             A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
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
      // EXTERNAL DGEMM, DLACPY, DLASET, DORGQR, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT

      EPS = dlamch( 'Epsilon' );

      // Copy the first k columns of the factorization to the array Q

      dlaset('Full', M, N, ROGUE, ROGUE, Q, LDA );
      dlacpy('Lower', M-1, K, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA );

      // Generate the first n columns of the matrix Q

     srnamc.SRNAMT = 'DORGQR';
      dorgqr(M, N, K, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy R(1:n,1:k)

      dlaset('Full', N, K, ZERO, ZERO, R, LDA );
      dlacpy('Upper', N, K, AF, LDA, R, LDA );

      // Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k)

      dgemm('Transpose', 'No transpose', N, K, M, -ONE, Q, LDA, A, LDA, ONE, R, LDA );

      // Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = dlange( '1', M, K, A, LDA, RWORK );
      RESID = dlange( '1', N, K, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / (max( 1, M )).toDouble() ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q'*Q

      dlaset('Full', N, N, ZERO, ONE, R, LDA );
      dsyrk('Upper', 'Transpose', N, M, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = dlansy( '1', 'Upper', N, R, LDA, RWORK );

      RESULT[2] = ( RESID / (max( 1, M )).toDouble() ) / EPS;

      }

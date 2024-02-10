      void sqrt02(M, N, K, A, AF, Q, R, LDA, TAU, final Array<double> WORK, final int LWORK, final Array<double> RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M, N;
      double               A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
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
      // EXTERNAL SGEMM, SLACPY, SLASET, SORGQR, SSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT

      EPS = SLAMCH( 'Epsilon' );

      // Copy the first k columns of the factorization to the array Q

      slaset('Full', M, N, ROGUE, ROGUE, Q, LDA );
      slacpy('Lower', M-1, K, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA );

      // Generate the first n columns of the matrix Q

     srnamc.SRNAMT = 'SORGQR';
      sorgqr(M, N, K, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy R(1:n,1:k)

      slaset('Full', N, K, ZERO, ZERO, R, LDA );
      slacpy('Upper', N, K, AF, LDA, R, LDA );

      // Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k)

      sgemm('Transpose', 'No transpose', N, K, M, -ONE, Q, LDA, A, LDA, ONE, R, LDA );

      // Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = SLANGE( '1', M, K, A, LDA, RWORK );
      RESID = SLANGE( '1', N, K, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / REAL( max( 1, M ) ) ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q'*Q

      slaset('Full', N, N, ZERO, ONE, R, LDA );
      ssyrk('Upper', 'Transpose', N, M, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = SLANSY( '1', 'Upper', N, R, LDA, RWORK );

      RESULT[2] = ( RESID / REAL( max( 1, M ) ) ) / EPS;

      }

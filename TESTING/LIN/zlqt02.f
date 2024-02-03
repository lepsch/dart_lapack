      SUBROUTINE ZLQT02( M, N, K, A, AF, Q, L, LDA, TAU, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             RESULT( * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), TAU( * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         ROGUE
      const              ROGUE = ( -1.0D+10, -1.0D+10 ) ;
      // ..
      // .. Local Scalars ..
      int                INFO;
      double             ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE, ZLANSY;
      // EXTERNAL DLAMCH, ZLANGE, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZHERK, ZLACPY, ZLASET, ZUNGLQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      EPS = DLAMCH( 'Epsilon' )

      // Copy the first k rows of the factorization to the array Q

      zlaset('Full', M, N, ROGUE, ROGUE, Q, LDA );
      zlacpy('Upper', K, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA );

      // Generate the first n columns of the matrix Q

      SRNAMT = 'ZUNGLQ'
      zunglq(M, N, K, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy L(1:k,1:m)

      zlaset('Full', K, M, DCMPLX( ZERO ), DCMPLX( ZERO ), L, LDA );
      zlacpy('Lower', K, M, AF, LDA, L, LDA );

      // Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)'

      zgemm('No transpose', 'Conjugate transpose', K, M, N, DCMPLX( -ONE ), A, LDA, Q, LDA, DCMPLX( ONE ), L, LDA );

      // Compute norm( L - A*Q' ) / ( N * norm(A) * EPS ) .

      ANORM = ZLANGE( '1', K, N, A, LDA, RWORK )
      RESID = ZLANGE( '1', K, M, L, LDA, RWORK )
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, N ) ) ) / ANORM ) / EPS
      } else {
         RESULT( 1 ) = ZERO
      }

      // Compute I - Q*Q'

      zlaset('Full', M, M, DCMPLX( ZERO ), DCMPLX( ONE ), L, LDA );
      zherk('Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, L, LDA );

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = ZLANSY( '1', 'Upper', M, L, LDA, RWORK )

      RESULT( 2 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / EPS

      RETURN

      // End of ZLQT02

      }

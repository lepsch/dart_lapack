      void drqt02(M, N, K, A, AF, Q, R, LDA, TAU, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double             ROGUE;
      const              ROGUE = -1.0e+10 ;
      // ..
      // .. Local Scalars ..
      int                INFO;
      double             ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, DLANGE, DLANSY;
      // EXTERNAL DLAMCH, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DLASET, DORGRQ, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( M == 0 || N == 0 || K == 0 ) {
         RESULT( 1 ) = ZERO;
         RESULT( 2 ) = ZERO;
         return;
      }

      EPS = DLAMCH( 'Epsilon' );

      // Copy the last k rows of the factorization to the array Q

      dlaset('Full', M, N, ROGUE, ROGUE, Q, LDA );
      if (K < N) dlacpy( 'Full', K, N-K, AF( M-K+1, 1 ), LDA, Q( M-K+1, 1 ), LDA );
      IF( K > 1 ) dlacpy( 'Lower', K-1, K-1, AF( M-K+2, N-K+1 ), LDA, Q( M-K+2, N-K+1 ), LDA );

      // Generate the last n rows of the matrix Q

      SRNAMT = 'DORGRQ';
      dorgrq(M, N, K, Q, LDA, TAU( M-K+1 ), WORK, LWORK, INFO );

      // Copy R(m-k+1:m,n-m+1:n)

      dlaset('Full', K, M, ZERO, ZERO, R( M-K+1, N-M+1 ), LDA );
      dlacpy('Upper', K, K, AF( M-K+1, N-K+1 ), LDA, R( M-K+1, N-K+1 ), LDA );

      // Compute R(m-k+1:m,n-m+1:n) - A(m-k+1:m,1:n) * Q(n-m+1:n,1:n)'

      dgemm('No transpose', 'Transpose', K, M, N, -ONE, A( M-K+1, 1 ), LDA, Q, LDA, ONE, R( M-K+1, N-M+1 ), LDA );

      // Compute norm( R - A*Q' ) / ( N * norm(A) * EPS ) .

      ANORM = DLANGE( '1', K, N, A( M-K+1, 1 ), LDA, RWORK );
      RESID = DLANGE( '1', K, M, R( M-K+1, N-M+1 ), LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = ( ( RESID / DBLE( max( 1, N ) ) ) / ANORM ) / EPS;
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute I - Q*Q'

      dlaset('Full', M, M, ZERO, ONE, R, LDA );
      dsyrk('Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = DLANSY( '1', 'Upper', M, R, LDA, RWORK );

      RESULT( 2 ) = ( RESID / DBLE( max( 1, N ) ) ) / EPS;

      return;
      }

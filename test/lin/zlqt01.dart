      void zlqt01(M, N, A, AF, Q, L, LDA, TAU, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             RESULT( * ), RWORK( * );
      Complex         A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         ROGUE;
      const              ROGUE = ( -1.0e+10, -1.0e+10 ) ;
      // ..
      // .. Local Scalars ..
      int                INFO, MINMN;
      double             ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, ZLANGE, ZLANSY;
      // EXTERNAL DLAMCH, ZLANGE, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGELQF, ZGEMM, ZHERK, ZLACPY, ZLASET, ZUNGLQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      MINMN = min( M, N );
      EPS = DLAMCH( 'Epsilon' );

      // Copy the matrix A to the array AF.

      zlacpy('Full', M, N, A, LDA, AF, LDA );

      // Factorize the matrix A in the array AF.

      SRNAMT = 'ZGELQF';
      zgelqf(M, N, AF, LDA, TAU, WORK, LWORK, INFO );

      // Copy details of Q

      zlaset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      if (N > 1) zlacpy( 'Upper', M, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA );

      // Generate the n-by-n matrix Q

      SRNAMT = 'ZUNGLQ';
      zunglq(N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy L

      zlaset('Full', M, N, DCMPLX( ZERO ), DCMPLX( ZERO ), L, LDA );
      zlacpy('Lower', M, N, AF, LDA, L, LDA );

      // Compute L - A*Q'

      zgemm('No transpose', 'Conjugate transpose', M, N, N, DCMPLX( -ONE ), A, LDA, Q, LDA, DCMPLX( ONE ), L, LDA );

      // Compute norm( L - Q'*A ) / ( N * norm(A) * EPS ) .

      ANORM = ZLANGE( '1', M, N, A, LDA, RWORK );
      RESID = ZLANGE( '1', M, N, L, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / DBLE( max( 1, N ) ) ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q*Q'

      zlaset('Full', N, N, DCMPLX( ZERO ), DCMPLX( ONE ), L, LDA );
      zherk('Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, L, LDA );

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = ZLANSY( '1', 'Upper', N, L, LDA, RWORK );

      RESULT[2] = ( RESID / DBLE( max( 1, N ) ) ) / EPS;

      return;
      }

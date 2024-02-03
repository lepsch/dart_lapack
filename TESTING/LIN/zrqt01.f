      void zrqt01(M, N, A, AF, Q, R, LDA, TAU, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double             RESULT( * ), RWORK( * );
      Complex         A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), TAU( * ), WORK( LWORK );
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
      // EXTERNAL ZGEMM, ZGERQF, ZHERK, ZLACPY, ZLASET, ZUNGRQ
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

      SRNAMT = 'ZGERQF';
      zgerqf(M, N, AF, LDA, TAU, WORK, LWORK, INFO );

      // Copy details of Q

      zlaset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      if ( M <= N ) {
         if (M > 0 && M < N) zlacpy( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA );
         IF( M > 1 ) zlacpy( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA, Q( N-M+2, N-M+1 ), LDA );
      } else {
         if (N > 1) zlacpy( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA, Q( 2, 1 ), LDA );
      }

      // Generate the n-by-n matrix Q

      SRNAMT = 'ZUNGRQ';
      zungrq(N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy R

      zlaset('Full', M, N, DCMPLX( ZERO ), DCMPLX( ZERO ), R, LDA );
      if ( M <= N ) {
         if (M > 0) zlacpy( 'Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ), LDA );
      } else {
         if (M > N && N > 0) zlacpy( 'Full', M-N, N, AF, LDA, R, LDA );
         IF( N > 0 ) zlacpy( 'Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ), LDA );
      }

      // Compute R - A*Q'

      zgemm('No transpose', 'Conjugate transpose', M, N, N, DCMPLX( -ONE ), A, LDA, Q, LDA, DCMPLX( ONE ), R, LDA );

      // Compute norm( R - Q'*A ) / ( N * norm(A) * EPS ) .

      ANORM = ZLANGE( '1', M, N, A, LDA, RWORK );
      RESID = ZLANGE( '1', M, N, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = ( ( RESID / DBLE( max( 1, N ) ) ) / ANORM ) / EPS;
      } else {
         RESULT( 1 ) = ZERO;
      }

      // Compute I - Q*Q'

      zlaset('Full', N, N, DCMPLX( ZERO ), DCMPLX( ONE ), R, LDA );
      zherk('Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = ZLANSY( '1', 'Upper', N, R, LDA, RWORK );

      RESULT( 2 ) = ( RESID / DBLE( max( 1, N ) ) ) / EPS;

      return;
      }

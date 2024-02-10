      void zqlt01(M, N, A, AF, Q, L, LDA, TAU, final Array<double> WORK, final int LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LWORK, M, N;
      double             RESULT( * ), RWORK( * );
      Complex         A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex         ROGUE;
      const              ROGUE = ( -1.0e+10, -1.0e+10 ) ;
      int                INFO, MINMN;
      double             ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- double             DLAMCH, ZLANGE, ZLANSY;
      // EXTERNAL DLAMCH, ZLANGE, ZLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZGEQLF, ZHERK, ZLACPY, ZLASET, ZUNGQL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX, MIN
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT

      MINMN = min( M, N );
      EPS = dlamch( 'Epsilon' );

      // Copy the matrix A to the array AF.

      zlacpy('Full', M, N, A, LDA, AF, LDA );

      // Factorize the matrix A in the array AF.

     srnamc.SRNAMT = 'ZGEQLF';
      zgeqlf(M, N, AF, LDA, TAU, WORK, LWORK, INFO );

      // Copy details of Q

      zlaset('Full', M, M, ROGUE, ROGUE, Q, LDA );
      if ( M >= N ) {
         if (N < M && N > 0) zlacpy( 'Full', M-N, N, AF, LDA, Q( 1, M-N+1 ), LDA );
         IF( N > 1 ) zlacpy( 'Upper', N-1, N-1, AF( M-N+1, 2 ), LDA, Q( M-N+1, M-N+2 ), LDA );
      } else {
         if (M > 1) zlacpy( 'Upper', M-1, M-1, AF( 1, N-M+2 ), LDA, Q( 1, 2 ), LDA );
      }

      // Generate the m-by-m matrix Q

     srnamc.SRNAMT = 'ZUNGQL';
      zungql(M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy L

      zlaset('Full', M, N, DCMPLX( ZERO ), DCMPLX( ZERO ), L, LDA );
      if ( M >= N ) {
         if (N > 0) zlacpy( 'Lower', N, N, AF( M-N+1, 1 ), LDA, L( M-N+1, 1 ), LDA );
      } else {
         if (N > M && M > 0) zlacpy( 'Full', M, N-M, AF, LDA, L, LDA );
         IF( M > 0 ) zlacpy( 'Lower', M, M, AF( 1, N-M+1 ), LDA, L( 1, N-M+1 ), LDA );
      }

      // Compute L - Q'*A

      zgemm('Conjugate transpose', 'No transpose', M, N, M, DCMPLX( -ONE ), Q, LDA, A, LDA, DCMPLX( ONE ), L, LDA );

      // Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = ZLANGE( '1', M, N, A, LDA, RWORK );
      RESID = ZLANGE( '1', M, N, L, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / (max( 1, M )).toDouble() ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q'*Q

      zlaset('Full', M, M, DCMPLX( ZERO ), DCMPLX( ONE ), L, LDA );
      zherk('Upper', 'Conjugate transpose', M, M, -ONE, Q, LDA, ONE, L, LDA );

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = ZLANSY( '1', 'Upper', M, L, LDA, RWORK );

      RESULT[2] = ( RESID / (max( 1, M )).toDouble() ) / EPS;

      }

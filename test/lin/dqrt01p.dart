      void dqrt01p(final int M, final int N, final int A, final int AF, final int Q, final int R, final int LDA, final int TAU, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final int RESULT,) {
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

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
      // EXTERNAL DGEMM, DGEQRFP, DLACPY, DLASET, DORGQR, DSYRK
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

     srnamc.SRNAMT = 'DGEQRFP';
      dgeqrfp(M, N, AF, LDA, TAU, WORK, LWORK, INFO );

      // Copy details of Q

      dlaset('Full', M, M, ROGUE, ROGUE, Q, LDA );
      dlacpy('Lower', M-1, N, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA );

      // Generate the m-by-m matrix Q

     srnamc.SRNAMT = 'DORGQR';
      dorgqr(M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy R

      dlaset('Full', M, N, ZERO, ZERO, R, LDA );
      dlacpy('Upper', M, N, AF, LDA, R, LDA );

      // Compute R - Q'*A

      dgemm('Transpose', 'No transpose', M, N, M, -ONE, Q, LDA, A, LDA, ONE, R, LDA );

      // Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = dlange( '1', M, N, A, LDA, RWORK );
      RESID = dlange( '1', M, N, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / (max( 1, M )).toDouble() ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q'*Q

      dlaset('Full', M, M, ZERO, ONE, R, LDA );
      dsyrk('Upper', 'Transpose', M, M, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = dlansy( '1', 'Upper', M, R, LDA, RWORK );

      RESULT[2] = ( RESID / (max( 1, M )).toDouble() ) / EPS;

      }

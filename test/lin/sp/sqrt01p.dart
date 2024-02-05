      void sqrt01p(M, N, A, AF, Q, R, LDA, TAU, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      double               A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ), WORK( LWORK );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      double               ROGUE;
      const              ROGUE = -1.0e+10 ;
      // ..
      // .. Local Scalars ..
      int                INFO, MINMN;
      double               ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGE, SLANSY;
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SGEQRFP, SLACPY, SLASET, SORGQR, SSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Executable Statements ..

      MINMN = min( M, N );
      EPS = SLAMCH( 'Epsilon' );

      // Copy the matrix A to the array AF.

      slacpy('Full', M, N, A, LDA, AF, LDA );

      // Factorize the matrix A in the array AF.

     srnamc.SRNAMT = 'SGEQRFP';
      sgeqrfp(M, N, AF, LDA, TAU, WORK, LWORK, INFO );

      // Copy details of Q

      slaset('Full', M, M, ROGUE, ROGUE, Q, LDA );
      slacpy('Lower', M-1, N, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA );

      // Generate the m-by-m matrix Q

     srnamc.SRNAMT = 'SORGQR';
      sorgqr(M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy R

      slaset('Full', M, N, ZERO, ZERO, R, LDA );
      slacpy('Upper', M, N, AF, LDA, R, LDA );

      // Compute R - Q'*A

      sgemm('Transpose', 'No transpose', M, N, M, -ONE, Q, LDA, A, LDA, ONE, R, LDA );

      // Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = SLANGE( '1', M, N, A, LDA, RWORK );
      RESID = SLANGE( '1', M, N, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / REAL( max( 1, M ) ) ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q'*Q

      slaset('Full', M, M, ZERO, ONE, R, LDA );
      ssyrk('Upper', 'Transpose', M, M, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = SLANSY( '1', 'Upper', M, R, LDA, RWORK );

      RESULT[2] = ( RESID / REAL( max( 1, M ) ) ) / EPS;

      return;
      }

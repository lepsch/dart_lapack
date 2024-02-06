      void cqrt01(M, N, A, AF, Q, R, LDA, TAU, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LWORK, M, N;
      double               RESULT( * ), RWORK( * );
      Complex            A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            ROGUE;
      const              ROGUE = ( -1.0e+10, -1.0e+10 ) ;
      int                INFO, MINMN;
      double               ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, CLANSY, SLAMCH;
      // EXTERNAL CLANGE, CLANSY, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CGEQRF, CHERK, CLACPY, CLASET, CUNGQR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT

      MINMN = min( M, N );
      EPS = SLAMCH( 'Epsilon' );

      // Copy the matrix A to the array AF.

      clacpy('Full', M, N, A, LDA, AF, LDA );

      // Factorize the matrix A in the array AF.

     srnamc.SRNAMT = 'CGEQRF';
      cgeqrf(M, N, AF, LDA, TAU, WORK, LWORK, INFO );

      // Copy details of Q

      claset('Full', M, M, ROGUE, ROGUE, Q, LDA );
      clacpy('Lower', M-1, N, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA );

      // Generate the m-by-m matrix Q

     srnamc.SRNAMT = 'CUNGQR';
      cungqr(M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy R

      claset('Full', M, N, CMPLX( ZERO ), CMPLX( ZERO ), R, LDA );
      clacpy('Upper', M, N, AF, LDA, R, LDA );

      // Compute R - Q'*A

      cgemm('Conjugate transpose', 'No transpose', M, N, M, CMPLX( -ONE ), Q, LDA, A, LDA, CMPLX( ONE ), R, LDA );

      // Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = CLANGE( '1', M, N, A, LDA, RWORK );
      RESID = CLANGE( '1', M, N, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / REAL( max( 1, M ) ) ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q'*Q

      claset('Full', M, M, CMPLX( ZERO ), CMPLX( ONE ), R, LDA );
      cherk('Upper', 'Conjugate transpose', M, M, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = CLANSY( '1', 'Upper', M, R, LDA, RWORK );

      RESULT[2] = ( RESID / REAL( max( 1, M ) ) ) / EPS;

      return;
      }

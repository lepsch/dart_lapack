      void cqrt02(M, N, K, A, AF, Q, R, LDA, TAU, WORK, LWORK, RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                K, LDA, LWORK, M, N;
      double               RESULT( * ), RWORK( * );
      Complex            A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), TAU( * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            ROGUE;
      const              ROGUE = ( -1.0e+10, -1.0e+10 ) ;
      int                INFO;
      double               ANORM, EPS, RESID;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, CLANSY, SLAMCH;
      // EXTERNAL CLANGE, CLANSY, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHERK, CLACPY, CLASET, CUNGQR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, REAL
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT

      EPS = SLAMCH( 'Epsilon' );

      // Copy the first k columns of the factorization to the array Q

      claset('Full', M, N, ROGUE, ROGUE, Q, LDA );
      clacpy('Lower', M-1, K, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA );

      // Generate the first n columns of the matrix Q

     srnamc.SRNAMT = 'CUNGQR';
      cungqr(M, N, K, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy R(1:n,1:k)

      claset('Full', N, K, CMPLX( ZERO ), CMPLX( ZERO ), R, LDA );
      clacpy('Upper', N, K, AF, LDA, R, LDA );

      // Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k)

      cgemm('Conjugate transpose', 'No transpose', N, K, M, CMPLX( -ONE ), Q, LDA, A, LDA, CMPLX( ONE ), R, LDA );

      // Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = CLANGE( '1', M, K, A, LDA, RWORK );
      RESID = CLANGE( '1', N, K, R, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / REAL( max( 1, M ) ) ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q'*Q

      claset('Full', N, N, CMPLX( ZERO ), CMPLX( ONE ), R, LDA );
      cherk('Upper', 'Conjugate transpose', N, M, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = CLANSY( '1', 'Upper', N, R, LDA, RWORK );

      RESULT[2] = ( RESID / REAL( max( 1, M ) ) ) / EPS;

      return;
      }

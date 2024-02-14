      void crqt02(final int M, final int N, final int K, final int A, final int AF, final int Q, final int R, final int LDA, final int TAU, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final int RESULT,) {
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

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
      // EXTERNAL CGEMM, CHERK, CLACPY, CLASET, CUNGRQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, REAL
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT

      // Quick return if possible

      if ( M == 0 || N == 0 || K == 0 ) {
         RESULT[1] = ZERO;
         RESULT[2] = ZERO;
         return;
      }

      EPS = SLAMCH( 'Epsilon' );

      // Copy the last k rows of the factorization to the array Q

      claset('Full', M, N, ROGUE, ROGUE, Q, LDA );
      if (K < N) clacpy( 'Full', K, N-K, AF( M-K+1, 1 ), LDA, Q( M-K+1, 1 ), LDA );
      IF( K > 1 ) clacpy( 'Lower', K-1, K-1, AF( M-K+2, N-K+1 ), LDA, Q( M-K+2, N-K+1 ), LDA );

      // Generate the last n rows of the matrix Q

     srnamc.SRNAMT = 'CUNGRQ';
      cungrq(M, N, K, Q, LDA, TAU( M-K+1 ), WORK, LWORK, INFO );

      // Copy R(m-k+1:m,n-m+1:n)

      claset('Full', K, M, CMPLX( ZERO ), CMPLX( ZERO ), R( M-K+1, N-M+1 ), LDA )       CALL CLACPY( 'Upper', K, K, AF( M-K+1, N-K+1 ), LDA, R( M-K+1, N-K+1 ), LDA );

      // Compute R(m-k+1:m,n-m+1:n) - A(m-k+1:m,1:n) * Q(n-m+1:n,1:n)'

      cgemm('No transpose', 'Conjugate transpose', K, M, N, CMPLX( -ONE ), A( M-K+1, 1 ), LDA, Q, LDA, CMPLX( ONE ), R( M-K+1, N-M+1 ), LDA );

      // Compute norm( R - A*Q' ) / ( N * norm(A) * EPS ) .

      ANORM = CLANGE( '1', K, N, A( M-K+1, 1 ), LDA, RWORK );
      RESID = CLANGE( '1', K, M, R( M-K+1, N-M+1 ), LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / REAL( max( 1, N ) ) ) / ANORM ) / EPS;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute I - Q*Q'

      claset('Full', M, M, CMPLX( ZERO ), CMPLX( ONE ), R, LDA );
      cherk('Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, R, LDA );

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = CLANSY( '1', 'Upper', M, R, LDA, RWORK );

      RESULT[2] = ( RESID / REAL( max( 1, N ) ) ) / EPS;

      }

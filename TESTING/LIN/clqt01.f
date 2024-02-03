      SUBROUTINE CLQT01( M, N, A, AF, Q, L, LDA, TAU, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               RESULT( * ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), TAU( * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0, ONE = 1.0 ;
      COMPLEX            ROGUE
      const              ROGUE = ( -1.0e+10, -1.0e+10 ) ;
      // ..
      // .. Local Scalars ..
      int                INFO, MINMN;
      REAL               ANORM, EPS, RESID
      // ..
      // .. External Functions ..
      REAL               CLANGE, CLANSY, SLAMCH
      // EXTERNAL CLANGE, CLANSY, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGELQF, CGEMM, CHERK, CLACPY, CLASET, CUNGLQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      MINMN = MIN( M, N )
      EPS = SLAMCH( 'Epsilon' )

      // Copy the matrix A to the array AF.

      clacpy('Full', M, N, A, LDA, AF, LDA );

      // Factorize the matrix A in the array AF.

      SRNAMT = 'CGELQF'
      cgelqf(M, N, AF, LDA, TAU, WORK, LWORK, INFO );

      // Copy details of Q

      claset('Full', N, N, ROGUE, ROGUE, Q, LDA );
      if (N > 1) CALL CLACPY( 'Upper', M, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA );

      // Generate the n-by-n matrix Q

      SRNAMT = 'CUNGLQ'
      cunglq(N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO );

      // Copy L

      claset('Full', M, N, CMPLX( ZERO ), CMPLX( ZERO ), L, LDA );
      clacpy('Lower', M, N, AF, LDA, L, LDA );

      // Compute L - A*Q'

      cgemm('No transpose', 'Conjugate transpose', M, N, N, CMPLX( -ONE ), A, LDA, Q, LDA, CMPLX( ONE ), L, LDA );

      // Compute norm( L - Q'*A ) / ( N * norm(A) * EPS ) .

      ANORM = CLANGE( '1', M, N, A, LDA, RWORK )
      RESID = CLANGE( '1', M, N, L, LDA, RWORK )
      if ( ANORM > ZERO ) {
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, N ) ) ) / ANORM ) / EPS
      } else {
         RESULT( 1 ) = ZERO
      }

      // Compute I - Q*Q'

      claset('Full', N, N, CMPLX( ZERO ), CMPLX( ONE ), L, LDA );
      cherk('Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, L, LDA );

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = CLANSY( '1', 'Upper', N, L, LDA, RWORK )

      RESULT( 2 ) = ( RESID / REAL( MAX( 1, N ) ) ) / EPS

      RETURN

      // End of CLQT01

      }

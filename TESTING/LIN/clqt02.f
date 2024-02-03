      SUBROUTINE CLQT02( M, N, K, A, AF, Q, L, LDA, TAU, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                K, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               RESULT( * ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDA, * ), L( LDA, * ), Q( LDA, * ), TAU( * ), WORK( LWORK )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            ROGUE
      PARAMETER          ( ROGUE = ( -1.0E+10, -1.0E+10 ) )
      // ..
      // .. Local Scalars ..
      int                INFO;
      REAL               ANORM, EPS, RESID
      // ..
      // .. External Functions ..
      REAL               CLANGE, CLANSY, SLAMCH
      // EXTERNAL CLANGE, CLANSY, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CHERK, CLACPY, CLASET, CUNGLQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, REAL
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..
*
      EPS = SLAMCH( 'Epsilon' )
*
      // Copy the first k rows of the factorization to the array Q
*
      CALL CLASET( 'Full', M, N, ROGUE, ROGUE, Q, LDA )
      CALL CLACPY( 'Upper', K, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA )
*
      // Generate the first n columns of the matrix Q
*
      SRNAMT = 'CUNGLQ'
      CALL CUNGLQ( M, N, K, Q, LDA, TAU, WORK, LWORK, INFO )
*
      // Copy L(1:k,1:m)
*
      CALL CLASET( 'Full', K, M, CMPLX( ZERO ), CMPLX( ZERO ), L, LDA )
      CALL CLACPY( 'Lower', K, M, AF, LDA, L, LDA )
*
      // Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)'
*
      CALL CGEMM( 'No transpose', 'Conjugate transpose', K, M, N, CMPLX( -ONE ), A, LDA, Q, LDA, CMPLX( ONE ), L, LDA )
*
      // Compute norm( L - A*Q' ) / ( N * norm(A) * EPS ) .
*
      ANORM = CLANGE( '1', K, N, A, LDA, RWORK )
      RESID = CLANGE( '1', K, M, L, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, N ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
      // Compute I - Q*Q'
*
      CALL CLASET( 'Full', M, M, CMPLX( ZERO ), CMPLX( ONE ), L, LDA )
      CALL CHERK( 'Upper', 'No transpose', M, N, -ONE, Q, LDA, ONE, L, LDA )
*
      // Compute norm( I - Q*Q' ) / ( N * EPS ) .
*
      RESID = CLANSY( '1', 'Upper', M, L, LDA, RWORK )
*
      RESULT( 2 ) = ( RESID / REAL( MAX( 1, N ) ) ) / EPS
*
      RETURN
*
      // End of CLQT02
*
      END

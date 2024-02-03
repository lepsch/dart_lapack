      SUBROUTINE CRQT01( M, N, A, AF, Q, R, LDA, TAU, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               RESULT( * ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDA, * ), Q( LDA, * ), R( LDA, * ), TAU( * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      COMPLEX            ROGUE
      PARAMETER          ( ROGUE = ( -1.0E+10, -1.0E+10 ) )
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
      // EXTERNAL CGEMM, CGERQF, CHERK, CLACPY, CLASET, CUNGRQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN, REAL
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      MINMN = MIN( M, N )
      EPS = SLAMCH( 'Epsilon' )

      // Copy the matrix A to the array AF.

      CALL CLACPY( 'Full', M, N, A, LDA, AF, LDA )

      // Factorize the matrix A in the array AF.

      SRNAMT = 'CGERQF'
      CALL CGERQF( M, N, AF, LDA, TAU, WORK, LWORK, INFO )

      // Copy details of Q

      CALL CLASET( 'Full', N, N, ROGUE, ROGUE, Q, LDA )
      IF( M.LE.N ) THEN
         IF( M.GT.0 .AND. M.LT.N ) CALL CLACPY( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA )          IF( M.GT.1 ) CALL CLACPY( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA, Q( N-M+2, N-M+1 ), LDA )
      ELSE
         IF( N.GT.1 ) CALL CLACPY( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA, Q( 2, 1 ), LDA )
      END IF

      // Generate the n-by-n matrix Q

      SRNAMT = 'CUNGRQ'
      CALL CUNGRQ( N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO )

      // Copy R

      CALL CLASET( 'Full', M, N, CMPLX( ZERO ), CMPLX( ZERO ), R, LDA )
      IF( M.LE.N ) THEN
         IF( M.GT.0 ) CALL CLACPY( 'Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ), LDA )
      ELSE
         IF( M.GT.N .AND. N.GT.0 ) CALL CLACPY( 'Full', M-N, N, AF, LDA, R, LDA )          IF( N.GT.0 ) CALL CLACPY( 'Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ), LDA )
      END IF

      // Compute R - A*Q'

      CALL CGEMM( 'No transpose', 'Conjugate transpose', M, N, N, CMPLX( -ONE ), A, LDA, Q, LDA, CMPLX( ONE ), R, LDA )

      // Compute norm( R - Q'*A ) / ( N * norm(A) * EPS ) .

      ANORM = CLANGE( '1', M, N, A, LDA, RWORK )
      RESID = CLANGE( '1', M, N, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, N ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF

      // Compute I - Q*Q'

      CALL CLASET( 'Full', N, N, CMPLX( ZERO ), CMPLX( ONE ), R, LDA )
      CALL CHERK( 'Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, R, LDA )

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = CLANSY( '1', 'Upper', N, R, LDA, RWORK )

      RESULT( 2 ) = ( RESID / REAL( MAX( 1, N ) ) ) / EPS

      RETURN

      // End of CRQT01

      END

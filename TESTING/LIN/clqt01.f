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
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            ROGUE
      const              ROGUE = ( -1.0E+10, -1.0E+10 ) ;
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
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Executable Statements ..

      MINMN = MIN( M, N )
      EPS = SLAMCH( 'Epsilon' )

      // Copy the matrix A to the array AF.

      CALL CLACPY( 'Full', M, N, A, LDA, AF, LDA )

      // Factorize the matrix A in the array AF.

      SRNAMT = 'CGELQF'
      CALL CGELQF( M, N, AF, LDA, TAU, WORK, LWORK, INFO )

      // Copy details of Q

      CALL CLASET( 'Full', N, N, ROGUE, ROGUE, Q, LDA )
      IF( N.GT.1 ) CALL CLACPY( 'Upper', M, N-1, AF( 1, 2 ), LDA, Q( 1, 2 ), LDA )

      // Generate the n-by-n matrix Q

      SRNAMT = 'CUNGLQ'
      CALL CUNGLQ( N, N, MINMN, Q, LDA, TAU, WORK, LWORK, INFO )

      // Copy L

      CALL CLASET( 'Full', M, N, CMPLX( ZERO ), CMPLX( ZERO ), L, LDA )
      CALL CLACPY( 'Lower', M, N, AF, LDA, L, LDA )

      // Compute L - A*Q'

      CALL CGEMM( 'No transpose', 'Conjugate transpose', M, N, N, CMPLX( -ONE ), A, LDA, Q, LDA, CMPLX( ONE ), L, LDA )

      // Compute norm( L - Q'*A ) / ( N * norm(A) * EPS ) .

      ANORM = CLANGE( '1', M, N, A, LDA, RWORK )
      RESID = CLANGE( '1', M, N, L, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, N ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF

      // Compute I - Q*Q'

      CALL CLASET( 'Full', N, N, CMPLX( ZERO ), CMPLX( ONE ), L, LDA )
      CALL CHERK( 'Upper', 'No transpose', N, N, -ONE, Q, LDA, ONE, L, LDA )

      // Compute norm( I - Q*Q' ) / ( N * EPS ) .

      RESID = CLANSY( '1', 'Upper', N, L, LDA, RWORK )

      RESULT( 2 ) = ( RESID / REAL( MAX( 1, N ) ) ) / EPS

      RETURN

      // End of CLQT01

      }

      SUBROUTINE CQLT01( M, N, A, AF, Q, L, LDA, TAU, WORK, LWORK, RWORK, RESULT )

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
      // EXTERNAL CGEMM, CGEQLF, CHERK, CLACPY, CLASET, CUNGQL
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

      SRNAMT = 'CGEQLF'
      CALL CGEQLF( M, N, AF, LDA, TAU, WORK, LWORK, INFO )

      // Copy details of Q

      CALL CLASET( 'Full', M, M, ROGUE, ROGUE, Q, LDA )
      IF( M.GE.N ) THEN
         IF( N.LT.M .AND. N.GT.0 ) CALL CLACPY( 'Full', M-N, N, AF, LDA, Q( 1, M-N+1 ), LDA )          IF( N.GT.1 ) CALL CLACPY( 'Upper', N-1, N-1, AF( M-N+1, 2 ), LDA, Q( M-N+1, M-N+2 ), LDA )
      ELSE
         IF( M.GT.1 ) CALL CLACPY( 'Upper', M-1, M-1, AF( 1, N-M+2 ), LDA, Q( 1, 2 ), LDA )
      END IF

      // Generate the m-by-m matrix Q

      SRNAMT = 'CUNGQL'
      CALL CUNGQL( M, M, MINMN, Q, LDA, TAU, WORK, LWORK, INFO )

      // Copy L

      CALL CLASET( 'Full', M, N, CMPLX( ZERO ), CMPLX( ZERO ), L, LDA )
      IF( M.GE.N ) THEN
         IF( N.GT.0 ) CALL CLACPY( 'Lower', N, N, AF( M-N+1, 1 ), LDA, L( M-N+1, 1 ), LDA )
      ELSE
         IF( N.GT.M .AND. M.GT.0 ) CALL CLACPY( 'Full', M, N-M, AF, LDA, L, LDA )          IF( M.GT.0 ) CALL CLACPY( 'Lower', M, M, AF( 1, N-M+1 ), LDA, L( 1, N-M+1 ), LDA )
      END IF

      // Compute L - Q'*A

      CALL CGEMM( 'Conjugate transpose', 'No transpose', M, N, M, CMPLX( -ONE ), Q, LDA, A, LDA, CMPLX( ONE ), L, LDA )

      // Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .

      ANORM = CLANGE( '1', M, N, A, LDA, RWORK )
      RESID = CLANGE( '1', M, N, L, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, M ) ) ) / ANORM ) / EPS
      ELSE
         RESULT( 1 ) = ZERO
      END IF

      // Compute I - Q'*Q

      CALL CLASET( 'Full', M, M, CMPLX( ZERO ), CMPLX( ONE ), L, LDA )
      CALL CHERK( 'Upper', 'Conjugate transpose', M, M, -ONE, Q, LDA, ONE, L, LDA )

      // Compute norm( I - Q'*Q ) / ( M * EPS ) .

      RESID = CLANSY( '1', 'Upper', M, L, LDA, RWORK )

      RESULT( 2 ) = ( RESID / REAL( MAX( 1, M ) ) ) / EPS

      RETURN

      // End of CQLT01

      }

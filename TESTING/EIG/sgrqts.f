      SUBROUTINE SGRQTS( M, P, N, A, AF, Q, R, LDA, TAUA, B, BF, Z, T, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, P, N;
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), AF( LDA, * ), R( LDA, * ), Q( LDA, * ), B( LDB, * ), BF( LDB, * ), T( LDB, * ), Z( LDB, * ), BWK( LDB, * ), TAUA( * ), TAUB( * ), RESULT( 4 ), RWORK( * ), WORK( LWORK )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      REAL               ROGUE
      PARAMETER          ( ROGUE = -1.0E+10 )
*     ..
*     .. Local Scalars ..
      int                INFO;
      REAL               ANORM, BNORM, ULP, UNFL, RESID
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE, SLANSY
      EXTERNAL           SLAMCH, SLANGE, SLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMM, SGGRQF, SLACPY, SLASET, SORGQR, SORGRQ, SSYRK
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      ULP = SLAMCH( 'Precision' )
      UNFL = SLAMCH( 'Safe minimum' )
*
*     Copy the matrix A to the array AF.
*
      CALL SLACPY( 'Full', M, N, A, LDA, AF, LDA )
      CALL SLACPY( 'Full', P, N, B, LDB, BF, LDB )
*
      ANORM = MAX( SLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
      BNORM = MAX( SLANGE( '1', P, N, B, LDB, RWORK ), UNFL )
*
*     Factorize the matrices A and B in the arrays AF and BF.
*
      CALL SGGRQF( M, P, N, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO )
*
*     Generate the N-by-N matrix Q
*
      CALL SLASET( 'Full', N, N, ROGUE, ROGUE, Q, LDA )
      IF( M.LE.N ) THEN
         IF( M.GT.0 .AND. M.LT.N ) CALL SLACPY( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA )          IF( M.GT.1 ) CALL SLACPY( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA, Q( N-M+2, N-M+1 ), LDA )
      ELSE
         IF( N.GT.1 ) CALL SLACPY( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA, Q( 2, 1 ), LDA )
      END IF
      CALL SORGRQ( N, N, MIN( M, N ), Q, LDA, TAUA, WORK, LWORK, INFO )
*
*     Generate the P-by-P matrix Z
*
      CALL SLASET( 'Full', P, P, ROGUE, ROGUE, Z, LDB )
      IF( P.GT.1 ) CALL SLACPY( 'Lower', P-1, N, BF( 2,1 ), LDB, Z( 2,1 ), LDB )
      CALL SORGQR( P, P, MIN( P,N ), Z, LDB, TAUB, WORK, LWORK, INFO )
*
*     Copy R
*
      CALL SLASET( 'Full', M, N, ZERO, ZERO, R, LDA )
      IF( M.LE.N )THEN
         CALL SLACPY( 'Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ), LDA )
      ELSE
         CALL SLACPY( 'Full', M-N, N, AF, LDA, R, LDA )
         CALL SLACPY( 'Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ), LDA )
      END IF
*
*     Copy T
*
      CALL SLASET( 'Full', P, N, ZERO, ZERO, T, LDB )
      CALL SLACPY( 'Upper', P, N, BF, LDB, T, LDB )
*
*     Compute R - A*Q'
*
      CALL SGEMM( 'No transpose', 'Transpose', M, N, N, -ONE, A, LDA, Q, LDA, ONE, R, LDA )
*
*     Compute norm( R - A*Q' ) / ( MAX(M,N)*norm(A)*ULP ) .
*
      RESID = SLANGE( '1', M, N, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / REAL(MAX(1,M,N) ) ) / ANORM ) / ULP
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute T*Q - Z'*B
*
      CALL SGEMM( 'Transpose', 'No transpose', P, N, P, ONE, Z, LDB, B, LDB, ZERO, BWK, LDB )       CALL SGEMM( 'No transpose', 'No transpose', P, N, N, ONE, T, LDB, Q, LDA, -ONE, BWK, LDB )
*
*     Compute norm( T*Q - Z'*B ) / ( MAX(P,N)*norm(A)*ULP ) .
*
      RESID = SLANGE( '1', P, N, BWK, LDB, RWORK )
      IF( BNORM.GT.ZERO ) THEN
         RESULT( 2 ) = ( ( RESID / REAL( MAX( 1,P,M ) ) )/BNORM ) / ULP
      ELSE
         RESULT( 2 ) = ZERO
      END IF
*
*     Compute I - Q*Q'
*
      CALL SLASET( 'Full', N, N, ZERO, ONE, R, LDA )
      CALL SSYRK( 'Upper', 'No Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA )
*
*     Compute norm( I - Q'*Q ) / ( N * ULP ) .
*
      RESID = SLANSY( '1', 'Upper', N, R, LDA, RWORK )
      RESULT( 3 ) = ( RESID / REAL( MAX( 1,N ) ) ) / ULP
*
*     Compute I - Z'*Z
*
      CALL SLASET( 'Full', P, P, ZERO, ONE, T, LDB )
      CALL SSYRK( 'Upper', 'Transpose', P, P, -ONE, Z, LDB, ONE, T, LDB )
*
*     Compute norm( I - Z'*Z ) / ( P*ULP ) .
*
      RESID = SLANSY( '1', 'Upper', P, T, LDB, RWORK )
      RESULT( 4 ) = ( RESID / REAL( MAX( 1,P ) ) ) / ULP
*
      RETURN
*
*     End of SGRQTS
*
      END

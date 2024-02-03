      SUBROUTINE DGRQTS( M, P, N, A, AF, Q, R, LDA, TAUA, B, BF, Z, T, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, N, P;
*     ..
*     .. Array Arguments ..
      double             A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), BWK( LDB, * ), Q( LDA, * ), R( LDA, * ), RESULT( 4 ), RWORK( * ), T( LDB, * ), TAUA( * ), TAUB( * ), WORK( LWORK ), Z( LDB, * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      double             ROGUE;
      PARAMETER          ( ROGUE = -1.0D+10 )
*     ..
*     .. Local Scalars ..
      int                INFO;
      double             ANORM, BNORM, RESID, ULP, UNFL;
*     ..
*     .. External Functions ..
      double             DLAMCH, DLANGE, DLANSY;
      EXTERNAL           DLAMCH, DLANGE, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGGRQF, DLACPY, DLASET, DORGQR, DORGRQ, DSYRK
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      ULP = DLAMCH( 'Precision' )
      UNFL = DLAMCH( 'Safe minimum' )
*
*     Copy the matrix A to the array AF.
*
      CALL DLACPY( 'Full', M, N, A, LDA, AF, LDA )
      CALL DLACPY( 'Full', P, N, B, LDB, BF, LDB )
*
      ANORM = MAX( DLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
      BNORM = MAX( DLANGE( '1', P, N, B, LDB, RWORK ), UNFL )
*
*     Factorize the matrices A and B in the arrays AF and BF.
*
      CALL DGGRQF( M, P, N, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO )
*
*     Generate the N-by-N matrix Q
*
      CALL DLASET( 'Full', N, N, ROGUE, ROGUE, Q, LDA )
      IF( M.LE.N ) THEN
         IF( M.GT.0 .AND. M.LT.N ) CALL DLACPY( 'Full', M, N-M, AF, LDA, Q( N-M+1, 1 ), LDA )          IF( M.GT.1 ) CALL DLACPY( 'Lower', M-1, M-1, AF( 2, N-M+1 ), LDA, Q( N-M+2, N-M+1 ), LDA )
      ELSE
         IF( N.GT.1 ) CALL DLACPY( 'Lower', N-1, N-1, AF( M-N+2, 1 ), LDA, Q( 2, 1 ), LDA )
      END IF
      CALL DORGRQ( N, N, MIN( M, N ), Q, LDA, TAUA, WORK, LWORK, INFO )
*
*     Generate the P-by-P matrix Z
*
      CALL DLASET( 'Full', P, P, ROGUE, ROGUE, Z, LDB )
      IF( P.GT.1 ) CALL DLACPY( 'Lower', P-1, N, BF( 2, 1 ), LDB, Z( 2, 1 ), LDB )
      CALL DORGQR( P, P, MIN( P, N ), Z, LDB, TAUB, WORK, LWORK, INFO )
*
*     Copy R
*
      CALL DLASET( 'Full', M, N, ZERO, ZERO, R, LDA )
      IF( M.LE.N ) THEN
         CALL DLACPY( 'Upper', M, M, AF( 1, N-M+1 ), LDA, R( 1, N-M+1 ), LDA )
      ELSE
         CALL DLACPY( 'Full', M-N, N, AF, LDA, R, LDA )
         CALL DLACPY( 'Upper', N, N, AF( M-N+1, 1 ), LDA, R( M-N+1, 1 ), LDA )
      END IF
*
*     Copy T
*
      CALL DLASET( 'Full', P, N, ZERO, ZERO, T, LDB )
      CALL DLACPY( 'Upper', P, N, BF, LDB, T, LDB )
*
*     Compute R - A*Q'
*
      CALL DGEMM( 'No transpose', 'Transpose', M, N, N, -ONE, A, LDA, Q, LDA, ONE, R, LDA )
*
*     Compute norm( R - A*Q' ) / ( MAX(M,N)*norm(A)*ULP ) .
*
      RESID = DLANGE( '1', M, N, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M, N ) ) ) / ANORM ) / ULP
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute T*Q - Z'*B
*
      CALL DGEMM( 'Transpose', 'No transpose', P, N, P, ONE, Z, LDB, B, LDB, ZERO, BWK, LDB )       CALL DGEMM( 'No transpose', 'No transpose', P, N, N, ONE, T, LDB, Q, LDA, -ONE, BWK, LDB )
*
*     Compute norm( T*Q - Z'*B ) / ( MAX(P,N)*norm(A)*ULP ) .
*
      RESID = DLANGE( '1', P, N, BWK, LDB, RWORK )
      IF( BNORM.GT.ZERO ) THEN
         RESULT( 2 ) = ( ( RESID / DBLE( MAX( 1, P, M ) ) ) / BNORM ) / ULP
      ELSE
         RESULT( 2 ) = ZERO
      END IF
*
*     Compute I - Q*Q'
*
      CALL DLASET( 'Full', N, N, ZERO, ONE, R, LDA )
      CALL DSYRK( 'Upper', 'No Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA )
*
*     Compute norm( I - Q'*Q ) / ( N * ULP ) .
*
      RESID = DLANSY( '1', 'Upper', N, R, LDA, RWORK )
      RESULT( 3 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / ULP
*
*     Compute I - Z'*Z
*
      CALL DLASET( 'Full', P, P, ZERO, ONE, T, LDB )
      CALL DSYRK( 'Upper', 'Transpose', P, P, -ONE, Z, LDB, ONE, T, LDB )
*
*     Compute norm( I - Z'*Z ) / ( P*ULP ) .
*
      RESID = DLANSY( '1', 'Upper', P, T, LDB, RWORK )
      RESULT( 4 ) = ( RESID / DBLE( MAX( 1, P ) ) ) / ULP
*
      RETURN
*
*     End of DGRQTS
*
      END

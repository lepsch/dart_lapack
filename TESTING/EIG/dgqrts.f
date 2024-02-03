      SUBROUTINE DGQRTS( N, M, P, A, AF, Q, R, LDA, TAUA, B, BF, Z, T, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )
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
      // EXTERNAL DLAMCH, DLANGE, DLANSY
*     ..
*     .. External Subroutines ..
      // EXTERNAL DGEMM, DGGQRF, DLACPY, DLASET, DORGQR, DORGRQ, DSYRK
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
      CALL DLACPY( 'Full', N, M, A, LDA, AF, LDA )
      CALL DLACPY( 'Full', N, P, B, LDB, BF, LDB )
*
      ANORM = MAX( DLANGE( '1', N, M, A, LDA, RWORK ), UNFL )
      BNORM = MAX( DLANGE( '1', N, P, B, LDB, RWORK ), UNFL )
*
*     Factorize the matrices A and B in the arrays AF and BF.
*
      CALL DGGQRF( N, M, P, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO )
*
*     Generate the N-by-N matrix Q
*
      CALL DLASET( 'Full', N, N, ROGUE, ROGUE, Q, LDA )
      CALL DLACPY( 'Lower', N-1, M, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA )
      CALL DORGQR( N, N, MIN( N, M ), Q, LDA, TAUA, WORK, LWORK, INFO )
*
*     Generate the P-by-P matrix Z
*
      CALL DLASET( 'Full', P, P, ROGUE, ROGUE, Z, LDB )
      IF( N.LE.P ) THEN
         IF( N.GT.0 .AND. N.LT.P ) CALL DLACPY( 'Full', N, P-N, BF, LDB, Z( P-N+1, 1 ), LDB )          IF( N.GT.1 ) CALL DLACPY( 'Lower', N-1, N-1, BF( 2, P-N+1 ), LDB, Z( P-N+2, P-N+1 ), LDB )
      ELSE
         IF( P.GT.1 ) CALL DLACPY( 'Lower', P-1, P-1, BF( N-P+2, 1 ), LDB, Z( 2, 1 ), LDB )
      END IF
      CALL DORGRQ( P, P, MIN( N, P ), Z, LDB, TAUB, WORK, LWORK, INFO )
*
*     Copy R
*
      CALL DLASET( 'Full', N, M, ZERO, ZERO, R, LDA )
      CALL DLACPY( 'Upper', N, M, AF, LDA, R, LDA )
*
*     Copy T
*
      CALL DLASET( 'Full', N, P, ZERO, ZERO, T, LDB )
      IF( N.LE.P ) THEN
         CALL DLACPY( 'Upper', N, N, BF( 1, P-N+1 ), LDB, T( 1, P-N+1 ), LDB )
      ELSE
         CALL DLACPY( 'Full', N-P, P, BF, LDB, T, LDB )
         CALL DLACPY( 'Upper', P, P, BF( N-P+1, 1 ), LDB, T( N-P+1, 1 ), LDB )
      END IF
*
*     Compute R - Q'*A
*
      CALL DGEMM( 'Transpose', 'No transpose', N, M, N, -ONE, Q, LDA, A, LDA, ONE, R, LDA )
*
*     Compute norm( R - Q'*A ) / ( MAX(M,N)*norm(A)*ULP ) .
*
      RESID = DLANGE( '1', N, M, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M, N ) ) ) / ANORM ) / ULP
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     Compute T*Z - Q'*B
*
      CALL DGEMM( 'No Transpose', 'No transpose', N, P, P, ONE, T, LDB, Z, LDB, ZERO, BWK, LDB )       CALL DGEMM( 'Transpose', 'No transpose', N, P, N, -ONE, Q, LDA, B, LDB, ONE, BWK, LDB )
*
*     Compute norm( T*Z - Q'*B ) / ( MAX(P,N)*norm(A)*ULP ) .
*
      RESID = DLANGE( '1', N, P, BWK, LDB, RWORK )
      IF( BNORM.GT.ZERO ) THEN
         RESULT( 2 ) = ( ( RESID / DBLE( MAX( 1, P, N ) ) ) / BNORM ) / ULP
      ELSE
         RESULT( 2 ) = ZERO
      END IF
*
*     Compute I - Q'*Q
*
      CALL DLASET( 'Full', N, N, ZERO, ONE, R, LDA )
      CALL DSYRK( 'Upper', 'Transpose', N, N, -ONE, Q, LDA, ONE, R, LDA )
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
*     End of DGQRTS
*
      END

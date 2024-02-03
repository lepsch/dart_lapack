      SUBROUTINE ZGQRTS( N, M, P, A, AF, Q, R, LDA, TAUA, B, BF, Z, T, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      double             RESULT( 4 ), RWORK( * );
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), BWK( LDB, * ), Q( LDA, * ), R( LDA, * ), T( LDB, * ), TAUA( * ), TAUB( * ), WORK( LWORK ), Z( LDB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         CROGUE
      PARAMETER          ( CROGUE = ( -1.0D+10, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      int                INFO;
      double             ANORM, BNORM, RESID, ULP, UNFL;
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE, ZLANHE;
      // EXTERNAL DLAMCH, ZLANGE, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZGGQRF, ZHERK, ZLACPY, ZLASET, ZUNGQR, ZUNGRQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      ULP = DLAMCH( 'Precision' )
      UNFL = DLAMCH( 'Safe minimum' )

      // Copy the matrix A to the array AF.

      CALL ZLACPY( 'Full', N, M, A, LDA, AF, LDA )
      CALL ZLACPY( 'Full', N, P, B, LDB, BF, LDB )

      ANORM = MAX( ZLANGE( '1', N, M, A, LDA, RWORK ), UNFL )
      BNORM = MAX( ZLANGE( '1', N, P, B, LDB, RWORK ), UNFL )

      // Factorize the matrices A and B in the arrays AF and BF.

      CALL ZGGQRF( N, M, P, AF, LDA, TAUA, BF, LDB, TAUB, WORK, LWORK, INFO )

      // Generate the N-by-N matrix Q

      CALL ZLASET( 'Full', N, N, CROGUE, CROGUE, Q, LDA )
      CALL ZLACPY( 'Lower', N-1, M, AF( 2, 1 ), LDA, Q( 2, 1 ), LDA )
      CALL ZUNGQR( N, N, MIN( N, M ), Q, LDA, TAUA, WORK, LWORK, INFO )

      // Generate the P-by-P matrix Z

      CALL ZLASET( 'Full', P, P, CROGUE, CROGUE, Z, LDB )
      IF( N.LE.P ) THEN
         IF( N.GT.0 .AND. N.LT.P ) CALL ZLACPY( 'Full', N, P-N, BF, LDB, Z( P-N+1, 1 ), LDB )          IF( N.GT.1 ) CALL ZLACPY( 'Lower', N-1, N-1, BF( 2, P-N+1 ), LDB, Z( P-N+2, P-N+1 ), LDB )
      ELSE
         IF( P.GT.1 ) CALL ZLACPY( 'Lower', P-1, P-1, BF( N-P+2, 1 ), LDB, Z( 2, 1 ), LDB )
      END IF
      CALL ZUNGRQ( P, P, MIN( N, P ), Z, LDB, TAUB, WORK, LWORK, INFO )

      // Copy R

      CALL ZLASET( 'Full', N, M, CZERO, CZERO, R, LDA )
      CALL ZLACPY( 'Upper', N, M, AF, LDA, R, LDA )

      // Copy T

      CALL ZLASET( 'Full', N, P, CZERO, CZERO, T, LDB )
      IF( N.LE.P ) THEN
         CALL ZLACPY( 'Upper', N, N, BF( 1, P-N+1 ), LDB, T( 1, P-N+1 ), LDB )
      ELSE
         CALL ZLACPY( 'Full', N-P, P, BF, LDB, T, LDB )
         CALL ZLACPY( 'Upper', P, P, BF( N-P+1, 1 ), LDB, T( N-P+1, 1 ), LDB )
      END IF

      // Compute R - Q'*A

      CALL ZGEMM( 'Conjugate transpose', 'No transpose', N, M, N, -CONE, Q, LDA, A, LDA, CONE, R, LDA )

      // Compute norm( R - Q'*A ) / ( MAX(M,N)*norm(A)*ULP ) .

      RESID = ZLANGE( '1', N, M, R, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M, N ) ) ) / ANORM ) / ULP
      ELSE
         RESULT( 1 ) = ZERO
      END IF

      // Compute T*Z - Q'*B

      CALL ZGEMM( 'No Transpose', 'No transpose', N, P, P, CONE, T, LDB, Z, LDB, CZERO, BWK, LDB )       CALL ZGEMM( 'Conjugate transpose', 'No transpose', N, P, N, -CONE, Q, LDA, B, LDB, CONE, BWK, LDB )

      // Compute norm( T*Z - Q'*B ) / ( MAX(P,N)*norm(A)*ULP ) .

      RESID = ZLANGE( '1', N, P, BWK, LDB, RWORK )
      IF( BNORM.GT.ZERO ) THEN
         RESULT( 2 ) = ( ( RESID / DBLE( MAX( 1, P, N ) ) ) / BNORM ) / ULP
      ELSE
         RESULT( 2 ) = ZERO
      END IF

      // Compute I - Q'*Q

      CALL ZLASET( 'Full', N, N, CZERO, CONE, R, LDA )
      CALL ZHERK( 'Upper', 'Conjugate transpose', N, N, -ONE, Q, LDA, ONE, R, LDA )

      // Compute norm( I - Q'*Q ) / ( N * ULP ) .

      RESID = ZLANHE( '1', 'Upper', N, R, LDA, RWORK )
      RESULT( 3 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / ULP

      // Compute I - Z'*Z

      CALL ZLASET( 'Full', P, P, CZERO, CONE, T, LDB )
      CALL ZHERK( 'Upper', 'Conjugate transpose', P, P, -ONE, Z, LDB, ONE, T, LDB )

      // Compute norm( I - Z'*Z ) / ( P*ULP ) .

      RESID = ZLANHE( '1', 'Upper', P, T, LDB, RWORK )
      RESULT( 4 ) = ( RESID / DBLE( MAX( 1, P ) ) ) / ULP

      RETURN

      // End of ZGQRTS

      END

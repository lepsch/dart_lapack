      SUBROUTINE ZGSVTS3( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V, LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK, LWORK, RWORK, RESULT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             ALPHA( * ), BETA( * ), RESULT( 6 ), RWORK( * );
      COMPLEX*16         A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), Q( LDQ, * ), R( LDR, * ), U( LDU, * ), V( LDV, * ), WORK( LWORK )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ZERO, ONE;
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) )
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, K, L;
      double             ANORM, BNORM, RESID, TEMP, ULP, ULPINV, UNFL;
      // ..
      // .. External Functions ..
      double             DLAMCH, ZLANGE, ZLANHE;
      // EXTERNAL DLAMCH, ZLANGE, ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, ZGEMM, ZGGSVD3, ZHERK, ZLACPY, ZLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..
*
      ULP = DLAMCH( 'Precision' )
      ULPINV = ONE / ULP
      UNFL = DLAMCH( 'Safe minimum' )
*
      // Copy the matrix A to the array AF.
*
      CALL ZLACPY( 'Full', M, N, A, LDA, AF, LDA )
      CALL ZLACPY( 'Full', P, N, B, LDB, BF, LDB )
*
      ANORM = MAX( ZLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
      BNORM = MAX( ZLANGE( '1', P, N, B, LDB, RWORK ), UNFL )
*
      // Factorize the matrices A and B in the arrays AF and BF.
*
      CALL ZGGSVD3( 'U', 'V', 'Q', M, N, P, K, L, AF, LDA, BF, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, LWORK, RWORK, IWORK, INFO )
*
      // Copy R
*
      DO 20 I = 1, MIN( K+L, M )
         DO 10 J = I, K + L
            R( I, J ) = AF( I, N-K-L+J )
   10    CONTINUE
   20 CONTINUE
*
      IF( M-K-L.LT.0 ) THEN
         DO 40 I = M + 1, K + L
            DO 30 J = I, K + L
               R( I, J ) = BF( I-K, N-K-L+J )
   30       CONTINUE
   40    CONTINUE
      END IF
*
      // Compute A:= U'*A*Q - D1*R
*
      CALL ZGEMM( 'No transpose', 'No transpose', M, N, N, CONE, A, LDA, Q, LDQ, CZERO, WORK, LDA )
*
      CALL ZGEMM( 'Conjugate transpose', 'No transpose', M, N, M, CONE, U, LDU, WORK, LDA, CZERO, A, LDA )
*
      DO 60 I = 1, K
         DO 50 J = I, K + L
            A( I, N-K-L+J ) = A( I, N-K-L+J ) - R( I, J )
   50    CONTINUE
   60 CONTINUE
*
      DO 80 I = K + 1, MIN( K+L, M )
         DO 70 J = I, K + L
            A( I, N-K-L+J ) = A( I, N-K-L+J ) - ALPHA( I )*R( I, J )
   70    CONTINUE
   80 CONTINUE
*
      // Compute norm( U'*A*Q - D1*R ) / ( MAX(1,M,N)*norm(A)*ULP ) .
*
      RESID = ZLANGE( '1', M, N, A, LDA, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M, N ) ) ) / ANORM ) / ULP
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
      // Compute B := V'*B*Q - D2*R
*
      CALL ZGEMM( 'No transpose', 'No transpose', P, N, N, CONE, B, LDB, Q, LDQ, CZERO, WORK, LDB )
*
      CALL ZGEMM( 'Conjugate transpose', 'No transpose', P, N, P, CONE, V, LDV, WORK, LDB, CZERO, B, LDB )
*
      DO 100 I = 1, L
         DO 90 J = I, L
            B( I, N-L+J ) = B( I, N-L+J ) - BETA( K+I )*R( K+I, K+J )
   90    CONTINUE
  100 CONTINUE
*
      // Compute norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP ) .
*
      RESID = ZLANGE( '1', P, N, B, LDB, RWORK )
      IF( BNORM.GT.ZERO ) THEN
         RESULT( 2 ) = ( ( RESID / DBLE( MAX( 1, P, N ) ) ) / BNORM ) / ULP
      ELSE
         RESULT( 2 ) = ZERO
      END IF
*
      // Compute I - U'*U
*
      CALL ZLASET( 'Full', M, M, CZERO, CONE, WORK, LDQ )
      CALL ZHERK( 'Upper', 'Conjugate transpose', M, M, -ONE, U, LDU, ONE, WORK, LDU )
*
      // Compute norm( I - U'*U ) / ( M * ULP ) .
*
      RESID = ZLANHE( '1', 'Upper', M, WORK, LDU, RWORK )
      RESULT( 3 ) = ( RESID / DBLE( MAX( 1, M ) ) ) / ULP
*
      // Compute I - V'*V
*
      CALL ZLASET( 'Full', P, P, CZERO, CONE, WORK, LDV )
      CALL ZHERK( 'Upper', 'Conjugate transpose', P, P, -ONE, V, LDV, ONE, WORK, LDV )
*
      // Compute norm( I - V'*V ) / ( P * ULP ) .
*
      RESID = ZLANHE( '1', 'Upper', P, WORK, LDV, RWORK )
      RESULT( 4 ) = ( RESID / DBLE( MAX( 1, P ) ) ) / ULP
*
      // Compute I - Q'*Q
*
      CALL ZLASET( 'Full', N, N, CZERO, CONE, WORK, LDQ )
      CALL ZHERK( 'Upper', 'Conjugate transpose', N, N, -ONE, Q, LDQ, ONE, WORK, LDQ )
*
      // Compute norm( I - Q'*Q ) / ( N * ULP ) .
*
      RESID = ZLANHE( '1', 'Upper', N, WORK, LDQ, RWORK )
      RESULT( 5 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / ULP
*
      // Check sorting
*
      CALL DCOPY( N, ALPHA, 1, RWORK, 1 )
      DO 110 I = K + 1, MIN( K+L, M )
         J = IWORK( I )
         IF( I.NE.J ) THEN
            TEMP = RWORK( I )
            RWORK( I ) = RWORK( J )
            RWORK( J ) = TEMP
         END IF
  110 CONTINUE
*
      RESULT( 6 ) = ZERO
      DO 120 I = K + 1, MIN( K+L, M ) - 1
         IF( RWORK( I ).LT.RWORK( I+1 ) ) RESULT( 6 ) = ULPINV
  120 CONTINUE
*
      RETURN
*
      // End of ZGSVTS3
*
      END

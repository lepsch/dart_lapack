      SUBROUTINE SGSVTS3( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V, LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               A( LDA, * ), AF( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), BF( LDB, * ), Q( LDQ, * ), R( LDR, * ), RESULT( 6 ), RWORK( * ), U( LDU, * ), V( LDV, * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, K, L;
      REAL               ANORM, BNORM, RESID, TEMP, ULP, ULPINV, UNFL
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLANGE, SLANSY
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMM, SGGSVD3, SLACPY, SLASET, SSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      ULP = SLAMCH( 'Precision' )
      ULPINV = ONE / ULP
      UNFL = SLAMCH( 'Safe minimum' )

      // Copy the matrix A to the array AF.

      slacpy('Full', M, N, A, LDA, AF, LDA );
      slacpy('Full', P, N, B, LDB, BF, LDB );

      ANORM = MAX( SLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
      BNORM = MAX( SLANGE( '1', P, N, B, LDB, RWORK ), UNFL )

      // Factorize the matrices A and B in the arrays AF and BF.

      sggsvd3('U', 'V', 'Q', M, N, P, K, L, AF, LDA, BF, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, LWORK, IWORK, INFO );

      // Copy R

      DO 20 I = 1, MIN( K+L, M )
         DO 10 J = I, K + L
            R( I, J ) = AF( I, N-K-L+J )
   10    CONTINUE
   20 CONTINUE

      if ( M-K-L.LT.0 ) {
         DO 40 I = M + 1, K + L
            DO 30 J = I, K + L
               R( I, J ) = BF( I-K, N-K-L+J )
   30       CONTINUE
   40    CONTINUE
      }

      // Compute A:= U'*A*Q - D1*R

      sgemm('No transpose', 'No transpose', M, N, N, ONE, A, LDA, Q, LDQ, ZERO, WORK, LDA );

      sgemm('Transpose', 'No transpose', M, N, M, ONE, U, LDU, WORK, LDA, ZERO, A, LDA );

      for (I = 1; I <= K; I++) { // 60
         DO 50 J = I, K + L
            A( I, N-K-L+J ) = A( I, N-K-L+J ) - R( I, J )
   50    CONTINUE
   60 CONTINUE

      DO 80 I = K + 1, MIN( K+L, M )
         DO 70 J = I, K + L
            A( I, N-K-L+J ) = A( I, N-K-L+J ) - ALPHA( I )*R( I, J )
   70    CONTINUE
   80 CONTINUE

      // Compute norm( U'*A*Q - D1*R ) / ( MAX(1,M,N)*norm(A)*ULP ) .

      RESID = SLANGE( '1', M, N, A, LDA, RWORK )

      if ( ANORM.GT.ZERO ) {
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, M, N ) ) ) / ANORM ) / ULP
      } else {
         RESULT( 1 ) = ZERO
      }

      // Compute B := V'*B*Q - D2*R

      sgemm('No transpose', 'No transpose', P, N, N, ONE, B, LDB, Q, LDQ, ZERO, WORK, LDB );

      sgemm('Transpose', 'No transpose', P, N, P, ONE, V, LDV, WORK, LDB, ZERO, B, LDB );

      for (I = 1; I <= L; I++) { // 100
         for (J = I; J <= L; J++) { // 90
            B( I, N-L+J ) = B( I, N-L+J ) - BETA( K+I )*R( K+I, K+J )
   90    CONTINUE
  100 CONTINUE

      // Compute norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP ) .

      RESID = SLANGE( '1', P, N, B, LDB, RWORK )
      if ( BNORM.GT.ZERO ) {
         RESULT( 2 ) = ( ( RESID / REAL( MAX( 1, P, N ) ) ) / BNORM ) / ULP
      } else {
         RESULT( 2 ) = ZERO
      }

      // Compute I - U'*U

      slaset('Full', M, M, ZERO, ONE, WORK, LDQ );
      ssyrk('Upper', 'Transpose', M, M, -ONE, U, LDU, ONE, WORK, LDU );

      // Compute norm( I - U'*U ) / ( M * ULP ) .

      RESID = SLANSY( '1', 'Upper', M, WORK, LDU, RWORK )
      RESULT( 3 ) = ( RESID / REAL( MAX( 1, M ) ) ) / ULP

      // Compute I - V'*V

      slaset('Full', P, P, ZERO, ONE, WORK, LDV );
      ssyrk('Upper', 'Transpose', P, P, -ONE, V, LDV, ONE, WORK, LDV );

      // Compute norm( I - V'*V ) / ( P * ULP ) .

      RESID = SLANSY( '1', 'Upper', P, WORK, LDV, RWORK )
      RESULT( 4 ) = ( RESID / REAL( MAX( 1, P ) ) ) / ULP

      // Compute I - Q'*Q

      slaset('Full', N, N, ZERO, ONE, WORK, LDQ );
      ssyrk('Upper', 'Transpose', N, N, -ONE, Q, LDQ, ONE, WORK, LDQ );

      // Compute norm( I - Q'*Q ) / ( N * ULP ) .

      RESID = SLANSY( '1', 'Upper', N, WORK, LDQ, RWORK )
      RESULT( 5 ) = ( RESID / REAL( MAX( 1, N ) ) ) / ULP

      // Check sorting

      scopy(N, ALPHA, 1, WORK, 1 );
      DO 110 I = K + 1, MIN( K+L, M )
         J = IWORK( I )
         if ( I.NE.J ) {
            TEMP = WORK( I )
            WORK( I ) = WORK( J )
            WORK( J ) = TEMP
         }
  110 CONTINUE

      RESULT( 6 ) = ZERO
      DO 120 I = K + 1, MIN( K+L, M ) - 1
         IF( WORK( I ).LT.WORK( I+1 ) ) RESULT( 6 ) = ULPINV
  120 CONTINUE

      RETURN

      // End of SGSVTS3

      }

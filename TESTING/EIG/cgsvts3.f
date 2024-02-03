      SUBROUTINE CGSVTS3( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V, LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               ALPHA( * ), BETA( * ), RESULT( 6 ), RWORK( * )
      COMPLEX            A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), Q( LDQ, * ), R( LDR, * ), U( LDU, * ), V( LDV, * ), WORK( LWORK )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, K, L;
      REAL               ANORM, BNORM, RESID, TEMP, ULP, ULPINV, UNFL
      // ..
      // .. External Functions ..
      REAL               CLANGE, CLANHE, SLAMCH
      // EXTERNAL CLANGE, CLANHE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CGGSVD3, CHERK, CLACPY, CLASET, SCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      ULP = SLAMCH( 'Precision' )
      ULPINV = ONE / ULP
      UNFL = SLAMCH( 'Safe minimum' )

      // Copy the matrix A to the array AF.

      clacpy('Full', M, N, A, LDA, AF, LDA );
      clacpy('Full', P, N, B, LDB, BF, LDB );

      ANORM = MAX( CLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
      BNORM = MAX( CLANGE( '1', P, N, B, LDB, RWORK ), UNFL )

      // Factorize the matrices A and B in the arrays AF and BF.

      cggsvd3('U', 'V', 'Q', M, N, P, K, L, AF, LDA, BF, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, LWORK, RWORK, IWORK, INFO );

      // Copy R

      DO 20 I = 1, MIN( K+L, M )
         for (J = I; J <= K + L; J++) { // 10
            R( I, J ) = AF( I, N-K-L+J )
         } // 10
      } // 20

      if ( M-K-L.LT.0 ) {
         for (I = M + 1; I <= K + L; I++) { // 40
            for (J = I; J <= K + L; J++) { // 30
               R( I, J ) = BF( I-K, N-K-L+J )
            } // 30
         } // 40
      }

      // Compute A:= U'*A*Q - D1*R

      cgemm('No transpose', 'No transpose', M, N, N, CONE, A, LDA, Q, LDQ, CZERO, WORK, LDA );

      cgemm('Conjugate transpose', 'No transpose', M, N, M, CONE, U, LDU, WORK, LDA, CZERO, A, LDA );

      for (I = 1; I <= K; I++) { // 60
         for (J = I; J <= K + L; J++) { // 50
            A( I, N-K-L+J ) = A( I, N-K-L+J ) - R( I, J )
         } // 50
      } // 60

      DO 80 I = K + 1, MIN( K+L, M )
         for (J = I; J <= K + L; J++) { // 70
            A( I, N-K-L+J ) = A( I, N-K-L+J ) - ALPHA( I )*R( I, J )
         } // 70
      } // 80

      // Compute norm( U'*A*Q - D1*R ) / ( MAX(1,M,N)*norm(A)*ULP ) .

      RESID = CLANGE( '1', M, N, A, LDA, RWORK )
      if ( ANORM.GT.ZERO ) {
         RESULT( 1 ) = ( ( RESID / REAL( MAX( 1, M, N ) ) ) / ANORM ) / ULP
      } else {
         RESULT( 1 ) = ZERO
      }

      // Compute B := V'*B*Q - D2*R

      cgemm('No transpose', 'No transpose', P, N, N, CONE, B, LDB, Q, LDQ, CZERO, WORK, LDB );

      cgemm('Conjugate transpose', 'No transpose', P, N, P, CONE, V, LDV, WORK, LDB, CZERO, B, LDB );

      for (I = 1; I <= L; I++) { // 100
         for (J = I; J <= L; J++) { // 90
            B( I, N-L+J ) = B( I, N-L+J ) - BETA( K+I )*R( K+I, K+J )
         } // 90
      } // 100

      // Compute norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP ) .

      RESID = CLANGE( '1', P, N, B, LDB, RWORK )
      if ( BNORM.GT.ZERO ) {
         RESULT( 2 ) = ( ( RESID / REAL( MAX( 1, P, N ) ) ) / BNORM ) / ULP
      } else {
         RESULT( 2 ) = ZERO
      }

      // Compute I - U'*U

      claset('Full', M, M, CZERO, CONE, WORK, LDQ );
      cherk('Upper', 'Conjugate transpose', M, M, -ONE, U, LDU, ONE, WORK, LDU );

      // Compute norm( I - U'*U ) / ( M * ULP ) .

      RESID = CLANHE( '1', 'Upper', M, WORK, LDU, RWORK )
      RESULT( 3 ) = ( RESID / REAL( MAX( 1, M ) ) ) / ULP

      // Compute I - V'*V

      claset('Full', P, P, CZERO, CONE, WORK, LDV );
      cherk('Upper', 'Conjugate transpose', P, P, -ONE, V, LDV, ONE, WORK, LDV );

      // Compute norm( I - V'*V ) / ( P * ULP ) .

      RESID = CLANHE( '1', 'Upper', P, WORK, LDV, RWORK )
      RESULT( 4 ) = ( RESID / REAL( MAX( 1, P ) ) ) / ULP

      // Compute I - Q'*Q

      claset('Full', N, N, CZERO, CONE, WORK, LDQ );
      cherk('Upper', 'Conjugate transpose', N, N, -ONE, Q, LDQ, ONE, WORK, LDQ );

      // Compute norm( I - Q'*Q ) / ( N * ULP ) .

      RESID = CLANHE( '1', 'Upper', N, WORK, LDQ, RWORK )
      RESULT( 5 ) = ( RESID / REAL( MAX( 1, N ) ) ) / ULP

      // Check sorting

      scopy(N, ALPHA, 1, RWORK, 1 );
      DO 110 I = K + 1, MIN( K+L, M )
         J = IWORK( I )
         if ( I.NE.J ) {
            TEMP = RWORK( I )
            RWORK( I ) = RWORK( J )
            RWORK( J ) = TEMP
         }
      } // 110

      RESULT( 6 ) = ZERO
      DO 120 I = K + 1, MIN( K+L, M ) - 1
         IF( RWORK( I ).LT.RWORK( I+1 ) ) RESULT( 6 ) = ULPINV
      } // 120

      RETURN

      // End of CGSVTS3

      }

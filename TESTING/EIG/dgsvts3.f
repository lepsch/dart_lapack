      SUBROUTINE DGSVTS3( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V, LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK, LWORK, RWORK, RESULT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             A( LDA, * ), AF( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), BF( LDB, * ), Q( LDQ, * ), R( LDR, * ), RESULT( 6 ), RWORK( * ), U( LDU, * ), V( LDV, * ), WORK( LWORK );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, INFO, J, K, L;
      double             ANORM, BNORM, RESID, TEMP, ULP, ULPINV, UNFL;
      // ..
      // .. External Functions ..
      double             DLAMCH, DLANGE, DLANSY;
      // EXTERNAL DLAMCH, DLANGE, DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMM, DGGSVD3, DLACPY, DLASET, DSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      ULP = DLAMCH( 'Precision' )
      ULPINV = ONE / ULP
      UNFL = DLAMCH( 'Safe minimum' )

      // Copy the matrix A to the array AF.

      dlacpy('Full', M, N, A, LDA, AF, LDA );
      dlacpy('Full', P, N, B, LDB, BF, LDB );

      ANORM = MAX( DLANGE( '1', M, N, A, LDA, RWORK ), UNFL )
      BNORM = MAX( DLANGE( '1', P, N, B, LDB, RWORK ), UNFL )

      // Factorize the matrices A and B in the arrays AF and BF.

      dggsvd3('U', 'V', 'Q', M, N, P, K, L, AF, LDA, BF, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, LWORK, IWORK, INFO );

      // Copy R

      DO 20 I = 1, MIN( K+L, M )
         for (J = I; J <= K + L; J++) { // 10
            R( I, J ) = AF( I, N-K-L+J )
         } // 10
      } // 20

      if ( M-K-L < 0 ) {
         for (I = M + 1; I <= K + L; I++) { // 40
            for (J = I; J <= K + L; J++) { // 30
               R( I, J ) = BF( I-K, N-K-L+J )
            } // 30
         } // 40
      }

      // Compute A:= U'*A*Q - D1*R

      dgemm('No transpose', 'No transpose', M, N, N, ONE, A, LDA, Q, LDQ, ZERO, WORK, LDA );

      dgemm('Transpose', 'No transpose', M, N, M, ONE, U, LDU, WORK, LDA, ZERO, A, LDA );

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

      RESID = DLANGE( '1', M, N, A, LDA, RWORK )

      if ( ANORM.GT.ZERO ) {
         RESULT( 1 ) = ( ( RESID / DBLE( MAX( 1, M, N ) ) ) / ANORM ) / ULP
      } else {
         RESULT( 1 ) = ZERO
      }

      // Compute B := V'*B*Q - D2*R

      dgemm('No transpose', 'No transpose', P, N, N, ONE, B, LDB, Q, LDQ, ZERO, WORK, LDB );

      dgemm('Transpose', 'No transpose', P, N, P, ONE, V, LDV, WORK, LDB, ZERO, B, LDB );

      for (I = 1; I <= L; I++) { // 100
         for (J = I; J <= L; J++) { // 90
            B( I, N-L+J ) = B( I, N-L+J ) - BETA( K+I )*R( K+I, K+J )
         } // 90
      } // 100

      // Compute norm( V'*B*Q - D2*R ) / ( MAX(P,N)*norm(B)*ULP ) .

      RESID = DLANGE( '1', P, N, B, LDB, RWORK )
      if ( BNORM.GT.ZERO ) {
         RESULT( 2 ) = ( ( RESID / DBLE( MAX( 1, P, N ) ) ) / BNORM ) / ULP
      } else {
         RESULT( 2 ) = ZERO
      }

      // Compute I - U'*U

      dlaset('Full', M, M, ZERO, ONE, WORK, LDQ );
      dsyrk('Upper', 'Transpose', M, M, -ONE, U, LDU, ONE, WORK, LDU );

      // Compute norm( I - U'*U ) / ( M * ULP ) .

      RESID = DLANSY( '1', 'Upper', M, WORK, LDU, RWORK )
      RESULT( 3 ) = ( RESID / DBLE( MAX( 1, M ) ) ) / ULP

      // Compute I - V'*V

      dlaset('Full', P, P, ZERO, ONE, WORK, LDV );
      dsyrk('Upper', 'Transpose', P, P, -ONE, V, LDV, ONE, WORK, LDV );

      // Compute norm( I - V'*V ) / ( P * ULP ) .

      RESID = DLANSY( '1', 'Upper', P, WORK, LDV, RWORK )
      RESULT( 4 ) = ( RESID / DBLE( MAX( 1, P ) ) ) / ULP

      // Compute I - Q'*Q

      dlaset('Full', N, N, ZERO, ONE, WORK, LDQ );
      dsyrk('Upper', 'Transpose', N, N, -ONE, Q, LDQ, ONE, WORK, LDQ );

      // Compute norm( I - Q'*Q ) / ( N * ULP ) .

      RESID = DLANSY( '1', 'Upper', N, WORK, LDQ, RWORK )
      RESULT( 5 ) = ( RESID / DBLE( MAX( 1, N ) ) ) / ULP

      // Check sorting

      dcopy(N, ALPHA, 1, WORK, 1 );
      DO 110 I = K + 1, MIN( K+L, M )
         J = IWORK( I )
         if ( I != J ) {
            TEMP = WORK( I )
            WORK( I ) = WORK( J )
            WORK( J ) = TEMP
         }
      } // 110

      RESULT( 6 ) = ZERO
      DO 120 I = K + 1, MIN( K+L, M ) - 1
         IF( WORK( I ) < WORK( I+1 ) ) RESULT( 6 ) = ULPINV
      } // 120

      RETURN

      // End of DGSVTS3

      }

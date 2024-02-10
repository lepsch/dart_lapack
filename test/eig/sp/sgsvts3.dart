      void sgsvts3(M, P, N, A, AF, LDA, B, BF, LDB, final Matrix<double> U, final int LDU, final Matrix<double> V, final int LDV, final Matrix<double> Q, final int LDQ, ALPHA, BETA, final Matrix<double> R, final int LDR, final Array<int> IWORK, final Array<double> WORK, final int LWORK, final Array<double> RWORK, RESULT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, N, P;
      int                IWORK( * );
      double               A( LDA, * ), AF( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), BF( LDB, * ), Q( LDQ, * ), R( LDR, * ), RESULT( 6 ), RWORK( * ), U( LDU, * ), V( LDV, * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      int                I, INFO, J, K, L;
      double               ANORM, BNORM, RESID, TEMP, ULP, ULPINV, UNFL;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLANGE, SLANSY;
      // EXTERNAL SLAMCH, SLANGE, SLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMM, SGGSVD3, SLACPY, SLASET, SSYRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL

      ULP = SLAMCH( 'Precision' );
      ULPINV = ONE / ULP;
      UNFL = SLAMCH( 'Safe minimum' );

      // Copy the matrix A to the array AF.

      slacpy('Full', M, N, A, LDA, AF, LDA );
      slacpy('Full', P, N, B, LDB, BF, LDB );

      ANORM = max( SLANGE( '1', M, N, A, LDA, RWORK ), UNFL );
      BNORM = max( SLANGE( '1', P, N, B, LDB, RWORK ), UNFL );

      // Factorize the matrices A and B in the arrays AF and BF.

      sggsvd3('U', 'V', 'Q', M, N, P, K, L, AF, LDA, BF, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, LWORK, IWORK, INFO );

      // Copy R

      for (I = 1; I <= min( K+L, M ); I++) { // 20
         for (J = I; J <= K + L; J++) { // 10
            R[I][J] = AF( I, N-K-L+J );
         } // 10
      } // 20

      if ( M-K-L < 0 ) {
         for (I = M + 1; I <= K + L; I++) { // 40
            for (J = I; J <= K + L; J++) { // 30
               R[I][J] = BF( I-K, N-K-L+J );
            } // 30
         } // 40
      }

      // Compute A:= U'*A*Q - D1*R

      sgemm('No transpose', 'No transpose', M, N, N, ONE, A, LDA, Q, LDQ, ZERO, WORK, LDA );

      sgemm('Transpose', 'No transpose', M, N, M, ONE, U, LDU, WORK, LDA, ZERO, A, LDA );

      for (I = 1; I <= K; I++) { // 60
         for (J = I; J <= K + L; J++) { // 50
            A[I][N-K-L+J] = A( I, N-K-L+J ) - R( I, J );
         } // 50
      } // 60

      for (I = K + 1; I <= min( K+L, M ); I++) { // 80
         for (J = I; J <= K + L; J++) { // 70
            A[I][N-K-L+J] = A( I, N-K-L+J ) - ALPHA( I )*R( I, J );
         } // 70
      } // 80

      // Compute norm( U'*A*Q - D1*R ) / ( max(1,M,N)*norm(A)*ULP ) .

      RESID = SLANGE( '1', M, N, A, LDA, RWORK );

      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / REAL( max( 1, M, N ) ) ) / ANORM ) / ULP;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute B := V'*B*Q - D2*R

      sgemm('No transpose', 'No transpose', P, N, N, ONE, B, LDB, Q, LDQ, ZERO, WORK, LDB );

      sgemm('Transpose', 'No transpose', P, N, P, ONE, V, LDV, WORK, LDB, ZERO, B, LDB );

      for (I = 1; I <= L; I++) { // 100
         for (J = I; J <= L; J++) { // 90
            B[I][N-L+J] = B( I, N-L+J ) - BETA( K+I )*R( K+I, K+J );
         } // 90
      } // 100

      // Compute norm( V'*B*Q - D2*R ) / ( max(P,N)*norm(B)*ULP ) .

      RESID = SLANGE( '1', P, N, B, LDB, RWORK );
      if ( BNORM > ZERO ) {
         RESULT[2] = ( ( RESID / REAL( max( 1, P, N ) ) ) / BNORM ) / ULP;
      } else {
         RESULT[2] = ZERO;
      }

      // Compute I - U'*U

      slaset('Full', M, M, ZERO, ONE, WORK, LDQ );
      ssyrk('Upper', 'Transpose', M, M, -ONE, U, LDU, ONE, WORK, LDU );

      // Compute norm( I - U'*U ) / ( M * ULP ) .

      RESID = SLANSY( '1', 'Upper', M, WORK, LDU, RWORK );
      RESULT[3] = ( RESID / REAL( max( 1, M ) ) ) / ULP;

      // Compute I - V'*V

      slaset('Full', P, P, ZERO, ONE, WORK, LDV );
      ssyrk('Upper', 'Transpose', P, P, -ONE, V, LDV, ONE, WORK, LDV );

      // Compute norm( I - V'*V ) / ( P * ULP ) .

      RESID = SLANSY( '1', 'Upper', P, WORK, LDV, RWORK );
      RESULT[4] = ( RESID / REAL( max( 1, P ) ) ) / ULP;

      // Compute I - Q'*Q

      slaset('Full', N, N, ZERO, ONE, WORK, LDQ );
      ssyrk('Upper', 'Transpose', N, N, -ONE, Q, LDQ, ONE, WORK, LDQ );

      // Compute norm( I - Q'*Q ) / ( N * ULP ) .

      RESID = SLANSY( '1', 'Upper', N, WORK, LDQ, RWORK );
      RESULT[5] = ( RESID / REAL( max( 1, N ) ) ) / ULP;

      // Check sorting

      scopy(N, ALPHA, 1, WORK, 1 );
      for (I = K + 1; I <= min( K+L, M ); I++) { // 110
         J = IWORK( I );
         if ( I != J ) {
            TEMP = WORK( I );
            WORK[I] = WORK( J );
            WORK[J] = TEMP;
         }
      } // 110

      RESULT[6] = ZERO;
      for (I = K + 1; I <= min( K+L, M ) - 1; I++) { // 120
         if[WORK( I ) < WORK( I+1 ) ) RESULT( 6] = ULPINV;
      } // 120

      }

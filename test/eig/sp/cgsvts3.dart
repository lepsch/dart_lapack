      void cgsvts3(final int M, final int P, final int N, final int A, final int AF, final int LDA, final int B, final int BF, final int LDB, final Matrix<double> U_, final int LDU, final Matrix<double> V_, final int LDV, final Matrix<double> Q_, final int LDQ, final int ALPHA, final int BETA, final Matrix<double> R_, final int LDR, final Array<int> IWORK_, final Array<double> WORK_, final int LWORK, final Array<double> RWORK_, final int RESULT,) {
  final U = U_.dim();
  final V = V_.dim();
  final Q = Q_.dim();
  final R = R_.dim();
  final IWORK = IWORK_.dim();
  final WORK = WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, N, P;
      int                IWORK( * );
      double               ALPHA( * ), BETA( * ), RESULT( 6 ), RWORK( * );
      Complex            A( LDA, * ), AF( LDA, * ), B( LDB, * ), BF( LDB, * ), Q( LDQ, * ), R( LDR, * ), U( LDU, * ), V( LDV, * ), WORK( LWORK );
      // ..

      double               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                I, INFO, J, K, L;
      double               ANORM, BNORM, RESID, TEMP, ULP, ULPINV, UNFL;
      // ..
      // .. External Functions ..
      //- REAL               CLANGE, CLANHE, SLAMCH;
      // EXTERNAL CLANGE, CLANHE, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CGGSVD3, CHERK, CLACPY, CLASET, SCOPY
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL

      ULP = SLAMCH( 'Precision' );
      ULPINV = ONE / ULP;
      UNFL = SLAMCH( 'Safe minimum' );

      // Copy the matrix A to the array AF.

      clacpy('Full', M, N, A, LDA, AF, LDA );
      clacpy('Full', P, N, B, LDB, BF, LDB );

      ANORM = max( CLANGE( '1', M, N, A, LDA, RWORK ), UNFL );
      BNORM = max( CLANGE( '1', P, N, B, LDB, RWORK ), UNFL );

      // Factorize the matrices A and B in the arrays AF and BF.

      cggsvd3('U', 'V', 'Q', M, N, P, K, L, AF, LDA, BF, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, LWORK, RWORK, IWORK, INFO );

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

      cgemm('No transpose', 'No transpose', M, N, N, CONE, A, LDA, Q, LDQ, CZERO, WORK, LDA );

      cgemm('Conjugate transpose', 'No transpose', M, N, M, CONE, U, LDU, WORK, LDA, CZERO, A, LDA );

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

      RESID = CLANGE( '1', M, N, A, LDA, RWORK );
      if ( ANORM > ZERO ) {
         RESULT[1] = ( ( RESID / REAL( max( 1, M, N ) ) ) / ANORM ) / ULP;
      } else {
         RESULT[1] = ZERO;
      }

      // Compute B := V'*B*Q - D2*R

      cgemm('No transpose', 'No transpose', P, N, N, CONE, B, LDB, Q, LDQ, CZERO, WORK, LDB );

      cgemm('Conjugate transpose', 'No transpose', P, N, P, CONE, V, LDV, WORK, LDB, CZERO, B, LDB );

      for (I = 1; I <= L; I++) { // 100
         for (J = I; J <= L; J++) { // 90
            B[I][N-L+J] = B( I, N-L+J ) - BETA( K+I )*R( K+I, K+J );
         } // 90
      } // 100

      // Compute norm( V'*B*Q - D2*R ) / ( max(P,N)*norm(B)*ULP ) .

      RESID = CLANGE( '1', P, N, B, LDB, RWORK );
      if ( BNORM > ZERO ) {
         RESULT[2] = ( ( RESID / REAL( max( 1, P, N ) ) ) / BNORM ) / ULP;
      } else {
         RESULT[2] = ZERO;
      }

      // Compute I - U'*U

      claset('Full', M, M, CZERO, CONE, WORK, LDQ );
      cherk('Upper', 'Conjugate transpose', M, M, -ONE, U, LDU, ONE, WORK, LDU );

      // Compute norm( I - U'*U ) / ( M * ULP ) .

      RESID = CLANHE( '1', 'Upper', M, WORK, LDU, RWORK );
      RESULT[3] = ( RESID / REAL( max( 1, M ) ) ) / ULP;

      // Compute I - V'*V

      claset('Full', P, P, CZERO, CONE, WORK, LDV );
      cherk('Upper', 'Conjugate transpose', P, P, -ONE, V, LDV, ONE, WORK, LDV );

      // Compute norm( I - V'*V ) / ( P * ULP ) .

      RESID = CLANHE( '1', 'Upper', P, WORK, LDV, RWORK );
      RESULT[4] = ( RESID / REAL( max( 1, P ) ) ) / ULP;

      // Compute I - Q'*Q

      claset('Full', N, N, CZERO, CONE, WORK, LDQ );
      cherk('Upper', 'Conjugate transpose', N, N, -ONE, Q, LDQ, ONE, WORK, LDQ );

      // Compute norm( I - Q'*Q ) / ( N * ULP ) .

      RESID = CLANHE( '1', 'Upper', N, WORK, LDQ, RWORK );
      RESULT[5] = ( RESID / REAL( max( 1, N ) ) ) / ULP;

      // Check sorting

      scopy(N, ALPHA, 1, RWORK, 1 );
      for (I = K + 1; I <= min( K+L, M ); I++) { // 110
         J = IWORK( I );
         if ( I != J ) {
            TEMP = RWORK( I );
            RWORK[I] = RWORK( J );
            RWORK[J] = TEMP;
         }
      } // 110

      RESULT[6] = ZERO;
      for (I = K + 1; I <= min( K+L, M ) - 1; I++) { // 120
         if[RWORK( I ) < RWORK( I+1 ) ) RESULT( 6] = ULPINV;
      } // 120

      }

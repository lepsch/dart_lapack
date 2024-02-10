      void clarcm(M, N, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, final Matrix<double> C, final int LDC, RWORK ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, LDC, M, N;
      double               A( LDA, * ), RWORK( * );
      Complex            B( LDB, * ), C( LDC, * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I, J, L;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, CMPLX, REAL
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM

      // Quick return if possible.

      if( ( M == 0 ) || ( N == 0 ) ) return;

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            RWORK[( J-1 )*M+I] = double( B( I, J ) );
         } // 10
      } // 20

      L = M*N + 1;
      sgemm('N', 'N', M, N, M, ONE, A, LDA, RWORK, M, ZERO, RWORK( L ), M );
      for (J = 1; J <= N; J++) { // 40
         for (I = 1; I <= M; I++) { // 30
            C[I][J] = RWORK( L+( J-1 )*M+I-1 );
         } // 30
      } // 40

      for (J = 1; J <= N; J++) { // 60
         for (I = 1; I <= M; I++) { // 50
            RWORK[( J-1 )*M+I] = AIMAG( B( I, J ) );
         } // 50
      } // 60
      sgemm('N', 'N', M, N, M, ONE, A, LDA, RWORK, M, ZERO, RWORK( L ), M );
      for (J = 1; J <= N; J++) { // 80
         for (I = 1; I <= M; I++) { // 70
            C[I][J] = CMPLX( double( C( I, J ) ), RWORK( L+( J-1 )*M+I-1 ) );
         } // 70
      } // 80

      }

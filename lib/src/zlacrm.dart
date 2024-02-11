      void zlacrm(final int M, final int N, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, final Matrix<double> C, final int LDC, final Array<double> RWORK,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                LDA, LDB, LDC, M, N;
      double             B( LDB, * ), RWORK( * );
      Complex         A( LDA, * ), C( LDC, * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                I, J, L;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, DIMAG
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM

      // Quick return if possible.

      if( ( M == 0 ) || ( N == 0 ) ) return;

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            RWORK[( J-1 )*M+I] = (A( I, J )).toDouble();
         } // 10
      } // 20

      L = M*N + 1;
      dgemm('N', 'N', M, N, N, ONE, RWORK, M, B, LDB, ZERO, RWORK( L ), M );
      for (J = 1; J <= N; J++) { // 40
         for (I = 1; I <= M; I++) { // 30
            C[I][J] = RWORK( L+( J-1 )*M+I-1 );
         } // 30
      } // 40

      for (J = 1; J <= N; J++) { // 60
         for (I = 1; I <= M; I++) { // 50
            RWORK[( J-1 )*M+I] = DIMAG( A( I, J ) );
         } // 50
      } // 60
      dgemm('N', 'N', M, N, N, ONE, RWORK, M, B, LDB, ZERO, RWORK( L ), M );
      for (J = 1; J <= N; J++) { // 80
         for (I = 1; I <= M; I++) { // 70
            C[I][J] = DCMPLX( (C( I, J )).toDouble(), RWORK( L+( J-1 )*M+I-1 ) );
         } // 70
      } // 80

      }

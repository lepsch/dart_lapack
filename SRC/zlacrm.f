      SUBROUTINE ZLACRM( M, N, A, LDA, B, LDB, C, LDC, RWORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LDC, M, N;
      // ..
      // .. Array Arguments ..
      double             B( LDB, * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), C( LDC, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D0, ZERO = 0.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, L;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, DIMAG
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM
      // ..
      // .. Executable Statements ..

      // Quick return if possible.

      IF( ( M == 0 ) || ( N == 0 ) ) RETURN

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            RWORK( ( J-1 )*M+I ) = DBLE( A( I, J ) )
         } // 10
      } // 20

      L = M*N + 1
      dgemm('N', 'N', M, N, N, ONE, RWORK, M, B, LDB, ZERO, RWORK( L ), M );
      for (J = 1; J <= N; J++) { // 40
         for (I = 1; I <= M; I++) { // 30
            C( I, J ) = RWORK( L+( J-1 )*M+I-1 )
         } // 30
      } // 40

      for (J = 1; J <= N; J++) { // 60
         for (I = 1; I <= M; I++) { // 50
            RWORK( ( J-1 )*M+I ) = DIMAG( A( I, J ) )
         } // 50
      } // 60
      dgemm('N', 'N', M, N, N, ONE, RWORK, M, B, LDB, ZERO, RWORK( L ), M );
      for (J = 1; J <= N; J++) { // 80
         for (I = 1; I <= M; I++) { // 70
            C( I, J ) = DCMPLX( DBLE( C( I, J ) ), RWORK( L+( J-1 )*M+I-1 ) )
         } // 70
      } // 80

      RETURN

      // End of ZLACRM

      }

      SUBROUTINE CLARCM( M, N, A, LDA, B, LDB, C, LDC, RWORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LDC, M, N;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), RWORK( * )
      COMPLEX            B( LDB, * ), C( LDC, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E0, ZERO = 0.0E0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, L;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, CMPLX, REAL
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM
      // ..
      // .. Executable Statements ..

      // Quick return if possible.

      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) RETURN

      for (J = 1; J <= N; J++) { // 20
         for (I = 1; I <= M; I++) { // 10
            RWORK( ( J-1 )*M+I ) = REAL( B( I, J ) )
         } // 10
      } // 20

      L = M*N + 1
      sgemm('N', 'N', M, N, M, ONE, A, LDA, RWORK, M, ZERO, RWORK( L ), M );
      for (J = 1; J <= N; J++) { // 40
         for (I = 1; I <= M; I++) { // 30
            C( I, J ) = RWORK( L+( J-1 )*M+I-1 )
         } // 30
      } // 40

      for (J = 1; J <= N; J++) { // 60
         for (I = 1; I <= M; I++) { // 50
            RWORK( ( J-1 )*M+I ) = AIMAG( B( I, J ) )
         } // 50
      } // 60
      sgemm('N', 'N', M, N, M, ONE, A, LDA, RWORK, M, ZERO, RWORK( L ), M );
      for (J = 1; J <= N; J++) { // 80
         for (I = 1; I <= M; I++) { // 70
            C( I, J ) = CMPLX( REAL( C( I, J ) ), RWORK( L+( J-1 )*M+I-1 ) )
         } // 70
      } // 80

      RETURN

      // End of CLARCM

      }

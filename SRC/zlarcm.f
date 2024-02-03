      SUBROUTINE ZLARCM( M, N, A, LDA, B, LDB, C, LDC, RWORK )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LDC, M, N;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), RWORK( * );
      COMPLEX*16         B( LDB, * ), C( LDC, * )
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

      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) RETURN

      DO 20 J = 1, N
         DO 10 I = 1, M
            RWORK( ( J-1 )*M+I ) = DBLE( B( I, J ) )
   10    CONTINUE
   20 CONTINUE

      L = M*N + 1
      dgemm('N', 'N', M, N, M, ONE, A, LDA, RWORK, M, ZERO, RWORK( L ), M );
      DO 40 J = 1, N
         DO 30 I = 1, M
            C( I, J ) = RWORK( L+( J-1 )*M+I-1 )
   30    CONTINUE
   40 CONTINUE

      DO 60 J = 1, N
         DO 50 I = 1, M
            RWORK( ( J-1 )*M+I ) = DIMAG( B( I, J ) )
   50    CONTINUE
   60 CONTINUE
      dgemm('N', 'N', M, N, M, ONE, A, LDA, RWORK, M, ZERO, RWORK( L ), M );
      DO 80 J = 1, N
         DO 70 I = 1, M
            C( I, J ) = DCMPLX( DBLE( C( I, J ) ), RWORK( L+( J-1 )*M+I-1 ) )
   70    CONTINUE
   80 CONTINUE

      RETURN

      // End of ZLARCM

      }

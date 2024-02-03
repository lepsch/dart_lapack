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

      DO 20 J = 1, N
         DO 10 I = 1, M
            RWORK( ( J-1 )*M+I ) = REAL( B( I, J ) )
   10    CONTINUE
   20 CONTINUE

      L = M*N + 1
      sgemm('N', 'N', M, N, M, ONE, A, LDA, RWORK, M, ZERO, RWORK( L ), M );
      DO 40 J = 1, N
         DO 30 I = 1, M
            C( I, J ) = RWORK( L+( J-1 )*M+I-1 )
   30    CONTINUE
   40 CONTINUE

      DO 60 J = 1, N
         DO 50 I = 1, M
            RWORK( ( J-1 )*M+I ) = AIMAG( B( I, J ) )
   50    CONTINUE
   60 CONTINUE
      sgemm('N', 'N', M, N, M, ONE, A, LDA, RWORK, M, ZERO, RWORK( L ), M );
      DO 80 J = 1, N
         DO 70 I = 1, M
            C( I, J ) = CMPLX( REAL( C( I, J ) ), RWORK( L+( J-1 )*M+I-1 ) )
   70    CONTINUE
   80 CONTINUE

      RETURN

      // End of CLARCM

      }

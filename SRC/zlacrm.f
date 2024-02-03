      SUBROUTINE ZLACRM( M, N, A, LDA, B, LDB, C, LDC, RWORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                LDA, LDB, LDC, M, N;
      // ..
      // .. Array Arguments ..
      double             B( LDB, * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), C( LDC, * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
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
*
      // Quick return if possible.
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) RETURN
*
      DO 20 J = 1, N
         DO 10 I = 1, M
            RWORK( ( J-1 )*M+I ) = DBLE( A( I, J ) )
   10    CONTINUE
   20 CONTINUE
*
      L = M*N + 1
      CALL DGEMM( 'N', 'N', M, N, N, ONE, RWORK, M, B, LDB, ZERO, RWORK( L ), M )
      DO 40 J = 1, N
         DO 30 I = 1, M
            C( I, J ) = RWORK( L+( J-1 )*M+I-1 )
   30    CONTINUE
   40 CONTINUE
*
      DO 60 J = 1, N
         DO 50 I = 1, M
            RWORK( ( J-1 )*M+I ) = DIMAG( A( I, J ) )
   50    CONTINUE
   60 CONTINUE
      CALL DGEMM( 'N', 'N', M, N, N, ONE, RWORK, M, B, LDB, ZERO, RWORK( L ), M )
      DO 80 J = 1, N
         DO 70 I = 1, M
            C( I, J ) = DCMPLX( DBLE( C( I, J ) ), RWORK( L+( J-1 )*M+I-1 ) )
   70    CONTINUE
   80 CONTINUE
*
      RETURN
*
      // End of ZLACRM
*
      END

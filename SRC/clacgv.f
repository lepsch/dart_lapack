      SUBROUTINE CLACGV( N, X, INCX )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INCX, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            X( * )
      // ..

* =====================================================================

      // .. Local Scalars ..
      int                I, IOFF;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG
      // ..
      // .. Executable Statements ..

      if ( INCX.EQ.1 ) {
         for (I = 1; I <= N; I++) { // 10
            X( I ) = CONJG( X( I ) )
   10    CONTINUE
      } else {
         IOFF = 1
         IF( INCX.LT.0 ) IOFF = 1 - ( N-1 )*INCX
         for (I = 1; I <= N; I++) { // 20
            X( IOFF ) = CONJG( X( IOFF ) )
            IOFF = IOFF + INCX
   20    CONTINUE
      }
      RETURN

      // End of CLACGV

      }

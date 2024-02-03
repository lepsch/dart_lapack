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
         DO 10 I = 1, N
            X( I ) = CONJG( X( I ) )
   10    CONTINUE
      } else {
         IOFF = 1
         IF( INCX.LT.0 ) IOFF = 1 - ( N-1 )*INCX
         DO 20 I = 1, N
            X( IOFF ) = CONJG( X( IOFF ) )
            IOFF = IOFF + INCX
   20    CONTINUE
      }
      RETURN

      // End of CLACGV

      }

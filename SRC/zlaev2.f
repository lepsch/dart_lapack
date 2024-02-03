      SUBROUTINE ZLAEV2( A, B, C, RT1, RT2, CS1, SN1 )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             CS1, RT1, RT2;
      COMPLEX*16         A, B, C, SN1
      // ..

* =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      double             ONE;
      const              ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      double             T;
      COMPLEX*16         W
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAEV2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG
      // ..
      // .. Executable Statements ..

      if ( ABS( B ).EQ.ZERO ) {
         W = ONE
      } else {
         W = DCONJG( B ) / ABS( B )
      }
      dlaev2(DBLE( A ), ABS( B ), DBLE( C ), RT1, RT2, CS1, T );
      SN1 = W*T
      RETURN

      // End of ZLAEV2

      }

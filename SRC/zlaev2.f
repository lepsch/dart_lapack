      void zlaev2(A, B, C, RT1, RT2, CS1, SN1 ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             CS1, RT1, RT2;
      Complex         A, B, C, SN1;
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0 ;
      double             ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      double             T;
      Complex         W;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAEV2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG
      // ..
      // .. Executable Statements ..

      if ( ABS( B ) == ZERO ) {
         W = ONE;
      } else {
         W = DCONJG( B ) / ABS( B );
      }
      dlaev2(DBLE( A ), ABS( B ), DBLE( C ), RT1, RT2, CS1, T );
      SN1 = W*T;
      return;

      // End of ZLAEV2

      }

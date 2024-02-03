      Complex     FUNCTION ZLADIV( X, Y );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      Complex         X, Y;
      // ..

// =====================================================================

      // .. Local Scalars ..
      double             ZI, ZR;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLADIV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, DIMAG
      // ..
      // .. Executable Statements ..

      dladiv(DBLE( X ), DIMAG( X ), DBLE( Y ), DIMAG( Y ), ZR, ZI );
      ZLADIV = DCMPLX( ZR, ZI );

      return;
      }

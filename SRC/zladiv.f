      COMPLEX*16     FUNCTION ZLADIV( X, Y );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX*16         X, Y;
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

      // End of ZLADIV

      }

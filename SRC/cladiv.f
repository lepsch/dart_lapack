      COMPLEX FUNCTION CLADIV( X, Y );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      COMPLEX            X, Y;
      // ..

// =====================================================================

      // .. Local Scalars ..
      REAL               ZI, ZR;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLADIV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC AIMAG, CMPLX, REAL
      // ..
      // .. Executable Statements ..

      sladiv(REAL( X ), AIMAG( X ), REAL( Y ), AIMAG( Y ), ZR, ZI );
      CLADIV = CMPLX( ZR, ZI );

      return;

      // End of CLADIV

      }

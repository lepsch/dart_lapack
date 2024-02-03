      void claipd(N, A, INDA, VINDA ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INDA, N, VINDA;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, IA, IXA;
      REAL               BIGNUM;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, REAL
      // ..
      // .. Executable Statements ..

      BIGNUM = SLAMCH( 'Epsilon' ) / SLAMCH( 'Safe minimum' );
      IA = 1;
      IXA = INDA;
      for (I = 1; I <= N; I++) { // 10
         A( IA ) = CMPLX( REAL( A( IA ) ), BIGNUM );
         IA = IA + IXA;
         IXA = IXA + VINDA;
      } // 10
      return;
      }

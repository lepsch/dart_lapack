      double second() {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

// =====================================================================

      // .. Local Scalars ..
      double               T1;
      double               TARRAY( 2 );
      // ..
      // .. Intrinsic Functions ..
      double               ETIME;
      // INTRINSIC ETIME

      T1 = ETIME( TARRAY );
      SECOND = TARRAY( 1 );
      }
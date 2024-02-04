      double second() {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

// =====================================================================

      // .. Local Scalars ..
      double               T1;
      // ..
      // .. Local Arrays ..
      double               TARRAY( 2 );
      // ..
      // .. Intrinsic Functions ..
      double               ETIME;
      // INTRINSIC ETIME
      // ..
      // .. Executable Statements ..

      T1 = ETIME( TARRAY );
      SECOND = TARRAY( 1 );
      return;
      }

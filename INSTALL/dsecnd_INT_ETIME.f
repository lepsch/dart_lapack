      double dsecnd() {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

// =====================================================================

      // .. Local Scalars ..
      REAL               T1;
      // ..
      // .. Local Arrays ..
      REAL               TARRAY( 2 );
      // ..
      // .. Intrinsic Functions ..
      REAL               ETIME;
      // INTRINSIC ETIME
      // ..
      // .. Executable Statements ..

      T1 = ETIME( TARRAY );
      DSECND = TARRAY( 1 );
      return;

      // End of DSECND

      }

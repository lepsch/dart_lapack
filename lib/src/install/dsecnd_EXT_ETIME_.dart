      double dsecnd() {

// -- LAPACK auxiliary routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
// =====================================================================

      // .. Local Scalars ..
      double               T1;
      // ..
      // .. Local Arrays ..
      double               TARRAY( 2 );
      // ..
      // .. External Functions ..
      //- REAL               ETIME_;
      // EXTERNAL ETIME_
      // ..
      // .. Executable Statements ..

      T1 = ETIME_( TARRAY );
      DSECND = TARRAY( 1 );
      return;
      }
      double second() {

// -- LAPACK auxiliary routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
// =====================================================================

      // .. Local Scalars ..
      double               T1;
      double               TARRAY( 2 );
      // ..
      // .. External Functions ..
      //- REAL               ETIME_;
      // EXTERNAL ETIME_

      T1 = ETIME_( TARRAY );
      SECOND = TARRAY( 1 );
      }

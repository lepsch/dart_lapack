      REAL FUNCTION SECOND( )

*  -- LAPACK auxiliary routine --
      // Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*=====================================================================

      // .. Local Scalars ..
      REAL               T1
      // ..
      // .. Local Arrays ..
      REAL               TARRAY( 2 )
      // ..
      // .. External Functions ..
      REAL               ETIME_
      // EXTERNAL ETIME_
      // ..
      // .. Executable Statements ..

      T1 = ETIME_( TARRAY )
      SECOND = TARRAY( 1 )
      RETURN

      // End of SECOND

      END

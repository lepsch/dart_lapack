      double           FUNCTION DSECND( );
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*
* =====================================================================
*
      // .. Local Scalars ..
      REAL               T1
      // ..
      // .. Local Arrays ..
      REAL               TARRAY( 2 )
      // ..
      // .. External Functions ..
      REAL               ETIME
      // EXTERNAL ETIME
      // ..
      // .. Executable Statements ..
*
      T1 = ETIME( TARRAY )
      DSECND = TARRAY( 1 )
      RETURN
*
      // End of DSECND
*
      END

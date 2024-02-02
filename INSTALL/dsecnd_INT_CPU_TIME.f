      DOUBLE PRECISION FUNCTION DSECND( )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
* =====================================================================
*
*     .. Local Scalars ..
*
      REAL T
*
* .. Intrinsic Functions ..
*
      INTRINSIC CPU_TIME
*
* .. Executable Statements .. *
*
      CALL CPU_TIME( T )
      DSECND = T
      RETURN
*
*     End of DSECND
*
      END
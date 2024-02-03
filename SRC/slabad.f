      SUBROUTINE SLABAD( SMALL, LARGE )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      REAL               LARGE, SMALL
      // ..
*
*  =====================================================================
*
      // .. Intrinsic Functions ..
      // INTRINSIC LOG10, SQRT
      // ..
      // .. Executable Statements ..
*
      // If it looks like we're on a Cray, take the square root of
      // SMALL and LARGE to avoid overflow and underflow problems.
*
       // IF( LOG10( LARGE ).GT.2000. ) THEN
          // SMALL = SQRT( SMALL )
          // LARGE = SQRT( LARGE )
       // END IF
*
      RETURN
*
      // End of SLABAD
*
      END

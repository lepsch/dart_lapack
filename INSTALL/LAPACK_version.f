      PROGRAM LAPACK_VERSION
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      int     MAJOR, MINOR, PATCH;
*     ..
*     .. External Subroutines ..
      // EXTERNAL ILAVER
*
      CALL ILAVER ( MAJOR, MINOR, PATCH )
      WRITE(*,*) "LAPACK ",MAJOR,".",MINOR,".",PATCH
*
      END

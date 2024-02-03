void main() {
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      int     MAJOR, MINOR, PATCH;
      // ..
      // .. External Subroutines ..
      // EXTERNAL ILAVER

      ilaver(MAJOR, MINOR, PATCH );
      WRITE(*,*) "LAPACK ",MAJOR,".",MINOR,".",PATCH

      }

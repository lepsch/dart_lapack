      SUBROUTINE DLASUM( TYPE, IOUNIT, IE, NRUN )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TYPE;
      int                IE, IOUNIT, NRUN;
      // ..

*  =====================================================================

      // .. Executable Statements ..

      if ( IE.GT.0 ) {
         WRITE( IOUNIT, FMT = 9999 )TYPE, ': ', IE, ' out of ', NRUN, ' tests failed to pass the threshold'
      } else {
         WRITE( IOUNIT, FMT = 9998 )'All tests for ', TYPE, ' passed the threshold ( ', NRUN, ' tests run)'
      }
 9999 FORMAT( 1X, A3, A2, I4, A8, I5, A35 )
 9998 FORMAT( / 1X, A14, A3, A24, I5, A11 )
      RETURN

      // End of DLASUM

      }

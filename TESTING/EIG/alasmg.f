      SUBROUTINE ALASMG( TYPE, NOUT, NFAIL, NRUN, NERRS )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TYPE;
      int                NFAIL, NOUT, NRUN, NERRS;
      // ..

* ======================================================================

      // .. Executable Statements ..

      if ( NFAIL.GT.0 ) {
         WRITE( NOUT, FMT = 9999 )TYPE, NFAIL, NRUN
      } else {
         WRITE( NOUT, FMT = 9998 )TYPE, NRUN
      }
      if ( NERRS.GT.0 ) {
         WRITE( NOUT, FMT = 9997 )NERRS
      }

 9999 FORMAT( 1X, A3, ': ', I6, ' out of ', I6,
     $      ' tests failed to pass the threshold' )
 9998 FORMAT( /1X, 'All tests for ', A3,
     $      ' routines passed the threshold ( ', I6, ' tests run)' )
 9997 FORMAT( 6X, I6, ' error messages recorded' )
      RETURN

      // End of ALASMG

      }

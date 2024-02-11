      void alasmg(final int TYPE, final int NOUT, final int NFAIL, final int NRUN, final int NERRS,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             TYPE;
      int                NFAIL, NOUT, NRUN, NERRS;
      // ..

// ======================================================================

      // .. Executable Statements ..

      if ( NFAIL > 0 ) {
         WRITE( NOUT, FMT = 9999 )TYPE, NFAIL, NRUN;
      } else {
         WRITE( NOUT, FMT = 9998 )TYPE, NRUN;
      }
      if ( NERRS > 0 ) {
         WRITE( NOUT, FMT = 9997 )NERRS;
      }

 9999 FORMAT(' ${.a3}: ${.i6} out of ${.i6} tests failed to pass the threshold' );
 9998 FORMAT('\n All tests for ${.a3} routines passed the threshold ( ${.i6} tests run)' );
 9997 FORMAT('${' ' * 6}${.i6} error messages recorded' );
      }
